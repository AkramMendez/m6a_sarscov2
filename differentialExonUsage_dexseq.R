library(data.table)
library(tidyverse)
library(tximport)
library(rnaseqDTU)
library(GenomicFeatures)
library(DRIMSeq)
library(stageR)
library(DEXSeq)
library(cowplot)


# DEXseq data preparation -------------------------------------------------


indir<-"~/MondalLab/BEA21P074_Roshan_Vaid/dexseq_counts"
countFiles<-list.files(indir,pattern = ".txt$",full.names = T)

flattenedFile= list.files(indir, pattern = "ensembl.gff$",full.names = T)
basename(flattenedFile)
sampleTable = data.frame(
  row.names = c("SA_i1","SA_i2","UK_i1","UK_i2","Vero_i1","Vero_i2","WU_i1","WU_i2"),
  condition= c("SA","SA","UK","UK","Vero","Vero","WU","WU"),
  libType=rep("single-end",8)
)

sampleTable$condition<-factor(sampleTable$condition)
sampleTable$condition<-relevel(sampleTable$condition,"Vero")

dxd_dataset <- DEXSeqDataSetFromHTSeq(countfiles = countFiles,
                                      sampleData = sampleTable,
                                      design = ~sample + exon + condition:exon,
                                      flattenedfile = flattenedFile)


# DEXseq Normalization ----------------------------------------------------


dxd_dataset = estimateSizeFactors(dxd_dataset,locfunc=shorth)
dxd_dataset = estimateDispersions(dxd_dataset)
#par(mar=c(1,1,1,1))
plotDispEsts(dxd_dataset)



# Test for differential exon usage: ---------------------------------------


dxd_dataset = testForDEU(dxd_dataset)

dxd_dataset = estimateExonFoldChanges(dxd_dataset, fitExpToVar = "condition")


dex_results <- DEXSeqResults(dxd_dataset)
save(dex_results,file = "~/MondalLab/BEA21P074_Roshan_Vaid/dexseq_counts/dex_results_by_strain.rda")
load(file = "~/MondalLab/BEA21P074_Roshan_Vaid/dexseq_counts/dex_results_by_strain.rda")
mcols(dex_results)$description

q.cutoff<-0.1


fwrite(dex_results %>% as.data.frame() %>% filter(padj < q.cutoff),"~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/tables/diffExonUsage_DEXseq_results_cov2_padj0.05.tsv", sep="\t")



#Impute empty gene names:
map_ensembl_to_entrezid<-getBM(filters ="ensembl_gene_id" , attributes = c("ensembl_gene_id","entrezgene_accession","hgnc_symbol","uniprot_gn_symbol","external_gene_name"),values =dex_filt$groupID ,mart = mart)
dex_filt<- left_join(dex_filt,map_ensembl_to_entrezid %>% dplyr::select(ensembl_gene_id,entrezgene_accession),by=c("groupID"="ensembl_gene_id")) %>% mutate(imputed_gene_name=if_else(entrezgene_accession=="",gene,entrezgene_accession)) %>% distinct()
#annotated_gene_symbols_df<- dex_filt_gene_symbol %>% mutate(gene=if_else(hgnc_symbol=="",uniprot_gn_symbol,hgnc_symbol)) %>% mutate(gene=if_else(gene=="",ensembl_gene_id,gene,ensembl_gene_id)) 

deuGenes_with_peaks_annotation<- left_join(dex_filt,peaks_table,by=c("imputed_gene_name"="peakAnn_gene_name")) %>% 
  filter(!is.na(name_macs2) &! is.na(condition) & short_annotation !="intergenic") %>% rename(condition_peak=condition)


# Number of peaks per condition in DEU genes: ----------------------------

deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak) %>% 
  dplyr::count() %>% ungroup()

#___________________Barplot Number of m6A peaks in DEU genes _______________________#
deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,groupID,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","wu","uk","sa")),n,fill=condition_peak)) + 
  geom_bar(stat = "identity") +
  xlab("") + ylab("Number of peaks") + labs(title="Number of peaks in DEU genes",fill="") + theme_light() + theme(aspect.ratio = 1.5) 



# Peaks per gene by condition (DEU genes): --------------------------------


countsPeaksPerGeneDEU <- deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,groupID,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak,groupID) %>% 
  dplyr::count()

countsPeaksPerGeneDEU

# Summary Peaks Per DEU Gene ----------------------------------------------
# Boxplot: Peaks per gene by condition (DEU genes) ------------------------


deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,groupID,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak,groupID) %>% 
  dplyr::count() %>% ungroup() %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","wu","uk","sa")),n,color=condition_peak)) + 
  geom_point(alpha=0.7,position = position_jitter(height = 0.05), size=0.75) +
  geom_boxplot(alpha=0.5,notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.5, colour="black", width=0.5)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + 
  theme(aspect.ratio = 1, panel.background = element_blank(), panel.grid = element_blank(),legend.position = "none") +
  ylab("Number of peaks per gene") +
  xlab("") + 
  labs(fill="",title="Peaks per gene",subtitle ="DEU genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","wu"="B.1","uk"="B.1.1.7","sa"="B.1.351")) +
  scale_y_continuous(breaks = seq(0,15), labels=seq(0,15)) + ylim(0,10) ->plotPeaksPerGene_DEU
plotPeaksPerGene_DEU

ggsave("~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/compiled_figures/peaksPerGene_DEUGenes.pdf",plot = plotPeaksPerGene_DEU,device = "pdf", width = 6,height = 8,dpi = 300)

  
#scale_y_continuous(breaks=c(0:10),labels=as.character(c(0:10)), limits = c(0.001,10), expand = c(0.001,0.5)) +


#Number of peaks per DEU gene:
  peaksPerDEUGene<- deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID,name_macs2) %>% 
    distinct() %>% 
    group_by(condition_peak,groupID) %>% 
    dplyr::count() %>% ungroup()


# Summary at GENE level per condition: ------------------------------------


  numDEUGenes_with_m6a<- deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID) %>% 
    distinct() %>% group_by(condition_peak) %>% dplyr::count() %>% as.data.frame()
  
  numDEUs_with_m6a_vero<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="vero")]
  numDEUs_with_m6a_wu<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="wu")]
  numDEUs_with_m6a_uk<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="uk")]
  numDEUs_with_m6a_sa<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="sa")]
  
  numPeaks_DEUgenes<-deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID,name_macs2) %>% 
    distinct() %>% 
    group_by(condition_peak) %>% 
    dplyr::count()

  #summary_peaks_in_DEU<- numPeaks_DEUgenes$n / c(numDEUs_with_m6a_sa,numDEUs_with_m6a_uk,numDEUs_with_m6a_vero,numDEUs_with_m6a_wu)
  

# Table Peaks per gene summary, DEU genes ---------------------------------

  summaryPeaksPerDEUGene <- peaksPerDEUGene %>% group_by(condition_peak) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
  summaryPeaksPerDEUGene$condition_peak <-factor(summaryPeaksPerDEUGene$condition_peak,levels = c("vero","wu","uk","sa"))
  summaryPeaksPerDEUGene %>% arrange(condition_peak) %>% kbl() %>% kable_classic(full_width=F)


# Number of Peaks by Annotation, DEU genes --------------------------------

  #Where are peaks from DEU genes located (exon, intron, etc):
  numPeaksByAnnotationDEUGenes<- deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID,name_macs2, short_annotation) %>% 
    distinct() %>% 
    group_by(condition_peak, short_annotation) %>% 
    dplyr::count() %>% ungroup()
  
  numPeaksByAnnotationDEUGenes %>% 
    ggplot(aes(fct_relevel(condition_peak,levels=c("vero","wu","uk","sa")),n,fill=short_annotation)) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    xlab("Condition") + ylab("Number of peaks per annotaion") + labs(title = "Number of peaks by annotation", subtitle = "DEU genes")
  
  #Number of DEU genes with direct intersection with m6A peaks:
  dex_filt_all_exonic_m6a_intersections %>% dplyr::select(condition_peaks,groupID) %>% group_by(condition_peaks) %>% dplyr::count()
  
  #Number of peaks in DEU genes with direct intersection of m6A:
  dex_filt_all_exonic_m6a_intersections %>% dplyr::select(condition_peaks,groupID,peak_name_peaks) %>% group_by(condition_peaks,groupID) %>% dplyr::count()
  
#DEU genes having m6A modification directly at the predicted DEU exons
#Number of peaks per gene:
  peaksPerDEUIntersectingExonm6A<-dex_filt_all_exonic_m6a_intersections %>% 
    dplyr::select(condition_peaks,groupID,peak_name_peaks) %>%
    distinct() %>%
    group_by(condition_peaks,groupID) %>% 
    dplyr::count() %>% ungroup()
  
  
# Table PeaksPerGene, DEUs with intersecting m6A at exon ------------------

  peaksPerDEUIntersectingExonm6A$condition_peaks <- factor(peaksPerDEUIntersectingExonm6A$condition_peaks,levels = c("Vero","WU","UK","SA"))
  peaksPerDEUIntersectingExonm6A %>% group_by(condition_peaks) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
  
#Number of peaks per exon:
peaksPerExonDEUIntersectingExonm6A<-dex_filt_all_exonic_m6a_intersections %>% 
  dplyr::select(condition_peaks,groupID,featureID,peak_name_peaks) %>%
    distinct() %>%
    group_by(condition_peaks,groupID,featureID) %>% 
    dplyr::count() %>% ungroup()

# Table Peaks Per Exon m6A intersecting Exon -----------------------------------

  peaksPerExonDEUIntersectingExonm6A$condition_peaks <- factor(peaksPerExonDEUIntersectingExonm6A$condition_peaks,levels = c("Vero","WU","UK","SA"))
  peaksPerExonDEUIntersectingExonm6A %>% group_by(condition_peaks) %>% 
    summarise(mean_num_peaks_per_gene=mean(n), median(n)) %>%
    arrange(condition_peaks) %>% kbl() %>% kable_classic(full_width=F)
  

  

# Boxplot peaks per gene in DEU genes intersecting m6A at exons -----------

  peaksPerGene_intersectingm6AtExon<- dex_filt_all_exonic_m6a_intersections %>%
  dplyr::select(condition_peaks,groupID,peak_name_peaks) %>% 
  distinct() %>% 
  group_by(condition_peaks,groupID) %>% 
  dplyr::count() %>%
  ungroup() 
  
  # Table Peaks Per Gene, DEU genes with Intersecting m6A at Exon -----------
  
  peaksPerGene_intersectingm6AtExon <- peaksPerGene_intersectingm6AtExon %>% group_by(condition_peaks) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
  peaksPerGene_intersectingm6AtExon$condition_peaks <-factor(summaryPeaksPerNonDEUGene$condition_peak,levels = c("Vero","WU","UK","SA"))
  peaksPerGene_intersectingm6AtExon %>% arrange(condition_peaks) %>% kbl() %>% kable_classic(full_width=F)
  
  
  peaksPerGene_intersectingm6AtExon %>% 
  ggplot(aes(fct_relevel(condition_peaks,levels=c("Vero","WU","UK","SA")),n,color=condition_peaks)) + 
  geom_point(alpha=0.1,position = position_jitterdodge()) +
  geom_boxplot(alpha=0.5,notch = F, notchwidth = 0.75, varwidth = T,outlier.size = 0.5, colour="black")  +
  stat_summary(fun = "mean", color="red")+
  theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank()) +
  ylab("Number of peaks per gene") + xlab("") + labs(fill="", title="Peaks per gene",subtitle = "m6A overlapping with DEU exons") -> plot_numPeaksGeneDEU_intersecting_m6A
  plot_numPeaksGeneDEU_intersecting_m6A

  # Boxplot peaks per Exon in DEU genes intersecting m6A at exons -----------
  
  dex_filt_all_exonic_m6a_intersections %>%
    dplyr::select(condition_peaks,groupID,featureID,peak_name_peaks) %>% 
    distinct() %>% 
    group_by(condition_peaks,groupID,featureID) %>% 
    dplyr::count() %>%
    ungroup() %>% 
    ggplot(aes(fct_relevel(condition_peaks,levels=c("vero","WU","UK","SA")),n,color=condition_peaks)) + 
    geom_point(alpha=0.1,position = position_jitterdodge()) +
    geom_boxplot(alpha=0.5,notch = F, notchwidth = 0.75, varwidth = T,outlier.size = 0.5, colour="black")  +
    stat_summary(fun.y = "mean", color="red") +
    theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank()) +
    ylab("Number of peaks per exon") + xlab("") + labs(fill="", title="Peaks per exon",subtitle = "m6A overlapping with DEU exons")





# Peaks Per Gene , Non-DEU genes ------------------------------------------


peaksPerGene_nonDEU<- allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(!(deg_gene %in% unique(dex_filt$imputed_gene_name))) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% 
  ungroup() 

# Table, summary peaks per gene, Non-DEU genes ----------------------------

summaryPeaksPerNonDEUGene <- peaksPerGene_nonDEU %>% group_by(condition_peak) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
summaryPeaksPerNonDEUGene$condition_peak <-factor(summaryPeaksPerNonDEUGene$condition_peak,levels = c("vero","WU","UK","SA"))
summaryPeaksPerNonDEUGene %>% arrange(condition_peak) %>% kbl() %>% kable_classic(full_width=F)


# Boxplot Peaks Per Gene, Non-DEU genes -----------------------------------

peaksPerGene_nonDEU %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")),n,color=condition_peak)) + 
  #stat_boxplot(geom = "errorbar",width=0.2, size=0.5, colour="black") +
  geom_point(alpha=0.1,position = position_jitter(), size=0.1) +
  geom_boxplot(alpha=0.7,notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.1,outlier.alpha = 0.5, colour="black", width=0.5)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank(), legend.position = "none") +
  ylab("Number of peaks per gene") + xlab("") + labs(color="",fill="", title = "Number of peaks per gene",subtitle = "Non-DEU genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","WU"="B.1","UK"="B.1.1.7","SA"="B.1.351")) +
  scale_y_discrete(limits=c(0,15), labels=seq(0,15)) + ylim(0,15)->plotPeaksPerGene_nonDEU
plotPeaksPerGene_nonDEU


#Counting peaks in randomly selected genes from the total with a group size equal to the number of DEU genes:
numberDEUgenes_vero<- peaksPerDEUGene %>% filter(condition_peak=="vero") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()
numberDEUgenes_wu<- peaksPerDEUGene %>% filter(condition_peak=="wu") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()
numberDEUgenes_uk<- peaksPerDEUGene %>% filter(condition_peak=="uk") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()
numberDEUgenes_sa<- peaksPerDEUGene %>% filter(condition_peak=="sa") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()

randomGenes_vero<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
                      dplyr::select(deg_gene) %>% 
                      distinct() %>% unlist(),size=numberDEUgenes_vero, replace = F)

randomGenes_wu<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="WU")  %>% 
                            dplyr::select(deg_gene) %>% 
                            distinct() %>% unlist(),size=numberDEUgenes_wu, replace = F)
randomGenes_uk<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="UK")  %>% 
                            dplyr::select(deg_gene) %>% 
                            distinct() %>% unlist(),size=numberDEUgenes_uk, replace = F)
randomGenes_sa<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="SA")  %>% 
                            dplyr::select(deg_gene) %>% 
                            distinct() %>% unlist(),size=numberDEUgenes_sa, replace = F)


# Peaks  per Gene, Random Genes -------------------------------------------


randomPeakCounts_vero <- allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_vero)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts_wu<-allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="WU")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_wu)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts_uk<-allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="UK")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_uk)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts_sa<-allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="SA")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_sa)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts <-rbind(randomPeakCounts_vero,randomPeakCounts_wu,randomPeakCounts_uk,randomPeakCounts_sa)

# Table, summary peaks per gene, Random Genes -----------------------------

summaryrandomPeakCounts <- randomPeakCounts %>% group_by(condition_peak) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
summaryrandomPeakCounts$condition_peak <-factor(summaryrandomPeakCounts$condition_peak,levels = c("vero","WU","UK","SA"))
summaryrandomPeakCounts %>% arrange(condition_peak) %>% kbl() %>% kable_classic(full_width=F)


# Boxplot Peaks Per Gene, Random Genes ------------------------------------


randomPeakCounts %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")),n,color=condition_peak)) + 
  #stat_boxplot(geom = "errorbar",width=0.2, size=0.5) +
  geom_point(alpha=0.7,position = position_jitter(height = 0.05), size=0.75) +
  geom_boxplot(notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.5, colour="black", width=0.5, alpha=0.7)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank(), legend.position = "none") +
  ylab("Number of peaks per gene") + xlab("") + labs(fill="", title = "Number of peaks per gene",subtitle = "Random genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","WU"="B.1","UK"="B.1.1.7","SA"="B.1.351")) +
  scale_y_continuous(breaks=seq(0,15), labels=seq(0,15)) + ylim(0,15) -> plotPeaksPerGene_random
plotPeaksPerGene_random





# Peaks per gene, DEU vs Random -------------------------------------------
randomSampCountsDEU <- randomPeakCounts %>% filter(condition_peak=="vero") %>% mutate(type="Random")
deu_and_random<- rbind(countsPeaksPerGeneDEU %>% ungroup() %>% rename(deg_gene=groupID) %>% filter(condition_peak=="vero") %>% mutate(type="DEU"),randomSampCountsDEU)

deu_and_random %>% ggplot(aes(fct_relevel(type,levels=c("vero","wu","uk","sa")),n)) + 
  geom_point(alpha=0.7,position = position_jitter(height = 0.05), size=0.15, color="steelblue") +
  geom_boxplot(alpha=0.5,notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.5, colour="black", width=0.5)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + 
  theme(aspect.ratio = 1.25, panel.background = element_blank(), panel.grid = element_blank(),legend.position = "none") +
  ylab("Number of peaks per gene") +
  xlab("") + 
  labs(fill="",title="Peaks per gene",subtitle ="DEU genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","wu"="B.1","uk"="B.1.1.7","sa"="B.1.351")) +
  scale_y_continuous(breaks = seq(0,15), labels=seq(0,15)) + ylim(0,10) ->plotPeaksPerGene_DEU_random

plotPeaksPerGene_DEU_random




# Null distribution -------------------------------------------------------

peaksPerGeneVero <- allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

point_estimate_mean <- countsPeaksPerGeneDEU %>% filter(condition_peak=="vero") %>% 
  ungroup() %>% specify(response = n) %>% calculate(stat="mean")

point_estimate_median <- countsPeaksPerGeneDEU %>% filter(condition_peak=="vero") %>% 
  ungroup() %>% specify(response = n) %>% calculate(stat="median")

#Run once:
  getRandomGenes <- function(df,n, size){
    set.seed(123)
    sampdf <- NULL
    
    for(i in 1:n){
      
      randomGenes<-sample( df %>% filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
                                  dplyr::select(deg_gene) %>% 
                                  distinct() %>% unlist(),size=size, replace = F)
      
      df1 <- df %>% 
        filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
        dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
        distinct() %>% 
        filter(deg_gene %in% unique(randomGenes)) %>% 
        group_by(condition_peak,deg_gene) %>% 
        dplyr::count() %>% ungroup() %>% 
        mutate(replicate=i)
      
      sampdf <-rbind(sampdf,df1)
      df1<-NULL
    }
    
    return(sampdf)
  }

  
  randomPeakCounts_vero_samp1000 <- getRandomGenes(df=allgenes_and_peaks_retained_lost_gained,n=1000,size=numberDEUgenes_vero)
  
  #save(randomPeakCounts_vero_samp1000,file="~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/compiled_figures/tables/randomSampling1000rep_DEUgenestest_vero.Rda")
  load(file="~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/compiled_figures/tables/randomSampling1000rep_DEUgenestest_vero.Rda") 
  
  
  mean_random_dist <- randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(mean=mean(n)) %>% summarise(mean(mean)) %>% unlist()
  
  pdf("~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/compiled_figures/compilation/suppFig4F_new_numPeaksDEUgenes_vs_random.pdf", width = 4, height =3)

# Permutation test DEU genes
  randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(mean=mean(n)) %>% 
    ggplot(aes(mean)) + geom_histogram(bins = 100) + 
    geom_vline(xintercept = unlist(point_estimate_mean), color="red") +
    geom_vline(xintercept = mean_random_dist, color="blue") +
    theme_bw(base_size = 11) + 
    xlab("Mean number of peaks") +
    ylab("Count") +
    annotate(geom = "text", x=2.2,y=75,label="t-test \n p < 2.2e-16") +
    scale_y_continuous(breaks=seq(0,80,by=10)) +
    scale_x_continuous(breaks=seq(1.5,3,by=0.25))
  
  dev.off()
  
  randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(median=median(n)) %>% 
    ggplot(aes(median)) + geom_histogram(bins = 100) + geom_vline(xintercept = unlist(point_estimate_median), color="red") +
    theme_bw(base_size = 11) + xlab("Median number of peaks")
  
  t.test(randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(mean=mean(n)) %>% dplyr::select(mean) %>% unlist(), mu=point_estimate_mean$stat)
