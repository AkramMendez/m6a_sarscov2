#!/usr/bin/env Rscript

library("Rsubread")
library("DESeq2")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("data.table")
library("ggtext")

args = commandArgs(trailingOnly=TRUE)

#inputdir <- "/crex/proj/nb_storage/private/BEA21P074_Roshan_Vaid/mappings/covid19_infected_cells/nodup_uniq_alns/all"
#annotation_file<-"/crex/proj/nb_project/private/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
inputdir <- "/home/akram/MondalLab/BEA21P074_Roshan_Vaid/mappings/covid19_infected_cells/nodup_uniq_alns/all/"
annotation_file<-"/home/akram/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
#metadata_path <- args[2]
threads <- 1



setwd(inputdir)

bamfiles<-list.files(inputdir,pattern = ".chroms.nodup.uniq.sorted.bam$")
bamfiles<-grep("(Input)",bamfiles,perl = T,ignore.case = T,value = T)

#bamfiles2<-grep("(Input_2)",bamfiles,perl = T,ignore.case = T,value = T)

metadata<- NULL
metadata$sample<-as.factor(basename(bamfiles))
metadata$condition<-as.factor(c(rep("infected",6),rep("noninfected",2)))
metadata$condition<-relevel(metadata$condition,"noninfected")
#metadata$type<-as.factor(c("B.1","B.1.351","B.1.1.7","Vero"))
metadata$type<-as.factor(c("B.1","B.1","B.1.351","B.1.351","B.1.1.7","B.1.1.7","Vero","Vero"))
metadata$type<-relevel(metadata$type,"Vero")
metadata$run<-as.factor(rep(c("r1","r2"),4))
#metadata$group<-as.factor(paste(metadata$type,metadata$run,sep="_"))
#metadata$condition<-as.factor(c(rep("infected",3),rep("noninfected",1)))
#metadata$type<-as.factor(c("B.1","B.1.351","B.1.1.7","Vero"))
#metadata$run<-as.factor(rep(c("r1"),4))
metadata<-as.data.frame(metadata)
rownames(metadata)<-gsub("_S\\d+.*$","",basename(bamfiles),perl=TRUE)



#metadata<-fread(metadata_path, colClasses = "factor")

featureCounts_scov2<-featureCounts(bamfiles, annot.ext=annotation_file, isGTFAnnotationFile = TRUE, nthreads = threads)

counts_scov2<-featureCounts_scov2$counts
colnames(counts_scov2)<-gsub("_S\\d+.*$","",colnames(counts_scov2),perl=TRUE)
fwrite(as.data.frame(counts_scov2),"counts_infectedCellsAll_cov2_nodup_uniq.tsv",col.names=T,sep="\t", row.names =T)

#counts_scov2<-fread("/home/akram/MondalLab/BEA21P074_Roshan_Vaid/mappings/covid19_infected_cells/nodup_uniq_alns/all/counts_infectedCellsAll_cov2_nodup_uniq.tsv",header = T,sep="\t")

# Check that sample names match those from the design matrix metadata
if(!all(colnames(counts_scov2) %in% rownames(metadata))){
  stop("colnames(counts_scov2) not matching rownames(coldata_info)")
}



dds<-DESeqDataSetFromMatrix(countData = counts_scov2, colData=metadata, design= ~ type)

#design(dds)<- ~ type

dds<-DESeq(dds,modelMatrixType = "standard")

resultsNames(dds)

#results_all<-results(dds)
#wu_vs_vero<-results(dds,contrast = c("type","B.1","Vero"))

#_________________________________#
normalized_counts<- counts(dds, normalized=TRUE)

normalized_counts<-normalized_counts %>% as.data.frame() %>%  mutate(gene=rownames(normalized_counts))
fwrite(normalized_counts,"~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/tables/normalizedCounts_diffexpAnalysis.tsv",sep="\t",col.names = T)
#results_all<-results(dds)

wu_vs_vero<-results(dds,contrast = c("type","B.1","Vero"),alpha = 0.01) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_vero")
uk_vs_vero<-results(dds,contrast = c("type","B.1.1.7","Vero")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_vero")
sa_vs_vero<-results(dds,contrast = c("type","B.1.351","Vero")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="sa_vs_vero")
# wu_vs_uk<-results(dds,contrast = c("type","B.1","B.1.1.7")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_uk")
# wu_vs_sa<-results(dds,contrast = c("type","B.1","B.1.351")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_sa")
# uk_vs_sa<-results(dds,contrast = c("type","B.1.1.7","B.1.351")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_sa")

#Shrink LFC estimates:
wu_vs_vero_lfcShrink<-lfcShrink(dds,coef = "type_B.1_vs_Vero",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_vero")
uk_vs_vero_lfcShrink<-lfcShrink(dds,coef = "type_B.1.1.7_vs_Vero",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_vero")
sa_vs_vero_lfcShrink<-lfcShrink(dds,coef = "type_B.1.351_vs_Vero",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="sa_vs_vero")
#wu_vs_uk_lfcShrhink<-lfcShrink(dds,coef = "type_WU_vs_UK",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_uk")
#wu_vs_sa_lfcShrhink<-lfcShrink(dds,coef = "type_WU_vs_SA",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_sa")
#uk_vs_sa_lfcShrhink<-lfcShrink(dds,coef = "type_UK_vs_SA",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_sa")



res_all<-rbind(wu_vs_vero,uk_vs_vero,sa_vs_vero)

res_all_lfcShrink<-rbind(wu_vs_vero_lfcShrink,uk_vs_vero_lfcShrink,sa_vs_vero_lfcShrink)

fwrite(res_all,"/home/akram/MondalLab/BEA21P074_Roshan_Vaid/differentialExpression/diffexp_infected_vs_noninfected_by_strain.tsv",col.names = T,sep="\t")
fwrite(res_all_lfcShrink,"/home/akram/MondalLab/BEA21P074_Roshan_Vaid/differentialExpression/diffexp_infected_vs_noninfected_by_strain_lfcShrink.tsv",col.names = T,sep="\t")
#res_all_lfcShrink<-fread("/home/akram/MondalLab/BEA21P074_Roshan_Vaid/differentialExpression/diffexp_infected_vs_noninfected_by_strain_lfcShrink.tsv",header = T,sep="\t")

res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene)) %>% arrange(dplyr::desc(contrast),dplyr::desc(log2FoldChange),padj) %>% dplyr::select(contrast,gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct() %>% 
  openxlsx::write.xlsx(.,"~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/tables/Table_differentialExpressionAnalysis.xls")

deg_tables <- list(
  "B.1_vs_NonInfectedVero"= res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & contrast=="wu_vs_vero" & abs(log2FoldChange) >1 & padj < 0.01 &! is.na(padj)) %>% arrange(dplyr::desc(log2FoldChange),padj) %>%
  dplyr::select(gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct() ,
  "B.1.1.7_vs_NonInfectedVero"= res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & contrast=="uk_vs_vero" & abs(log2FoldChange) >1 & padj < 0.01 &! is.na(padj)) %>% arrange(dplyr::desc(log2FoldChange),padj) %>%
    dplyr::select(gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct(),
  "B.1.351_vs_NonInfectedVero"= res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & contrast=="sa_vs_vero" & abs(log2FoldChange) >1 & padj < 0.01 &! is.na(padj)) %>% arrange(dplyr::desc(log2FoldChange),padj) %>%
    dplyr::select(gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct()
)

openxlsx::write.xlsx(deg_tables,"~/MondalLab/BEA21P074_Roshan_Vaid/plots/manuscript_figures/tables/DiffExprGenes_by_strain.xls")
