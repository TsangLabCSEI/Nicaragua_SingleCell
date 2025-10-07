library(reshape2)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(purrr)
library(edgeR)
library(Seurat)
library(edgeR)
library(tidyverse)


#adult cohort pseudobulk normalization loop

var.threshold <- 1
pseudobulk_list_normalized <- list()

for (i in list.files(snakemake@input[["pbulk_folder"]],pattern="adults.*.rds")){
  cat(i,"-\n")
  path<-paste0(snakemake@input[["pbulk_folder"]],"/",i)
  eset <- readRDS(path)$SeuratProject
  geneExpr <- eset$counts
  sample.info <- eset$samples
  sample.info$sample.id <- paste0(sample.info$Donor_id,".",sample.info$Age)
  
  # filter out lowly expressed genes
  cat("Number of genes with no counts across samples:\n")
  print(table(rowSums(eset$counts==0)==nrow(sample.info)))
  keep.exprs <- filterByExpr(eset,group = sample.info$Sex, min.count = 2,min.prop=0.5)
  cat("Number of genes to keep by expr:",sum(keep.exprs),"\n")
  eset <- eset[keep.exprs,,keep.lib.sizes = F]
 
  # stage 1: normalize and correct
  eset <- calcNormFactors(eset)
  sample.info$Sex <- ifelse(is.na(sample.info$Sex),"Control",as.character(sample.info$Sex))
  # a model that accommodate the control samples
  tmp.model <- model.matrix(~ 0 + Batch,sample.info)
  vobj <- voom(eset ,design = tmp.model,plot = F)
  #title("",i)
  geneExpr <- vobj$E # normalized
  
  # PCA using baseline samples before correction
  baseline.samples <- sample.info
  baseline.geneExpr <- geneExpr[,rownames(baseline.samples)]
  baseline.geneExpr.pca <- prcomp(t(baseline.geneExpr),center = T,scale. = T)
  baseline.samples <- cbind(baseline.samples,baseline.geneExpr.pca$x[rownames(baseline.samples),1:2])
  print(ggplot(baseline.samples,aes(PC1,PC2)) + geom_point(aes(color=paste0(Batch),shape=Donor_id == "Control"),size=2,alpha=0.8) + 
          ggtitle(paste0("Before correction: ",i)))
  cat("Before correction:\n")
  print(car::Anova(lm(PC1 ~ n_barcodes + Batch,baseline.samples)))
  print(car::Anova(lm(PC2 ~ n_barcodes + Batch,baseline.samples)))
  
  # remove batch effect
  normalized.expr.batch.effect <- removeBatchEffect(geneExpr, batch = eset$samples$Batch, 
                                                    covariates = eset$samples$n_barcodes,
                                                    design = tmp.model)
  geneExpr <- normalized.expr.batch.effect
  eset$normalizedExpr <- geneExpr
  pseudobulk_list_normalized[[i]] <- eset
}
saveRDS(pseudobulk_list_normalized,"data/output/all_adult_pbulk_list_normalized.rds")
