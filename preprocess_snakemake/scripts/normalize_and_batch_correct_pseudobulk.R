library(reshape2)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(purrr)
library(edgeR)
library(Seurat)
library(edgeR)
library(tidyverse)


var.threshold <- 1

#initialize list to store normalized pseudobulk dta
pseudobulk_list<-readRDS(snakemake@input[["pbulk_obj"]])
pseudobulk_list_normalized <- list()

DGELISTS_OUT_PATH <- "data/output/NICA_batch2to4_new_timepoints_cleaned_normalized.rds"


#iterate over raw pseudobulk pooled data and run normalization and batch correction
for (i in names(pseudobulk_list)[which(!names(pseudobulk_list) %in% c("Platelet","Platelets","unlabeled","doublets"))]) {
  cat(i,"-\n")
  eset <- pseudobulk_list[[i]]
  geneExpr <- eset$counts
  sample.info <- eset$samples
  sample.info$sample.id <- paste0(sample.info$matched.individual,".",sample.info$matched.timepoint.age)
  
  # filter out lowly expressed genes
  cat("Number of genes with no counts across samples:\n")
  print(table(rowSums(eset$counts==0)==nrow(sample.info)))
  keep.exprs <- filterByExpr(eset,group = sample.info$gender, min.count = 2,min.prop=0.5)
  cat("Number of genes to keep by expr:",sum(keep.exprs),"\n")
  eset <- eset[keep.exprs,,keep.lib.sizes = F]
 
  # stage 1: normalize and correct
  eset <- calcNormFactors(eset)
  sample.info$gender <- ifelse(is.na(sample.info$gender),"Control",as.character(sample.info$gender))

  # a model that accommodate the control samples
  tmp.model <- model.matrix(~ 0 + gender,sample.info)
  vobj <- voom(eset ,design = tmp.model,plot = F)
  geneExpr <- vobj$E # normalized
  
  # PCA using baseline samples before correction
  baseline.samples <- sample.info
  baseline.geneExpr <- geneExpr[,rownames(baseline.samples)]
  baseline.geneExpr.pca <- prcomp(t(baseline.geneExpr),center = T,scale. = T)
  baseline.samples <- cbind(baseline.samples,baseline.geneExpr.pca$x[rownames(baseline.samples),1:2])
  print(ggplot(baseline.samples,aes(PC1,PC2)) + geom_point(aes(color=paste0(batch),shape=matched.individual == "Control"),size=2,alpha=0.8) + 
          ggtitle(paste0("Before correction: ",i)))
  cat("Before correction:\n")
  print(car::Anova(lm(PC1 ~ n_barcodes + batch,baseline.samples)))
  print(car::Anova(lm(PC2 ~ n_barcodes + batch,baseline.samples)))
  
  # remove batch effect
  normalized.expr.batch.effect <- removeBatchEffect(geneExpr, batch = eset$samples$batch, 
                                                    covariates = eset$samples$n_barcodes,
                                                    design = tmp.model)
  geneExpr <- normalized.expr.batch.effect
  eset$normalizedExpr <- geneExpr
  pseudobulk_list_normalized[[i]] <- eset

  # PCA using baseline samples after correction
  baseline.samples <- sample.info
  baseline.geneExpr <- geneExpr[,rownames(baseline.samples)]
  baseline.geneExpr.pca <- prcomp(t(baseline.geneExpr),center = T,scale. = T)
  baseline.samples <- cbind(baseline.samples,baseline.geneExpr.pca$x[rownames(baseline.samples),1:2])
  print(ggplot(baseline.samples,aes(PC1,PC2)) + geom_point(aes(color=paste0(batch),shape=matched.individual == "Control"),size=2,alpha=0.8) + 
          ggtitle(paste0("After correction: ",i)))
  cat("After correction:\n")
  print(car::Anova(lm(PC1 ~ n_barcodes + batch,baseline.samples)))
  print(car::Anova(lm(PC2 ~ n_barcodes + batch,baseline.samples)))
}

saveRDS(pseudobulk_list_normalized,DGELISTS_OUT_PATH)
