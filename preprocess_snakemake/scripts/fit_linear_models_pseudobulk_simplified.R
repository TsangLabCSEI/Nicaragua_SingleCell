library(limma)
library(tidyverse)
library(edgeR)
library(variancePartition)
library(fgsea)
library(dplyr)

pbulk <- readRDS(snakemake@input[["pbulk_obj"]]) #"data/output/NICA_batch2to4_new_timepoints_cleaned_normalized.rds") #snakemake@input[["pbulk_obj"]])

fit_list <- lapply(pbulk, function(dge){
  meta <- dge$samples %>%
    select(matched.individual, matched.timepoint.age) %>%
    rename(age = matched.timepoint.age, subject = matched.individual)
  
  design <- model.matrix(~age + subject, meta)
  
  v <- voom(dge, design = design)
  
  fit <- lmFit(v, design = design)
  
  fit <- eBayes(fit)
})

saveRDS(fit_list, "data/output/simplified_age_model_limma_fit_list.rds")

toptab_list <- lapply(fit_list, function(fit){
  topTable(fit, coef = "age", number = Inf) %>%
    rownames_to_column("gene")
})

geneset.list <- readRDS(snakemake@input[["gene_sets"]]) #"data/input/gene_sets/combined_gene_sets.RDS") #snakemake@input[["gene_sets"]])

fgsea_list <- lapply(toptab_list, function(toptab){
  fc <- toptab$logFC
  names(fc) <- toptab$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = fc,
                    minSize=15)
})

saveRDS(fgsea_list, "data/output/NICA_batch2to4_new_timepoints_limma_fgsea_simplified_age_results.rds")
