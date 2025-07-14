library(limma)
library(tidyverse)
library(edgeR)
library(variancePartition)
library(fgsea)

dir.create("outputs")

pbulk <- readRDS("NICA_batch2to4_new_timepoints_cleaned_normalized.rds")

fit_list <- lapply(pbulk, function(dge){
  meta <- dge$samples %>%
    select(matched.individual, matched.timepoint.age) %>%
    rename(age = matched.timepoint.age, subject = matched.individual)
  
  design <- model.matrix(~age + subject, meta)
  
  v <- voom(dge, design = design)
  
  fit <- lmFit(v, design = design)
  
  fit <- eBayes(fit)
})

saveRDS(fit_list, "outputs/limma_fit_list.rds")

toptab_list <- lapply(fit_list, function(fit){
  topTable(fit, coef = "age", number = Inf) %>%
    rownames_to_column("gene")
})

lapply(toptab_list, function(dat){
  dat %>% filter(adj.P.Val < .05) %>%
    pull(t) %>%
    sign() %>%
    table()
})

toptab_list[[1]] %>%
  filter(AveExpr > 10)

toptab_list[[2]] %>%
  filter(AveExpr > 10)

geneset.list <- readRDS("combined_gene_sets.RDS")

fgsea_list <- lapply(toptab_list, function(toptab){
  fc <- toptab$logFC
  names(fc) <- toptab$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = fc,
                    minSize=15)
})

saveRDS(fgsea_list, "outputs/NICA_batch2to4_new_timepoints_limma_fgsea_simplified_age_results.rds")
