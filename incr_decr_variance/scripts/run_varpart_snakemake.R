library(variancePartition)
library(BiocParallel)
library(SummarizedExperiment)
library(edgeR)
library(readr)
#library(doParallel)

#cl <- makeCluster(4)
#registerDoParallel(cl)
register(SerialParam())
registered()

# Input Paths -----------------------------------------------------
DGE.IN.PATH <- snakemake@input[[1]]
#DGE.IN.PATH <- "data/processed/limma_dgelist/dgelist.rds"

# Output Paths ----------------------------------------------------
VARPART.OUT.PATH <- snakemake@output[[1]]
META.OUT.PATH <- snakemake@output[[2]]
SUBJ.OUT.PATH <- snakemake@output[[3]]

# Read in data ----------------------------------------------------
dge <- readRDS(DGE.IN.PATH)

boot <- snakemake@params[["boot"]]
set.seed(boot)

# remove controls -------------------------------------------------
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

#new stuff added
n_subj_keep <- round(length(unique(dge$samples$Subject.ID)) * .8)
keep_subj <- sample(as.character(dge$samples$Subject.ID), size = n_subj_keep)
dge <- dge[, as.character(dge$samples$Subject.ID) %in% as.character(keep_subj)]

writeLines(as.character(keep_subj), SUBJ.OUT.PATH)

# make things factors ---------------------------------------------
meta <- dge$samples
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)

# Define formula --------------------------------------------------
form <-
  ~ (1|Subject.ID) + sex.numeric +  
  (1|RNA.isolation.Batch) + (1 | Lib_prep_batches) + (1|Year.Drawn) +
  Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent 


# Run variancePartition analysis ----------------------------------
design.for.voom <- model.matrix(~Lib_prep_batches + Age.months, data = meta)
v <- voom(dge, design.for.voom)
varPart <- fitExtractVarPartModel(v, form, meta)

saveRDS(varPart, VARPART.OUT.PATH)
write_tsv(meta, META.OUT.PATH)
