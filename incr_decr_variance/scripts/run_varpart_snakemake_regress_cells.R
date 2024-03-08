library(variancePartition)
library(BiocParallel)
library(SummarizedExperiment)
library(readr)
#library(doParallel)

#cl <- makeCluster(4)
#registerDoParallel(cl)
register(SerialParam())
registered()

# Input Paths -----------------------------------------------------
vst.IN.PATH <- snakemake@input[[1]]
#vst.IN.PATH <- "data/processed/limma_vstlist/vstlist.rds"

# Output Paths ----------------------------------------------------
VARPART.OUT.PATH <- snakemake@output[[1]]
META.OUT.PATH <- snakemake@output[[2]]
SUBJ.OUT.PATH <- snakemake@output[[3]]

# Read in data ----------------------------------------------------
vst <- readRDS(vst.IN.PATH)

boot <- snakemake@params[["boot"]]
set.seed(boot)

# remove controls -------------------------------------------------
vst <- vst[, vst$Sample.type == "case" & !is.na(vst$RIN)]

#new stuff added
n_subj_keep <- round(length(unique(vst$Subject.ID)) * .8)
keep_subj <- sample(as.character(vst$Subject.ID), size = n_subj_keep)
vst <- vst[, as.character(vst$Subject.ID) %in% as.character(keep_subj)]

writeLines(as.character(keep_subj), SUBJ.OUT.PATH)

# make things factors ---------------------------------------------
meta <- as.data.frame(colData(vst))
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
varPart <- fitExtractVarPartModel(assays(vst)[[1]], form, meta)

saveRDS(varPart, VARPART.OUT.PATH)
write_tsv(meta, META.OUT.PATH)
