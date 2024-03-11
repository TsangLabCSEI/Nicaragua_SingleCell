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

#new stuff added
n_subj_keep <- round(length(unique(vst$matched.individual)) * .8)
keep_subj <- sample(as.character(vst$matched.individual), size = n_subj_keep)
vst <- vst[, as.character(vst$matched.individual) %in% as.character(keep_subj)]

writeLines(as.character(keep_subj), SUBJ.OUT.PATH)

# make things factors ---------------------------------------------
meta <- as.data.frame(colData(vst))
meta$Subject.ID <- factor(meta$matched.individual)
meta$sex.numeric <- as.numeric(factor(meta$gender))
meta$Age.months <- factor(meta$matched.timepoint.age * 12)
meta$RNA.isolation.Batch <- factor(meta$batch)

# Define formula --------------------------------------------------
form <- ~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + Age.months

# Run variancePartition analysis ----------------------------------
varPart <- fitExtractVarPartModel(assays(vst)[[1]], form, meta)

saveRDS(varPart, VARPART.OUT.PATH)
write_tsv(meta, META.OUT.PATH)
