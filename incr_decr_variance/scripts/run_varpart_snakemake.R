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
dir.create(dirname(VARPART.OUT.PATH), recursive = TRUE)
dir.create(dirname(META.OUT.PATH), recursive = TRUE)
dir.create(dirname(SUBJ.OUT.PATH), recursive = TRUE)


# Read in data ----------------------------------------------------
dge <- readRDS(DGE.IN.PATH)

boot <- snakemake@params[["boot"]]
set.seed(boot)

#new stuff added
n_subj_keep <- 10 #round(length(unique(dge$samples$matched.individual)) * .8)
keep_subj <- sample(as.character(dge$samples$matched.individual), size = n_subj_keep)
dge <- dge[, as.character(dge$samples$matched.individual) %in% as.character(keep_subj)]

writeLines(as.character(keep_subj), SUBJ.OUT.PATH)

# make things factors ---------------------------------------------
meta <- dge$samples
meta$Subject.ID <- factor(meta$matched.individual)
meta$sex.numeric <- as.numeric(factor(meta$gender))
meta$Age.months <- as.numeric(meta$matched.timepoint.age * 12)
meta$RNA.isolation.Batch <- factor(meta$batch)

# Define formula --------------------------------------------------
form <- ~ (1|Subject.ID) + sex.numeric +  (1|RNA.isolation.Batch) + Age.months


# Run variancePartition analysis ----------------------------------
v <- voomWithDreamWeights( dge, form, meta, suppressWarnings=TRUE)
varPart <- fitExtractVarPartModel(v, form, meta)

saveRDS(varPart, VARPART.OUT.PATH)
write_tsv(meta, META.OUT.PATH)
