library(Seurat)
library(edgeR)
library(tidyverse)
source("scripts/pseudobulk_pooling_functions.R")


dataset.all <- readRDS(snakemake@input[["seurat_obj"]])
DGELISTS_OUT_PATH <- "data/output/NICA_batch1to4_new_timepoints_cleaned.rds"

#Parameters for pseudobulk pooling
LIBSIZE_FILTER <- 0
MIN_CELLS_PER_POOL <- 5
MIN_SAMPLES_PER_CELLTYPE <- 15

#Cell annotation column
CELL_ANNOTATION_COLUMN <- "manual.cid"

#read in data
dataset.all.filt<-dataset.all[,which(dataset.all$gender %in% c("M","F"))]

seurat_obj <- dataset.all.filt

meta <- seurat_obj@meta.data
meta$barcodes<-colnames(dataset.all.filt)

rna <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

sample_cols <- c("matched.individual", "matched.timepoint.age", "gender", "age", "batch")
sample_cols <- c(sample_cols, CELL_ANNOTATION_COLUMN)


print("sample_cols that aren't in metadata")
sample_cols[!sample_cols %in% colnames(meta)]

samples <- paste(meta$matched.individual, meta$matched.timepoint.age, meta$batch, sep = "_")

pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = meta[[CELL_ANNOTATION_COLUMN]], 
                                     meta = meta,
                                     samples = samples, 
                                     min_cells_per_pool = MIN_CELLS_PER_POOL, 
                                     min_samples_per_celltype = MIN_SAMPLES_PER_CELLTYPE,
                                     barcode_col_name = "barcodes", 
                                     sample_level_meta_cols = sample_cols,
                                     pooling_function = "sum",
                                     output_type = "DGEList")

pseudobulk_list <- lapply(pseudobulk_list, function(dge){
  dge[, dge$samples$lib.size > LIBSIZE_FILTER]
  
})

n_samples_per_celltype <- sapply(pseudobulk_list, ncol)
print("n_samples_per_celltype before final filtering")
print(n_samples_per_celltype)

pseudobulk_list <- pseudobulk_list[n_samples_per_celltype > MIN_SAMPLES_PER_CELLTYPE]

saveRDS(pseudobulk_list, DGELISTS_OUT_PATH)
