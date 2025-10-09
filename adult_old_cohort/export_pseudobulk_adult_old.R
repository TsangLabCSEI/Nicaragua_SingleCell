library(Seurat)
library(zellkonverter)
source("pseudobulk_pooling_functions.R")


#adult cohort raw count data loading example myecells
sce <- readH5AD("all_pbmcs_rna_harmony.h5ad")
sce.seurat <- as.Seurat(sce[,sce$Cluster_names=="Myeloid cells"], counts = "X", data = "X") # sce$Cluster_names=="TRAV1-2- CD8+ T cells"
sce.raw <- readH5AD("pbmc_gex_raw_with_var_obs.h5ad", use_hdf5 = TRUE)
sce.mye <- sce.raw[,sce.raw$Cluster_names=="Myeloid cells"]
assay(sce.mye, "X") <- as.matrix(assay(sce.mye, "X"))
sce.seurat.mye <- as.Seurat(sce.mye, counts = "X", data = "X")
sce.seurat.mye$Donor_id<-sce.seurat$Donor_id
sce.seurat.mye$Batch<-sce.seurat$Batch
sce.seurat.mye$Age<-sce.seurat$Age
sce.seurat.mye$Sex<-sce.seurat$Sex


#adult cohort pseudobulk loop example myecells
donors<-readRDS("age_subject_variance_donors_adult_old_aging_cohort.RDS")
for (batch in seq(1,20)){
    print(batch)
    d<-donors[[batch]]
    DGELISTS_OUT_PATH <- paste0("adults_myecells_",batch,".rds")

    #Parameters for pseudobulk pooling
    LIBSIZE_FILTER <- 0
    MIN_CELLS_PER_POOL <- 5
    MIN_SAMPLES_PER_CELLTYPE <- 15
    #Cell annotation column
    CELL_ANNOTATION_COLUMN <- "orig.ident"
    #read in data
    seurat_obj <- sce.seurat.mye[,sce.seurat.mye$Donor_id %in% d]
    meta <- seurat_obj@meta.data
    meta$barcodes<-colnames(seurat_obj)
    rna <- GetAssayData(seurat_obj, assay = "originalexp", slot = "counts")
    sample_cols <- c("Donor_id", "Age", "Sex", "Batch")
    sample_cols <- c(sample_cols, CELL_ANNOTATION_COLUMN)
    samples <- paste(meta$Donor_id, meta$Age, meta$Batch, sep = "_")
    
    pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = meta[[CELL_ANNOTATION_COLUMN]], 
                                         meta = meta,
                                         samples = samples, 
                                         min_cells_per_pool = MIN_CELLS_PER_POOL, 
                                         min_samples_per_celltype = MIN_SAMPLES_PER_CELLTYPE,
                                         barcode_col_name = "barcodes", 
                                         sample_level_meta_cols = sample_cols,
                                         #cell_level_meta_cols = cell_cols, 
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
}
