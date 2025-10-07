library(Seurat)
library(edgeR)
library(tidyverse)
library(lspline)
library(variancePartition)
library(stringr)

#adult cohort pseudobulk variance partition loop

all_tts_vpar<-list()
for (a in c(str_split_i(snakemake@input[["pbulk_obj"]],"\\.",1))){
  print(a)

  # Output Paths ----------------------------------------------------
  DGE.IN.PATH <- paste0(a,".rds")

  all_tts_bcells_vpar<-list() 
  all_tts_nkcells_vpar<-list()
  all_tts_myecells_vpar<-list()
  all_tts_tcellscd4_vpar<-list()
  all_tts_tcellscd8_vpar<-list()

  # Read in data ----------------------------------------------------
  dge_a <- readRDS(DGE.IN.PATH)
  for (c in names(dge_a)) { #
    print(c)
    dge <- dge_a[[c]]

    # remove controls -------------------------------------------------
    dge <- dge[, !is.na(dge$samples$Age)]

    # make things factors ---------------------------------------------
    meta <- dge$samples
    meta$Donor_id <- factor(meta$Donor_id)
    meta$Age <- factor(meta$Age)
    meta$Age <- as.numeric(meta$Age)*12
    meta$Batch <- factor(meta$Batch)

    print("fit simplified age model")
    form <- ~ (1|Donor_id) + Age + (1 | Batch)
    vobjDream = voomWithDreamWeights( dge, form, meta, suppressWarnings=TRUE)
    vp = fitExtractVarPartModel( vobjDream, form, meta)

    idx<-as.numeric(str_split_i(str_split_i(c,"_",-1),"\\.",1))
    if(grepl("bcells",c)){
      vp$cell.type<-"Bcells"
      all_tts_bcells_vpar[[idx]]<-vp
    } else if (grepl("nkcells",c)){
      vp$cell.type<-"NKcells"
      all_tts_nkcells_vpar[[idx]]<-vp
    } else if (grepl("mye",c)){
      vp$cell.type<-"Monocytes"
      all_tts_myecells_vpar[[idx]]<-vp
    } else if (grepl("cd4",c)){
      vp$cell.type<-"Tcells_CD4"
      all_tts_tcellscd4_vpar[[idx]]<-vp
    } else if (grepl("cd8",c)){
      vp$cell.type<-"Tcells_CD8"
      all_tts_tcellscd8_vpar[[idx]]<-vp   
    } 
  }
  all_tts_vpar[["Bcells"]]<-all_tts_bcells_vpar
  all_tts_vpar[["NKcells"]]<-all_tts_nkcells_vpar
  all_tts_vpar[["Monocytes"]]<-all_tts_myecells_vpar
  all_tts_vpar[["Tcells_CD4"]]<-all_tts_tcellscd4_vpar
  all_tts_vpar[["Tcells_CD8"]]<-all_tts_tcellscd8_vpar
}
saveRDS(all_tts_vpar,"data/output/age_subject_variance_VESlist_adult_old_aging_cohort.RDS")
