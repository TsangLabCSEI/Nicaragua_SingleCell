library(Seurat)
library(edgeR)
library(tidyverse)
library(lspline)
library(variancePartition)

#adult cohort pseudobulk variance partition loop

all_tts_vpar<-list()
for (a in c("all_adult_pbulk_list_normalized")){
  print(a)
  DGE.IN.PATH <- paste0(a,".rds")

  # Output Paths ----------------------------------------------------
  FIT.OUT.PATH <- "limma_fit_no_subj_elspline_8knot_fit_fsex.rds"
  SPLINE.DESIGN.OUT.PATH <- "limma_fit_no_subj_elspline_8knot_spline_design_fsex.rds"
  TT.OUT.PATH <- "limma_fit_no_subj_elspline_8knot_toptab_anova_fsex.csv"

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
    fitmm = dream( vobjDream, form, meta)
    tt <- topTable(fitmm, "Age", number = Inf)
    tt$gene<-rownames(tt)

    idx<-as.numeric(str_split_i(str_split_i(c,"_",-1),"\\.",1))
    if(grepl("bcells",c)){
      tt$cell.type<-"Bcells"
      all_tts_bcells_vpar[[idx]]<-tt
    } else if (grepl("nkcells",c)){
      tt$cell.type<-"NKcells"
      all_tts_nkcells_vpar[[idx]]<-tt
    } else if (grepl("mye",c)){
      tt$cell.type<-"Monocytes"
      all_tts_myecells_vpar[[idx]]<-tt
    } else if (grepl("cd4",c)){
      tt$cell.type<-"Tcells_CD4"
      all_tts_tcellscd4_vpar[[idx]]<-tt
    } else if (grepl("cd8",c)){
      tt$cell.type<-"Tcells_CD8"
      all_tts_tcellscd8_vpar[[idx]]<-tt   
    } 
  }
  all_tts_vpar[["Bcells"]]<-all_tts_bcells_vpar
  all_tts_vpar[["NKcells"]]<-all_tts_nkcells_vpar
  all_tts_vpar[["Monocytes"]]<-all_tts_myecells_vpar
  all_tts_vpar[["Tcells_CD4"]]<-all_tts_tcellscd4_vpar
  all_tts_vpar[["Tcells_CD8"]]<-all_tts_tcellscd8_vpar
}
saveRDS(all_tts_vpar,"age_subject_variance_VESlist_adult_old_aging_cohort.RDS")
