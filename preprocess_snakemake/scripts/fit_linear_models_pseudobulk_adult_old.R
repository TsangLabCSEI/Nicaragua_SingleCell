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

  all_tts_bcellsnaive_vpar<-list()
  all_tts_bcellsmem_vpar<-list() 
  all_tts_nkcellscd16hi_vpar<-list()
  all_tts_nkcellscd56hi_vpar<-list()
  all_tts_myecells_vpar<-list()
  all_tts_mononc_vpar<-list()
  all_tts_tcellscd4naive_vpar<-list()   
  all_tts_tcellscd4treg_vpar<-list()
  all_tts_tcellscd8naive_vpar<-list()
  all_tts_tcellscd8gzmb_vpar<-list()
  all_tts_tcellscd8temra_vpar<-list()
  all_tts_tcellscd8em_vpar<-list()
  all_tts_pdcs_vpar<-list()
  all_tts_cdcs_vpar<-list()


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
    if(grepl("bnaive",c)){
	vp$cell.type<-"B_Naive"
	all_tts_bcellsnaive_vpar[[idx]]<-vp
    } else if (grepl("bmem",c)){
	vp$cell.type<-"B_Mem"
	all_tts_bcellsmem_vpar[[idx]]<-vp
    } else if (grepl("nk16",c)){
	vp$cell.type<-"NK_CD16hi"
	all_tts_nkcellscd16hi_vpar[[idx]]<-vp
    } else if (grepl("nk56",c)){
	vp$cell.type<-"NK_CD56hiCD16lo"
	all_tts_nkcellscd56hi_vpar[[idx]]<-vp
    } else if (grepl("mye",c)){
	vp$cell.type<-"Mono_Classical"
	all_tts_myecells_vpar[[idx]]<-vp
    } else if (grepl("mono",c)){
	vp$cell.type<-"Mono_NonClassical"
	all_tts_mononc_vpar[[idx]]<-vp
    } else if (grepl("cd4naive",c)){
	vp$cell.type<-"CD4_Naive"
	all_tts_tcellscd4naive_vpar[[idx]]<-vp
    } else if (grepl("cd4treg",c)){
	vp$cell.type<-"CD4_Treg"
	all_tts_tcellscd4treg_vpar[[idx]]<-vp
    } else if (grepl("cd8naive",c)){
	vp$cell.type<-"CD8_Naive"
	all_tts_tcellscd8naive_vpar[[idx]]<-vp   
    } else if (grepl("cd8gzmb",c)){
	vp$cell.type<-"CD8_EM_GZMB+"
	all_tts_tcellscd8gzmb_vpar[[idx]]<-vp   
    } else if (grepl("cd8temra",c)){
	vp$cell.type<-"CD8_TEMRA"
	all_tts_tcellscd8temra_vpar[[idx]]<-vp   
    } else if (grepl("cd8em",c)){
	vp$cell.type<-"CD8_EM"
	all_tts_tcellscd8em_vpar[[idx]]<-vp   
    } else if (grepl("cdc",c)){
	vp$cell.type<-"cDC"
	all_tts_cdcs_vpar[[idx]]<-vp   
    } else if (grepl("pdc",c)){
	vp$cell.type<-"pDC"
	all_tts_pdcs_vpar[[idx]]<-vp   
    }
  }
  all_tts_vpar[["B_Naive"]]<-all_tts_bcellsnaive_vpar
  all_tts_vpar[["B_Mem"]]<-all_tts_bcellsmem_vpar
  all_tts_vpar[["NK_CD16hi"]]<-all_tts_nkcellscd16hi_vpar
  all_tts_vpar[["NK_CD56hiCD16lo"]]<-all_tts_nkcellscd56hi_vpar
  all_tts_vpar[["Mono_Classical"]]<-all_tts_myecells_vpar
  all_tts_vpar[["Mono_NonClassical"]]<-all_tts_mononc_vpar
  all_tts_vpar[["CD4_Naive"]]<-all_tts_tcellscd4naive_vpar
  all_tts_vpar[["CD4_Treg"]]<-all_tts_tcellscd4treg_vpar
  all_tts_vpar[["CD8_Naive"]]<-all_tts_tcellscd8naive_vpar
  all_tts_vpar[["CD8_EM_GZMB"]]<-all_tts_tcellscd8gzmb_vpar
  all_tts_vpar[["CD8_TEMRA"]]<-all_tts_tcellscd8temra_vpar
  all_tts_vpar[["CD8_EM"]]<-all_tts_tcellscd8em_vpar
  all_tts_vpar[["cDC"]]<-all_tts_cdcs_vpar
  all_tts_vpar[["pDC"]]<-all_tts_pdcs_vpar
}
saveRDS(all_tts_vpar,"data/output/age_subject_variance_VESlist_adult_old_aging_cohort.RDS")
