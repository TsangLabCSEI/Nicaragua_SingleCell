library(Seurat)
library(edgeR)
library(tidyverse)
library(lspline)
library(variancePartition)
library(stringr)

for (a in c(str_split_i(snakemake@input[["pbulk_obj"]],"\\.",1))){
  print(a)
  DGE.IN.PATH <- paste0(a,".rds")
  
  # Output Paths ----------------------------------------------------
  VAR.OUT.PATH <- paste0(a,"_fsex_vpar.rds")
  SPLINE.OUT.PATH <- paste0(a,"_fsex_spline.rds")

  all_tts_vpar<-list()
  all_tts<-list()
  all_tts_fitmm<-list()
  all_tts_spline<-list()

  # Read in data ----------------------------------------------------
  dge_a <- readRDS(DGE.IN.PATH)
  for (c in names(dge_a)) { 
    print(c)
    dge <- dge_a[[c]]
    
    # remove controls -------------------------------------------------
    dge <- dge[, dge$samples$matched.timepoint.age!=1]
    dge <- dge[, !is.na(dge$samples$age)]
    dge <- dge[, as.numeric(dge$samples$matched.timepoint.age)<=14]
    
    # make things factors ---------------------------------------------
    meta <- dge$samples
    meta$matched.individual <- factor(meta$matched.individual)
    meta$gender[which(meta$gender=="M")]<-1
    meta$gender[which(meta$gender=="F")]<-2
    meta$gender <- as.numeric(meta$gender)
    meta$matched.timepoint.age <- factor(meta$matched.timepoint.age)
    meta$age.months <- as.numeric(meta$matched.timepoint.age)*12
    meta$batch <- factor(meta$batch)
    
    X <- elspline(meta$age.months, n = 4)
    
    for(i in 1:4){
      meta[[paste0("X", i)]] <- X[,i]
    }
    
    form <- ~ (1|matched.individual) + X1 + X2 + X3 + X4 + (1 | batch) + gender 
    
    # estimate weights using linear mixed model of dream
    print("estimate weights using linear mixed model of dream")
    vobjDream = voomWithDreamWeights( dge, form, meta, suppressWarnings=TRUE)
    
    # Fit the dream model on each gene
    # By default, uses the Satterthwaite approximation for the hypothesis test
    print("fit extract variance partition model")
    vp = fitExtractVarPartModel( vobjDream, form, meta)
    
    print("fit dream model")
    fitmm = dream( vobjDream, form, meta)
    tt <- topTable(fitmm, which(startsWith(colnames(fitmm), "X")), number = Inf)
    
    all_tts_vpar[[c]] <- vp
    all_tts_spline[[c]]<-tt
    all_tts_fitmm[[c]]<-fitmm
  }
  saveRDS(all_tts_vpar,VAR.OUT.PATH)
  saveRDS(all_tts_spline,SPLINE.OUT.PATH)
}
