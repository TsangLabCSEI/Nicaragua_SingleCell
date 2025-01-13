library(vioplot)
library(fmsb)
library(multcomp)
library(RColorBrewer)
library(scales)


#' Compute Signature Stability
#'
#' This function computes the average single cell subject variance for a gene signature.
#' Input is a vector of gene symbols and it retrieves the distribution of its subject
#' variance to the distribution of all genes within each celltype as background.
#'
#' @param signature vector of gene symbols
#' @return A data frame of subject variance per gene and celltype
compute_signature_stability<-function(signature){
  subj_vars<-list()
  for (ct in unique(gene_stability_sc_avg$celltype)){
    df_svar<-gene_stability_sc_avg[gene_stability_sc_avg$celltype==ct,]
    subj_vars[[paste0(ct,"-background")]]<-cbind(df_svar$subject_variance_explained,ct,"background")
    df_svar<-gene_stability_sc_avg[gene_stability_sc_avg$celltype==ct,][signature,]
    df_svar<-df_svar[!is.na(df_svar$batch),]
    subj_vars[[paste0(ct,"-foreground")]]<-cbind(df_svar$subject_variance_explained,ct,"foreground")
  }
  df_svar_ct<-as.data.frame(do.call(rbind,subj_vars))
  colnames(df_svar_ct)<-c("subject_variance","celltype","subset")
  df_svar_ct$subject_variance<-as.numeric(df_svar_ct$subject_variance)
  return(df_svar_ct)
}


#' Plot Violins of Subject Variance
#'
#' This function takes a signature name and single cell subject variance data frame as
#' input and visualizes the distribution of its subject variance in contrast to the 
#' background in a violinplot.
#'
#' @param signature_name string of signature name
#' @param signature_stability_df data frame of subject variance per gene and celltype
#' @return violinplot
plot_signature_stability<-function(signature_name,signature_stability_df){
  vioplot(subject_variance~celltype, data=signature_stability_df[signature_stability_df$subset=="foreground",], col = "orange", plotCentre = "line", side = "left", horizontal=FALSE, las=2, xlab="", ylab=paste0(signature_name," subject variance"))
  vioplot(subject_variance~celltype, data=signature_stability_df[signature_stability_df$subset=="background",], col = "lightblue", plotCentre = "line", side = "right", horizontal=FALSE, las=2, add=TRUE, ylab=paste0(signature_name," subject variance"))
  legend("topright", legend=c(signature_name, "background"),fill=c("orange", "lightblue"), cex = 0.6)
}


#' Heatmap of Subject Variance per Signature
#'
#' This function takes a vector of signature names as input and 
#' visualizes the distribution of its subject variance in contrast 
#' to the background in a violinplot.
#'
#' @param signature_names vector of signature names
#' @return heatmap
plot_signature_stability_heatmap<-function(signature_names){
  stability_df_list<-list()
  for (sig in signature_names){
        df<-compute_signature_stability(baseline_signatures[[sig]])
        df$signature<-sig
        stability_df_list[[sig]]<-df
  }
  stability_df<-do.call(rbind,stability_df_list)
  stability_df<-stability_df[stability_df$subset=="foreground",]
  sig_subj_vars_all_fg<-as.data.frame(stability_df %>% group_by(celltype,signature) %>% summarize(median(subject_variance)))
  colnames(sig_subj_vars_all_fg)<-c("celltype","signature","subject_variance")
  celltypes<-c("bulk","B_Mem","B_Naive","CD4_Mem","CD4_Naive","CD8_Mem","CD8_Naive","ILC","MAIT","gdT_Vd1","gdT_Vd2","Mono_Classical","Mono_NonClassical","NK_CD16hi","NK_CD56hi","cDC","pDC")
  sig_subj_vars_all_fg$celltype<-factor(sig_subj_vars_all_fg$celltype, levels=celltypes)
  ggplot(data=sig_subj_vars_all_fg[sig_subj_vars_all_fg$celltype!="ILC",],aes(x=signature,y=celltype,fill=subject_variance))+geom_tile()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis_c(option = "magma")
}


#' Calculate Rate of Change Subject Variance 
#'
#' This function takes a signature name and a vector of signature gene symbols as
#' input and retrieves the estimated subject variance in sliding age windows of
#' size 3 in bulk and size 5 in single cell leading to the rate of change across
#' age. It returns a data frame of mean and median subject variance per window, 
#' gene and celltype.
#' 
#' @param signature_name string of signature name
#' @param signature vector of gene symbols
#' @return data frame of subject variance rate of change per gene and celltype
compute_signature_stability_rateofchange<-function(signature_name,signature){
  svar_sig_df_list_ct<-list()
  for (ct in names(gene_stability)){
    svar_df<-gene_stability[[ct]]
    svar_df_sig<-as.data.frame(svar_df[svar_df$gene %in% signature,])
    svar_df_sig$celltype<-ct
    svar_df_sig$signature<-signature_name
    svar_df_sig$n_genes<-length(unique(svar_df_sig$gene))
    svar_df_sig<-svar_df_sig[svar_df_sig$window_number<9,]
    svar_df_sig_window<-svar_df_sig %>% group_by(window_number) %>% summarise(mean_vexp_ct=mean(mean_vexp))
    svar_sig_df_list_ct[[ct]]<-svar_df_sig
  }
  svar_rate_sig_df<-do.call(rbind, svar_sig_df_list_ct)
  return(svar_rate_sig_df)
}


#' Calculate linear model and p-value
#'
#' This function takes a gene symbol and a subject variance rate of change data frame 
#' as input to subseting it for the specified gene. It then fits a linear regression
#' estimating the intercept and slope of subject variance and a p-value of significane
#' for the fit. These coefficients are returned as an output
#'
#' @param x string of gene name
#' @param df data frame of subject variance rate of change
#' @return linear model coefficients
calclmodp<-function(x,df){
  form <- mean_vexp ~ window_number
  lmod<-lm(form,df[df$gene==x,])
  coefm<-as.data.frame(t(coef(lmod)))
  res<-tryCatch({
    coefp<-summary(glht(lmod, linfct = c("window_number = 0")))$test$pvalues[1]
  }, error = function(e) {return(NA)})
  coefm$p<-res
  return(unlist(coefm))
}


#' Compute subject variance slope for a signature per celltype
#'
#' This function takes a signature name and a vector of signature gene symbols as
#' input and retrieves the estimated slope of subject variance across age windows.
#' It returns a data frame of the slope of subject variance per window gene and 
#' celltype.
#'
#' @param signature_name string of signature name
#' @param signature_gene_list vector of gene symbols
#' @return data frame of subject variance slope of rate of change
compute_signature_stability_slope<-function(signature_name,signature_gene_list){
  svar_sig_rate_df<-compute_signature_stability_rateofchange(signature_name,signature_gene_list)
  svar_sig_slope_list<-list()
  for (ct in unique(svar_sig_rate_df$celltype)){
      df<-svar_sig_rate_df[svar_sig_rate_df$celltype==ct,]
      df_lm<-as.data.frame(t(sapply(unique(df$gene),function(x){ calclmodp(x,df) })))
      df_lm$celltype<-ct
      df_lm$signature<-signature_name
      df_lm$gene<-rownames(df_lm)
      svar_sig_slope_list[[ct]]<-df_lm
  }
  df_stability_slope<-do.call(rbind,svar_sig_slope_list)
  return(df_stability_slope)
}


#' Plot Radar of Signature Variance Rate of Change Compostion per Celltype
#'
#' This function takes a data frame of the slope of subject variance per gene of
#' a gene signature and celltype as input and categorizes each gene into increasing 
#' (positive slope, p<=0.05) or decreasing (negative slope, p<=0.05) or stable 
#' (p>0.05). It then visualizes the fractions of these categories as radarplot per 
#' celltype.
#'
#' @param signature_stability_slope data frame of subject variance slope of rate of change
#' @return radarplot
plot_signature_stability_radar<-function(signature_stability_slope){
  incr_decr_list<-list()
  for (ct in unique(signature_stability_slope$celltype)){
    sig<-unique(signature_stability_slope$signature)
    decr<-sum(signature_stability_slope$window_number[((signature_stability_slope$celltype==ct)&(signature_stability_slope$p<=0.05))]<=0, na.rm = T)
    incr<-sum(signature_stability_slope$window_number[((signature_stability_slope$celltype==ct)&(signature_stability_slope$p<=0.05))]>0, na.rm = T)
    stab<-sum(abs(signature_stability_slope$p[(signature_stability_slope$celltype==ct)])>0.05, na.rm = T)
    incr_decr_list[[ct]]<-c(decr,incr,stab)
  }
  df_incr_decr<-do.call(rbind,incr_decr_list)
  colnames(df_incr_decr)<-c("decreasing","increasing","stable")
  df_incr_decr<-as.data.frame(t(df_incr_decr))
  for (name in names(df_incr_decr)){df_incr_decr[[name]]<-df_incr_decr[[name]]/colSums(df_incr_decr)[[name]]}
  data <- rbind(rep(1,dim(df_incr_decr)[2]) , rep(0,dim(df_incr_decr)[2]) , df_incr_decr)

  radarchart(data,pcol = hue_pal()(3))
  legend(x=1.1, y=1.35, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , cex=0.8, pt.cex=1.5,col=hue_pal()(3))
  title(sig)
}


#' Plot the Average Rate of Change of Signature Subject Variance per Celltype
#'
#' This function takes the data frame of subject variance rate of change for a gene
#' signature as input and visualizes the average subject variance trend across all
#' genes per celltype as a loess spline fit.
#'
#' @param signature_rateOfchange data frame of subject variance rate of change for gene signature
#' @param cell_subsets celltypes to be visualized
#' @return lineplot
plot_signature_agetrend<-function(signature_rateOfchange, cell_subsets){
  sig<-unique(signature_rateOfchange$signature)
  ggplot(data=signature_rateOfchange[signature_rateOfchange$celltype %in% cell_subsets,], aes(x=window_number,y=mean_vexp, group=celltype, colour=celltype))+ geom_smooth(method = "loess") + theme_bw()+ggtitle(paste0(sig," - ",paste(cell_subsets,collapse = " / ")))
}


#' Plot rateOfchange of signature subject variance per gene
#'
#' This function takes the data frame of subject variance rate of change for a gene
#' signature as input and visualizes the subject variance trend for each gene per 
#' celltype as a linear trend lines.
#'
#' @param signature_rateOfchange data frame of subject variance rate of change for gene signature
#' @param cell_subsets celltypes to be visualized
#' @return lineplot
plot_signature_gene_agetrend<-function(signature_rateOfchange, cell_subsets){
  sig<-unique(signature_rateOfchange$signature)
  ggplot(data=signature_rateOfchange[signature_rateOfchange$celltype %in% cell_subsets,], aes(x=window_number,y=mean_vexp, group=gene, colour=gene))+ geom_smooth(method = "loess", se = F) + theme_bw()+ggtitle(paste0(sig," - ",paste(cell_subsets,collapse = " / ")))
}
