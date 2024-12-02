library(vioplot)
library(fmsb)
library(multcomp)
library(RColorBrewer)
library(scales)


### compute average single cell subject variance for a gene signature
compute_signature_stability<-function(signature){
  subj_vars<-list()
  for (ct in unique(gene_stability_sc_avg$celltype)){
    df_svar<-gene_stability_sc_avg[gene_stability_sc_avg$celltype==ct,]
    subj_vars[[paste0(ct,"-background")]]<-df_svar$subject_variance_explained
    df_svar<-gene_stability_sc_avg[gene_stability_sc_avg$celltype==ct,][signature,]
    df_svar<-df_svar[!is.na(df_svar$batch),]
    subj_vars[[paste0(ct,"-foreground")]]<-df_svar$subject_variance_explained
  }
  df_svar_ct<-melt(do.call(cbind,subj_vars))
  df_svar_ct
  df_svar_ct$Var3<-sapply(df_svar_ct$Var2, function(x){str_split_i(x,"-",2)})
  df_svar_ct$Var2<-sapply(df_svar_ct$Var2, function(x){str_split_i(x,"-",1)})
  return(df_svar_ct)
}


#plot violins of subject variance foreground vs background
plot_signature_stability<-function(signature_name,signature_stability_df){
  vioplot(value~Var2, data=signature_stability_df[signature_stability_df$Var3=="foreground",], col = "orange", plotCentre = "line", side = "left", horizontal=FALSE, las=2, xlab="", ylab=paste0(signature_name," subject variance"))
  vioplot(value~Var2, data=signature_stability_df[signature_stability_df$Var3=="background",], col = "lightblue", plotCentre = "line", side = "right", horizontal=FALSE, las=2, add=TRUE, ylab=paste0(signature_name," subject variance"))
  legend("topright", legend=c(signature_name, "background"),fill=c("orange", "lightblue"), cex = 0.6)
}


### calculate subject variance rate of change
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
    svar_sig_df_list_ct[[paste0(sig,"_",ct)]]<-svar_df_sig
  }
  svar_rate_sig_df<-do.call(rbind, svar_sig_df_list_ct)
  return(svar_rate_sig_df)
}


### calculate subject variance slope
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


### compute subject variance slope for a signature per celltype
compute_signature_stability_slope<-function(signature_name,signature_gene_list){
  svar_sig_rate_df<-compute_signature_stability_rateofchange(signature_name,signature_gene_list)
  svar_sig_slope_list<-list()
  for (ct in unique(svar_sig_rate_df$celltype)){
      df<-svar_sig_rate_df[svar_sig_rate_df$celltype==ct,]
      form <- mean_vexp ~ window_number
      df_lm<-as.data.frame(t(sapply(unique(df$gene),function(x){ calclmodp(x,df) })))
      df_lm$celltype<-ct
      df_lm$signature<-signature_name
      df_lm$gene<-rownames(df_lm)
      svar_sig_slope_list[[ct]]<-df_lm
  }
  df_stability_slope<-do.call(rbind,svar_sig_slope_list)
  return(df_stability_slope)
}


###plot radar of signature variance rate of change compostion per celltype
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
  legend(x=1.1, y=1.3, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , cex=0.8, pt.cex=1.5,col=hue_pal()(3))
  title(sig)
}

###plot rateOfchange of signature subject variance per celltype
plot_signature_agetrend<-function(signature_rateOfchange, cell_subsets){
  sig<-unique(signature_rateOfchange$signature)
  ggplot(data=signature_rateOfchange[signature_rateOfchange$celltype %in% cell_subsets,], aes(x=window_number,y=mean_vexp, group=celltype, colour=celltype))+ geom_smooth(method = "loess") + theme_bw()+ggtitle(paste0(sig," - ",paste(cell_subsets,collapse = " / ")))
}

###plot rateOfchange of signature subject variance per gene
plot_signature_gene_agetrend<-function(signature_rateOfchange, cell_subsets){
  sig<-unique(signature_rateOfchange$signature)
  ggplot(data=signature_rateOfchange[signature_rateOfchange$celltype %in% cell_subsets,], aes(x=window_number,y=mean_vexp, group=gene, colour=gene))+ geom_smooth(method = "loess", se = F) + theme_bw()+ggtitle(paste0(sig," - ",paste(cell_subsets,collapse = " / ")))
}
