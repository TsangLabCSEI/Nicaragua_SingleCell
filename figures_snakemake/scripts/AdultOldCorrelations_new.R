library(ggplot2)
library(stats)
library(matrixStats)
library(dplyr)


### Aging cohort: Adult / Old correlations

#load bootstrapped variance partition of adult cohort
#adults: bootstraps[[1-10]], #old bootstraps[[11-20]]
nind_list<-readRDS(snakemake@input[["nind_list"]])
ultrastab<-readRDS(snakemake@input[["ultrastab"]])
ultrastab_clust_vpars<-readRDS(snakemake@input[["ultrastab_clust_vpars"]])

age_corr_list<-list()
age_corrp_list<-list()

for (ct in intersect(unique(ultrastab_clust_vpars$celltype),names(nind_list))){
  print(ct)
  ultra<-ultrastab[[ct]]
  for (age in c("Adult","Old")){
   print(age)
   pseudo.sce.list<-nind_list[[ct]]
   nind_markers<-rownames(pseudo.sce.list[[11]])
   #adults: bootstraps[[1-10]], #old bootstraps[[11-20]]
   if (age=="Adult"){
      idxs<-seq(1,10)
   } else if (age=="Old"){
      idxs<-seq(11,20)
   }
   df_nind<-as.data.frame(cbind(pseudo.sce.list[[idxs[1]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[2]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[3]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[4]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[5]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[6]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[7]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[8]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[9]]][nind_markers,"subject_variance_explained"],pseudo.sce.list[[idxs[10]]][nind_markers,"subject_variance_explained"]))
   rownames(df_nind)<-nind_markers
   colnames(df_nind)<-c("A","B","C","D","E","F","G","H","I","J")
   df_nind<-rowMedians(as.matrix(df_nind))
   names(df_nind)<-nind_markers
   df_nind<-df_nind[which(!is.na(df_nind))]
   df_nind_adult<-df_nind
    
   ultrastab_clust_vpars_ct<-ultrastab_clust_vpars[ultrastab_clust_vpars$celltype==ct,]
   ultrastab_clust_vpars_ct<-as.data.frame(ultrastab_clust_vpars_ct %>% group_by(gene) %>% summarize(max(subject_variance_explained)))
   colnames(ultrastab_clust_vpars_ct)<-c("gene","subject_variance_explained")
   rownames(ultrastab_clust_vpars_ct)<-ultrastab_clust_vpars_ct$gene
   ultrastab_clust_vpars_ct<-ultrastab_clust_vpars_ct[nind_markers,]
   ultrastab_clust_vpars_ct<-ultrastab_clust_vpars_ct[which(!is.na(ultrastab_clust_vpars_ct[,1])),]
   df_nind<-ultrastab_clust_vpars_ct$subject_variance_explained
   names(df_nind)<-ultrastab_clust_vpars_ct$gene
   df_nind_young<-df_nind
   df_nind_adult<-df_nind_adult[names(df_nind_young)]
   df_nind_young<-df_nind_young[!is.na(df_nind_adult)]
   df_nind_adult<-df_nind_adult[!is.na(df_nind_adult)]
   df_nind<-as.data.frame(cbind(df_nind_adult,df_nind_young))
   colnames(df_nind)<-c("adult","young")
   test<-cor.test(df_nind$adult, df_nind$young,method = "spearman")
   ptest<-cor.test(df_nind$adult, df_nind$young)
   mod <- lm(df_nind$adult ~ df_nind$young + 0)
   cf <- coef(mod)

   df_nind$age_group<-age
   df_nind$pval<-test$p.value
   df_nind$corr<-as.vector(test$estimate)
   df_nind$celltype<-ct
   df_nind$ultra<-"normal"
   df_nind$gene<-rownames(df_nind)
   df_nind$ultra[df_nind$gene %in% ultra]<-"ultrastab"
   age_corr_list[[paste0(ct,"-",age)]]<-df_nind
   age_corrp_list[[paste0(ct,"-",age)]]<-c(as.vector(test$estimate),test$p.value,ptest$estimate,as.numeric(cf[1]),ct,age)
  }
}
age_corr_list_df<-do.call(rbind,age_corr_list)
age_corrp_list_df<-as.data.frame(do.call(rbind,age_corrp_list))
colnames(age_corrp_list_df)<-c("rho","p-val","pearson","slope","celltype","age")
age_corrp_list_df<-age_corrp_list_df[age_corrp_list_df$celltype %in% c("Mono_Classical","CD8_EM","CD8_TEMRA","CD8_Naive"),]
age_corr_list_df$ct_age_grp<-paste0(age_corr_list_df$celltype,"-",age_corr_list_df$age_group)
age_corr_list_df.filt<-age_corr_list_df[age_corr_list_df$ct_age_grp %in% c("Mono_Classical-Adult","CD8_EM-Adult","Mono_Classical-Old","CD8_EM-Old"),]
age_corr_list_df.filt$ct_age_grp_fct<-factor(age_corr_list_df.filt$ct_age_grp,levels=c("Mono_Classical-Adult","CD8_EM-Adult","Mono_Classical-Old","CD8_EM-Old"))
ggplot(data=age_corr_list_df.filt[age_corr_list_df.filt$ultra=="normal",],aes(y=adult,x=young)) + geom_point(colour="black") + geom_density_2d(bins=150) + geom_abline(intercept = 0, slope = 1, color="blue", lwd=1) + geom_point(data=age_corr_list_df.filt[age_corr_list_df.filt$ultra=="ultrastab",],aes(y=adult,x=young), colour="red") + geom_smooth(method="lm",formula=y~x+0, se=F, lty=2, color="orange", lwd=1) + facet_wrap(vars(ct_age_grp_fct),nrow = 2) + theme_light()
ggsave("data/output/AdultOld_Age_Correlations.pdf")


age_corr_list_df.filt$delta<-age_corr_list_df.filt$adult-age_corr_list_df.filt$young
age_corr_list_df.filt<-age_corr_list_df.filt[order(age_corr_list_df.filt$delta),]
delta_adults_df<-rbind(head(age_corr_list_df.filt[which((age_corr_list_df.filt$ultra=="ultrastab")&(age_corr_list_df.filt$age_group=="Adult")),], 5),tail(age_corr_list_df.filt[which((age_corr_list_df.filt$ultra=="ultrastab")&(age_corr_list_df.filt$age_group=="Adult")),], 5))
delta_adults_df$gene<-factor(delta_adults_df$gene, levels = delta_adults_df$gene)
ggplot(data=delta_adults_df,aes(y=gene,x=delta, fill=celltype))+geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values=c("black","grey")) + xlab("delta VES [young-adult]") + ylab("") + ggtitle("Adult")
ggsave("data/output/Adult_Age_Correlations_TopGenes.pdf")

delta_olds_df<-rbind(head(age_corr_list_df.filt[which((age_corr_list_df.filt$ultra=="ultrastab")&(age_corr_list_df.filt$age_group=="Old")),], 5),tail(age_corr_list_df.filt[which((age_corr_list_df.filt$ultra=="ultrastab")&(age_corr_list_df.filt$age_group=="Old")),], 5))
delta_olds_df$gene<-factor(delta_olds_df$gene, levels = delta_olds_df$gene)
ggplot(data=delta_olds_df,aes(y=gene,x=delta, fill=celltype))+geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values=c("black","grey")) + xlab("delta VES [young-old]") + ylab("") + ggtitle("Old")
ggsave("data/output/Old_Age_Correlations_TopGenes.pdf")


### Aging cohort: Young / Adult / Old VES boxplots
print("boxplots")
#age_corr_list_df_ed<-readRDS(snakemake@input[["age_corr_list_df_ed"]])
age_corr_list_df_ed<-age_corr_list_df
age_corr_list_df_ed<-age_corr_list_df_ed[age_corr_list_df_ed$age_group=="Adult",]
age_corr_list_df_ed$adult<-age_corr_list_df_ed$young
age_corr_list_df_ed$young<-NULL
age_corr_list_df_ed<-age_corr_list_df_ed[c("adult","age_group","gene","celltype","ultra")]
age_corr_list_df_ed$age_group<-"Young"
age_corr_list_df_ed<-rbind(age_corr_list_df[c("adult","age_group","gene","celltype","ultra")],age_corr_list_df_ed)
age_corr_list_df_ed$age_group<-factor(age_corr_list_df_ed$age_group,levels=c("Young","Adult","Old"))
ggplot(data=age_corr_list_df_ed[(age_corr_list_df_ed$ultra=="ultrastab"),],aes(x=celltype,y=adult,fill=age_group))+geom_boxplot()+theme_bw() + geom_point(shape=1, position = position_jitterdodge(jitter.width = 0.1), alpha=0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("VES [ultrastable genes]")
ggsave("data/output/AdultOld_Age_Correlations_Boxplot.pdf",width = 10)
ggplot(data=age_corr_list_df_ed,aes(x=celltype,y=adult,fill=age_group))+geom_violin(scale="width",draw_quantiles = c(0.25, 0.5, 0.75))+theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("VES [all genes]")
ggsave("data/output/AdultOld_Age_Correlations_Boxplot_all.pdf",width = 10)

