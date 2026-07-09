library(stringr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(lme4)
library(lmerTest)
library(variancePartition)
library(matrixStats)
library(dplyr)


### read input
df_rq<-read.csv(snakemake@input[["RQ_data"]])
df_aq<-read.csv(snakemake@input[["AQ_data"]])
metadata<-read.csv(snakemake@input[["meta_data"]])
tech_meta<-read.csv(snakemake@input[["tech_data"]])

science_data<-read.csv(snakemake@input[["science_data"]])
science_data_age<-read.csv(snakemake@input[["science_age_data"]])
science_data_age<-science_data_age[science_data_age$term=="Age",]
AGEgenes_up<-science_data_age$Protein[(science_data_age$p.value<=0.001)&(science_data_age$statistic>0)] # GZMA, FGF21, VEGFD, (slightly FGF19,PGF) + NGF/ contrary EGF (VEGFA non-sig),PDGFA,PDGFB,TGFB1
AGEgenes_dn<-science_data_age$Protein[(science_data_age$p.value<=0.001)&(science_data_age$statistic<0)] # intersect and agree eg. IL2RA, IL22, CXCL8

bulkVES<-read.csv(snakemake@input[["supp9_data"]])
rownames(bulkVES)<-bulkVES$Gene


### reformat
df<-df_rq
df$SampleName_short<-str_split_i(df$SampleName,"P01_",-1)
df$SampleName_short<-as.vector(sapply(df$SampleName_short,function(x){gsub("[ATCDBVPX]+$", "", x)}))
df$SampleName_short<-as.vector(sapply(df$SampleName_short,function(x){str_replace(x,"_",".")}))
metadata$SampleName_short<-as.vector(sapply(toupper(metadata$Original.ID),function(x){gsub("[ATCDBVPX]+$", "", x)}))
metadata$Sex<-str_replace(metadata$Sex,"Male","M")
metadata$Sex<-str_replace(metadata$Sex,"Female","F")
dfm<-merge(df,metadata)
dfm$Subject.ID<-paste0("nica",dfm$Subject.ID)
dfm_grp<-as.data.frame(dfm %>% group_by(Subject.ID) %>% summarize(length(unique(Age))))
subjts<-dfm_grp$Subject.ID[dfm_grp$`length(unique(Age))`!=1]
dfm<-dfm[dfm$Subject.ID %in% subjts,]

overlap_genes_science<-intersect(unique(science_data$Protein),unique(df$Target))
science_data<-science_data[science_data$Protein %in% overlap_genes_science,]
science_data<-science_data[str_detect(science_data$Component,"subject"),]

dfm$Target_Sample<-paste0(dfm$Target,"_",dfm$SampleName)
dfm_tech<-merge(dfm,tech_meta,by = "Target_Sample")
#print(colnames(dfm_tech)[str_detect(colnames(dfm_tech),"NPQ")|str_detect(colnames(dfm_tech),"Plate")])
dfm<-dfm_tech


### cohort overview
age_df<-dfm[c("Subject.ID","Age","Sex")]
age_df<-age_df %>% group_by(Subject.ID) %>% mutate(min=min(Age),max=max(Age)) %>% ungroup()
age_df<-age_df[order(age_df$Age,age_df$max-age_df$min, decreasing = F),]
age_df$idx<-as.numeric(factor(age_df$Subject.ID,levels=unique(age_df$Subject.ID)))
age_df %>%
ggplot(aes(x=Age, y = idx, col=Sex, group=idx)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(breaks = seq(1,20,1)) +
  xlab("Followed age range")
ggsave("data/output/proteomics_cohort_overview.pdf")

### variance partition
fit_gene_wise <- function(dge){
  data<-dge[c("NPQ","Sex","Age","Subject.ID","PlateID","Sample_QC_Detectability","Sample_QC_ICReads","Sample_QC_NumReads")]
  data$Sex<-as.numeric(as.factor(data$Sex))
  myformula <- as.formula("NPQ ~ Age + Sex + Sample_QC_Detectability + Sample_QC_ICReads + Sample_QC_NumReads + (1|Subject.ID) + (1|PlateID)")
  tryCatch(expr = {
  fit <- lmer(myformula, data = data)
  calcVarPart(fit)
  },error = function(e) {
  return(NA)})
}

gene_list<-list()
for (gene in unique(dfm$Target)){
    varpart_result<-fit_gene_wise(dfm[dfm$Target==gene,])
    gene_list[[gene]]<-varpart_result
}
df_varp<-do.call(rbind,gene_list)
df_varp<-as.data.frame(df_varp)

### visualize variance partition
df_varp_mltd<-reshape2::melt(df_varp)
ggplot(df_varp_mltd,aes(x=variable,y=value,fill=variable))+geom_violin(scale = "width")+theme_bw()
ggsave("data/output/proteomics_variance_partition.pdf")

### correlate with science paper
df_varp_science<-as.data.frame(df_varp[overlap_genes_science,])
df_varp_science$varExp_subj_science<-science_data$Variance
colnames(df_varp_science)<-c("Plate","varExp_subj_nica","Age","Sex","Sample_QC_Detectability","Sample_QC_ICReads","Sample_QC_NumReads","Residuals","varExp_subj_science")
df_varp_science$GeneID<-rownames(df_varp_science)
ggplot(df_varp_science,aes(x=varExp_subj_nica,y=varExp_subj_science,label=GeneID))+geom_point()+ geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+geom_smooth(method="lm")+ geom_label_repel(data=df_varp_science[(df_varp_science$varExp_subj_nica>0.5)&(df_varp_science$varExp_subj_science>0.25),])+ stat_cor(method = "spearman") +theme_bw()
ggsave("data/output/proteomics_science_correlation.pdf")

### visualize top genes with highest variance explained by subject
df_varp<-do.call(rbind,gene_list)
df_varp<-as.data.frame(df_varp)
df_varp$gene<-rownames(df_varp)
df_varp <- df_varp %>% mutate(var_bio_sum = Sex + Age + Subject.ID + Residuals) %>% mutate(ves_renorm = Subject.ID / var_bio_sum, vex_renorm= Sex / var_bio_sum, vea_renorm= Age / var_bio_sum, res_renorm=Residuals / var_bio_sum)
df_varp<-df_varp[,c("ves_renorm","vea_renorm","vex_renorm","res_renorm","gene")]
genes<-unique(head(df_varp[order(df_varp$ves_renorm, decreasing = T),]$gene,30))
df_varp<-df_varp[(df_varp$gene %in% genes),]
colnames(df_varp)<-c("Subject","Age","Sex","Residuals","gene")
df_varp<-reshape2::melt(df_varp,id.vars="gene")
df_varp_subj<-df_varp[df_varp$variable=="Subject",]
df_varp_subj<-df_varp_subj[order(df_varp_subj$value, decreasing = T),]
df_varp<-rbind(df_varp_subj,df_varp[df_varp$variable!="Subject",])
df_varp$variable<-factor(df_varp$variable,levels=c("Residuals","Age","Sex","Subject"))

stable_nulisa_genes<-df_varp_subj$gene
ggplot(data=df_varp,aes(y=gene,x=value,group=variable,fill=variable))+geom_bar(stat='identity',colour="black",size=0.25)+ scale_y_discrete(limits = rev(stable_nulisa_genes))+ theme_classic() + labs(title="Sorted by Variance-explained by Subject (VES)") + scale_fill_manual(values=c("white","#619CFF","#00BA38","#F8766D"))
ggsave("data/output/proteomics_variance_partition_genes.pdf")

