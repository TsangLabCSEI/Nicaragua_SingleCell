library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(edgeR)
library(limma)
library(stringr)


# load pbulk and vpars tables
dge <- readRDS(snakemake@input[["bulk"]])
bulkVES <- read.csv(snakemake@input[["vpars"]])
df_gex <- read.csv(snakemake@input[["euclids"]],sep=",",row.names = 1)
sdrf <- read.csv(snakemake@input[["euclids_meta"]],sep="\t")


# min max scaling function
min_max_scale <- function(x) {
	    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

### correlate with NICA snapshot VES
dge_cpm<-cpm(dge$counts, log = F)
dfNICA<-as.data.frame(t(dge_cpm[,which(dge$samples$Age.years==6)]))
dfNICA$individual<-as.vector(dge$samples$Subject.ID)[which(dge$samples$Age.years==6)]
dfNICA_var<-as.data.frame(cbind(log(sapply(dfNICA %>% dplyr::select(-individual), var)),log(sapply(dfNICA %>% dplyr::select(-individual), mean))))
colnames(dfNICA_var)<-c("variance","mean")
model <- lm(variance ~ mean, data = dfNICA_var)
dfNICA_var$variance<-min_max_scale(as.vector(residuals(model)))

rownames(bulkVES)<-bulkVES$Gene
bulkVES<-bulkVES[intersect(rownames(dfNICA_var),bulkVES$Gene),]
dfNICA_var<-dfNICA_var[intersect(rownames(dfNICA_var),bulkVES$Gene),]

mean_expr<-rowMeans(dge_cpm)
mean_expr<-mean_expr[intersect(rownames(dfNICA_var),bulkVES$Gene)]
bulkVES$mean_blk<-log10(mean_expr)
model <- lm(Subject.ID ~ mean_blk, data = bulkVES)
bulkVES$Subject.ID<-min_max_scale(as.vector(residuals(model)))

VES_vs_var_age6=as.data.frame(cbind(bulkVES$Subject.ID,dfNICA_var$variance))
colnames(VES_vs_var_age6)<-c("bulkVES","Subject.Variance_Age6")
VES_vs_var_age6$gene<-rownames(dfNICA_var)
ggplot(VES_vs_var_age6,aes(x=bulkVES,y=Subject.Variance_Age6,label=gene))+geom_point()+ geom_abline(intercept = 0, slope = 1, color="blue", lwd=1) + geom_label_repel(data=VES_vs_var_age6[(VES_vs_var_age6$bulkVES>0.5)&(VES_vs_var_age6$Subject.Variance_Age6>0.6),],max.overlaps = 100)+ theme_bw()+stat_cor(method = "spearman")
ggsave("data/output/snapshot_NICAage6_corr.pdf")


### EUCLIDS
smpls<-sdrf$Source.Name[(sdrf$Characteristics.disease.=="normal")&(sdrf$Characteristics.age.>=24)&(sdrf$Factor.Value.cohort.=="EUCLIDS")]
ages<-sdrf$Characteristics.age.[(sdrf$Characteristics.disease.=="normal")&(sdrf$Characteristics.age.>=24)&(sdrf$Factor.Value.cohort.=="EUCLIDS")]
sexs<-sdrf$Characteristics.sex.[(sdrf$Characteristics.disease.=="normal")&(sdrf$Characteristics.age.>=24)&(sdrf$Factor.Value.cohort.=="EUCLIDS")]
df_meta<-as.data.frame(cbind(smpls,ages,sexs))
colnames(df_meta)<-c("Sample","Age.months","Sex")
df_meta$Age.months<-as.numeric(df_meta$Age.months)
df_meta$Sample<-str_replace(df_meta$Sample," ",".")
dfEUCmean<-log10(rowMeans(df_gex[,df_meta$Sample]))
dfEUCvar<-log10(apply(df_gex[,df_meta$Sample], 1, var, na.rm = TRUE))
dfEUC<-as.data.frame(cbind(dfEUCmean,dfEUCvar))
colnames(dfEUC)<-c("Mean","Var")
model <- lm(Var ~ Mean, data = dfEUC)
dfEUC$EUCLIDS_GEX_Var<-as.vector(residuals(model))
dfEUC$gene<-rownames(dfEUC)
dfEUC<-dfEUC[intersect(dfEUC$gene,bulkVES$Gene),]

mean_expr<-rowMeans(dge_cpm)
rownames(bulkVES)<-bulkVES$Gene
bulkVES<-bulkVES[intersect(dfEUC$gene,bulkVES$Gene),]
mean_expr<-mean_expr[intersect(dfEUC$gene,bulkVES$Gene)]
bulkVES$mean_blk<-log10(mean_expr)
model <- lm(Subject.ID ~ mean_blk, data = bulkVES)
bulkVES$Subject.ID<-as.vector(residuals(model))

dfEUCbulkVES<-cbind(bulkVES,dfEUC)
dfEUCbulkVES$EUCLICS_GEX_Var_minmax<-min_max_scale(dfEUCbulkVES$EUCLIDS_GEX_Var)
dfEUCbulkVES$Subject.ID<-min_max_scale(dfEUCbulkVES$Subject.ID)
ggplot(dfEUCbulkVES,aes(x=Subject.ID,y=EUCLICS_GEX_Var_minmax,label=gene))+geom_point()+ geom_abline(intercept = 0, slope = 1, color="blue", lwd=1) + geom_label_repel(data=dfEUCbulkVES[(dfEUCbulkVES$Subject.ID>0.65)&(dfEUCbulkVES$EUCLICS_GEX_Var_minmax>0.65),],max.overlaps = 100)+ theme_bw()+stat_cor(method = "spearman")
#ggplot(dfEUCbulkVES,aes(x=rank(Subject.ID),y=rank(EUCLICS_GEX_Var_minmax),label=gene))+geom_point()+ geom_abline(intercept = 0, slope = 1, color="blue", lwd=1)+geom_density_2d(bins=100)
ggsave("data/output/snapshot_EUCLIDS_corr.pdf")
