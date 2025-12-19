# set seed
set.seed(4)

# load libraries
library(dplyr)
library(ggplot2)
library(stringr)

# load pbulk and vpars tables
dge <- readRDS(snakemake@input[["data"]]) 
vpars <- readRDS(snakemake@input[["vpars"]])

### Correlation analysis of mean expression and scVES per gene and celltype
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
corr_list<-list()
for (ct in names(dge)){
  dge_ct<-dge[[ct]]
  df<-rowMeans(dge_ct$normalizedExpr)
  genes<-names(df)
  vpars_ct<-vpars[vpars$celltype==ct,]
  genes<-intersect(genes,vpars_ct$gene)
  vpars_ct<-vpars_ct[vpars_ct$gene %in% genes,]
  vpars_ct<-vpars_ct[!duplicated(vpars_ct$gene),]
  rownames(vpars_ct)<-vpars_ct$gene
  genes<-intersect(genes,vpars_ct$gene)
  df_corr<-cbind(df[genes],vpars_ct[genes,]$subject_variance_explained)
  colnames(df_corr)<-c("mean_geneExpr","scVES")
  df_corr<-as.data.frame(df_corr)
  df_corr<-df_corr[!is.na(df_corr$scVES),]
  df_corr$mean_geneExpr<-range01(as.numeric(df_corr$mean_geneExpr))
  df_corr$scVES<-as.numeric(df_corr$scVES)
  df_corr$celltype<-ct
  corr_list[[ct]]<-df_corr
}
corrs<-do.call(rbind,corr_list)

ggplot(data=corrs,aes(y=scVES,x=mean_geneExpr)) + geom_point(colour="black",pch=1) + geom_density_2d(bins=25) + geom_abline(intercept = 0, slope = 1) + geom_smooth(method="lm",formula=y~x+0, se=F, lty=2, color="orange", lwd=1) + facet_wrap(vars(celltype),nrow = 5) + theme_light() + ylim(0,1.2) + stat_cor(label.x = 0.05, label.y = 1.15) + stat_regline_equation(label.x = 0.05, label.y = 1.01)

ggsave("data/output/GEX_VES_correlations.pdf")

