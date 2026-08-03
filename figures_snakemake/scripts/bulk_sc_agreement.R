library(reshape2)
library(ggplot2)
library(tibble)
library(limma)
library(stringr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggrepel)

### VES correlation
VESgenes<-readRDS(snakemake@input[["ves_genes"]])
scvpars<-readRDS(snakemake@input[["vpars"]])
scvpars<-scvpars %>% filter(celltype %in% names(which(table(scvpars$celltype)>10000)))
scvpars<-scvpars %>% group_by(gene) %>% summarise(celltype_max = celltype[which.max(subject_variance_explained)], scVESmax=max(subject_variance_explained),scVESmedian=median(subject_variance_explained),scVESmean=mean(subject_variance_explained)) %>% as.data.frame()
rownames(scvpars)<-scvpars$gene

vpars_bulk<-read.csv(snakemake@input[["bulk_ves"]])
rownames(vpars_bulk)<-vpars_bulk$Gene
scvpars$bulkVES<-vpars_bulk[scvpars$gene,]$Subject.ID
scvpars<-scvpars[!is.na(scvpars$bulkVES),]
scvpars$VES<-""
scvpars$VES[scvpars$gene %in% VESgenes$persistent]<-"persistent"
ggplot(data = scvpars, aes(x = bulkVES, y = scVESmax,label=gene, colour=VES)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", lwd = 1) +
  geom_density_2d(bins = 25) +
  theme_bw()+
  stat_cor(method = "spearman")+
  geom_label_repel(data=scvpars[(scvpars$scVESmax<0.2)&(scvpars$bulkVES>0.2),],max.overlaps = 100,show.legend = F)

ggsave("data/output/bulk_vs_sc_VES_correlation.pdf")

scvpars$VES_level<-"VES < 0.1"
scvpars$VES_level[scvpars$bulkVES>0.1]<-"0.1 < VES > 0.25"
scvpars$VES_level[scvpars$bulkVES>0.25]<-"0.25 < VES > 0.5"
scvpars$VES_level[scvpars$bulkVES>0.5]<-"VES > 0.5"

ggplot(data = scvpars, aes(x = bulkVES, y = scVESmax,label=gene, colour=VES_level)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", lwd = 1) +
  theme_bw()+
  stat_cor(method = "spearman")

ggsave("data/output/bulk_vs_sc_VES_correlation_segmented.pdf")
