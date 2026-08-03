library(fgsea)
library(reshape2)
library(ggplot2)
library(tibble)
library(limma)
library(stringr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ComplexHeatmap)
library(ggrepel)
library(ggVennDiagram)

set.seed(0)

gene_sets<-readRDS(snakemake@input[["gene_sets"]])
ves_genes<-readRDS(snakemake@input[["ves_genes"]])
age_gene_df<-read.csv(snakemake@input[["age_trends"]])

age_gene_df$age_mean_coeff <- apply(age_gene_df[paste0("X",seq(1,8))], 1, function(x) x[which.max(abs(x))])
age_gene_df$age_trend<-"Age-increasing"
age_gene_df$age_trend[age_gene_df$age_mean_coeff<0]<-"Age-decreasing"
age_corr_genes<-age_gene_df$gene[(age_gene_df$adj.P.Val<=0.05)]

datasets <- list(
  "Persistent-VES" = ves_genes$persistent,
  "Age-correlated" = age_corr_genes
)

ggVennDiagram(datasets) + scale_fill_gradient(low = "white", high = "white") + coord_flip() + theme(legend.position = "none")
ggsave("data/output/Age_and_VES_Venn.pdf")

ves_age_corr_overlap<-list()
for (gs in names(ves_genes)){
  ves_age_corr_overlap[[paste0(gs,"_agecorr")]]<-intersect(ves_genes[[gs]],age_corr_genes)
}
ves_age_corr_overlap[["Age-increasing"]]<-intersect(ves_age_corr_overlap$persistent_agecorr,age_gene_df$gene[(age_gene_df$age_trend=="Age-increasing")])
ves_age_corr_overlap[["Age-decreasing"]]<-intersect(ves_age_corr_overlap$persistent_agecorr,age_gene_df$gene[(age_gene_df$age_trend=="Age-decreasing")])


df_up_dwn<-table(age_gene_df$age_trend[age_gene_df$gene %in% ves_age_corr_overlap$persistent_agecorr]) #%in% unique(unlist(enrichment_df$overlapGenes[enrichment_df$padj<=0.01
ggplot(melt(df_up_dwn),aes(x=Var1,y=value,fill=Var1))+geom_bar(stat = "identity")+theme_bw()+ylab("Leading-edge, persistent VES genes")+ theme(legend.position = "none") +xlab("")
ggsave("data/output/Age_and_VES_bars.pdf")

bg <- unique(c(unique(as.vector(unlist(gene_sets$reactome))),unique(as.vector(unlist(gene_sets$go.bp))),unique(as.vector(unlist(gene_sets$kegg))),unique(as.vector(unlist(gene_sets$btms))),unique(age_gene_df$gene)))

enrichment_list<-list()
for (gs in c("Age-increasing","Age-decreasing")){
  for (p in names(gene_sets)){
  # 2. Run Overrepresentation Analysis
  fg<-ves_age_corr_overlap[[gs]]
  enrichment_list[[paste0(gs,"_",p)]]<-fora(pathways = gene_sets[[p]], genes = fg, universe = bg)
  enrichment_list[[paste0(gs,"_",p)]]<-enrichment_list[[paste0(gs,"_",p)]][enrichment_list[[paste0(gs,"_",p)]]$padj<=0.05]
  enrichment_list[[paste0(gs,"_",p)]]$gene_set<-gs
  enrichment_list[[paste0(gs,"_",p)]]$pway<-p
  }
}

enrichment_df<-do.call(rbind,enrichment_list)
enrichment_df$nlog_pval<-(-log10(enrichment_df$pval))
enrichment_df<-enrichment_df[enrichment_df$padj<=0.01,]
#enrichment_df[,str_detect(enrichment_df$gene_set,"persist")]

ggplot(enrichment_df, aes(x = gene_set, y = pathway)) +
  geom_point(aes(size = -log10(padj), color = nlog_pval, shape = padj < .05)) +
  scale_color_gradient2(low = "blue",mid = "white",high = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Genesets that change with age and are persistent in VES")

ggsave("data/output/Age_and_VES_Pathways.pdf")

