library(fgsea)
library(reshape2)
library(ggplot2)
library(tibble)
library(limma)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
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


### Age * Sex interactions

summ_list<-readRDS(snakemake@input[["bulk_windowed_analysis"]])
summ <- summ_list$sex_div_subj_age_sex_resid

fit_gene_wise <- function(dge){
  result<-list()
  data<-dge[c("mean_vexp","window_number","window")]
  myformula <- as.formula("mean_vexp ~ window_number")
  tryCatch(expr = {
    fit <- lm(myformula, data = data)
    data$vexp_fitted<-predict(fit,data)
    res <- summary(fit)
    result[["pval"]] <- as.vector(res$coefficients["window_number",][4])
    result[["t"]] <- as.vector(res$coefficients["window_number",][3])
    result[["prdxns"]] <- data
    return(result)
  },error = function(e) {
  return(NA)
  })
}

gene_list<-list()
gene_list_t<-list()
gene_list_fitted_expr<-list()
for (gene in unique(summ$gene)){
  varpart_result<-fit_gene_wise(summ[summ$gene==gene,])
  gene_list[[gene]]<-varpart_result[["pval"]]
  gene_list_t[[gene]]<-varpart_result[["t"]]
  gene_list_fitted_expr[[gene]]<-varpart_result[["prdxns"]]
  gene_list_fitted_expr[[gene]]$gene<-gene
}

sex_puberty_hits<-c("RPS4Y1", "EIF1AY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "TMSB4Y", "PRKY", "CD99", "PSMC3IP", "INHBB", "SPAG1", "TP63", "PRKX","XIST", "GPER1")

dfm_sel<-do.call(rbind,gene_list_fitted_expr[sex_puberty_hits])
dfm_sel <- dfm_sel %>%
  group_by(gene) %>%
  mutate(delta_vexp = mean_vexp - mean_vexp[window_number == 1])

my_matrix<-as.data.frame(pivot_wider(dfm_sel, id_cols = c("window_number"), names_from = gene, values_from = "delta_vexp",values_fn = mean))
rownames(my_matrix)<-my_matrix$window_number
my_matrix$window_number<-NULL
my_matrix<-t(my_matrix)

row_cluster <- hclust(dist(my_matrix))
row_order <- row_cluster$order
# 2. Reorder factor levels in long data
dfm_sel$gene <- factor(dfm_sel$gene, levels = rownames(my_matrix)[row_order])

ggplot(dfm_sel[dfm_sel$gene %in% sex_puberty_hits,],aes(x=window_number,y=gene,fill=delta_vexp))+geom_tile()+theme_bw()+scale_fill_gradient2(low = "blue",mid = "white",high = "red", name="delta VE-Sex")+ggtitle("Variance-explained by sex per sliding age-window")+ylab("Genes")
ggsave("data/output/Age_and_Sex_interaction.pdf")


### Age * Age interactions / VEA per window

age_gene_df<-read.csv(snakemake@input[["bulk_vea"]])                                                                    

summ <- summ_list$age_div_subj_age_gender_resid
age_corr_genes<-age_gene_df$Gene[age_gene_df$Age.months>=0.25]
summ <- summ[summ$gene %in% age_corr_genes,]
summ_agg<-summ %>% group_by(window_number) %>% summarize(mean_vea=mean(mean_vexp),var_vea=var(mean_vexp),quantile(mean_vexp,0.75))
colnames(summ_agg)<-c("window_number","mean_vea","var_vea","upper_quantile")

ggplot(summ,aes(x=window_number,mean_vexp,group=window_number))+geom_boxplot()+geom_hline(yintercept = 0.012,linetype = "dashed",color="red")+geom_signif(comparisons = list(c(1,3), c(1,7), c(1,12)),y_position = c(0.1, 0.105, 0.11),map_signif_level = TRUE, textsize = 4)+theme_bw()+ggtitle("Variance-explained by age per sliding age-window")+ylab("VE-Age")
ggsave("data/output/VEAge_per_window.pdf")

