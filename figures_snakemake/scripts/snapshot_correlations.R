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
#add supplement6
ageTraj <- read.csv(snakemake@input[["agetraj"]])

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


### Finnish cohort
gse <- getGEO("GSE30211", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]] #first one is Affy data

# Access expression data and phenotype data
exprs_data <- exprs(eset)
pheno_data <- pData(eset)

symbols <- AnnotationDbi::select(hgu219.db,
                                 keys = rownames(exprs_data),
                                 columns = c("SYMBOL", "ENTREZID"),
                                 keytype = "PROBEID"
)

#if multiple gene names/IDs map to single probe, just take the first gene name
symbols <- symbols %>% distinct(PROBEID, .keep_all = T)
symbols_vec <- setNames(symbols$SYMBOL, symbols$PROBEID)

#collapse multiple probes by taking median
exprs_data <- as.data.frame(exprs_data)

exprs_data$probe <- rownames(exprs_data)
exprs_data$symbol <- symbols_vec[rownames(exprs_data)]
exprs_data <- relocate(exprs_data, probe, symbol, .before = GSM1063599)

exprs_data_collapse <- dplyr::group_by(exprs_data, symbol) %>% dplyr::summarize(across(starts_with("GSM"), median))

#remove missing values
exprs_data_collapse <- filter(exprs_data_collapse, !is.na(symbol) & !(symbol == ""))

exprs_data_collapse <- tibble::column_to_rownames(exprs_data_collapse, "symbol")
exprs_mat <- as.matrix(exprs_data_collapse)


pheno_data_short<-pheno_data[c("age at sample (months):ch1","gender:ch1","supplementary_file")]
pheno_data_short$sample_name<-str_split_i(str_split_i(str_split_i(pheno_data_short$supplementary_file,"suppl",2),"_",2),"\\.",1)
pheno_data_short$sample_rep<-str_split_i(str_split_i(str_split_i(pheno_data_short$supplementary_file,"suppl",2),"_",3),"\\.",1)
pheno_data_short<-pheno_data_short[str_detect(pheno_data_short$sample_name,"Control"),]
pheno_data_short$`age at sample (months):ch1`<-as.numeric(pheno_data_short$`age at sample (months):ch1`)
pheno_data_short<-pheno_data_short[pheno_data_short$`age at sample (months):ch1`>=6,]
pheno_data_short<-pheno_data_short[pheno_data_short$sample_name %in% names(table(pheno_data_short$sample_name)[which(table(pheno_data_short$sample_name)>1)]),]
exprs_data_flt<-exprs_mat[,rownames(pheno_data_short)]

pheno_data_short_plt<-pheno_data_short
pheno_data_short_plt$age.years<-pheno_data_short_plt$`age at sample (months):ch1`/12
ggplot(pheno_data_short_plt,aes(x=age.years,y=sample_name,color=`gender:ch1`))+geom_point()+theme_bw()+scale_x_continuous(breaks = seq(1,20,2))+theme(legend.title = element_blank())

fit_gene_wise <- function(dge){
  data<-dge[c("GEX","Sex","Age","Subject")]
  data$Sex<-as.numeric(as.factor(data$Sex))
  myformula <- as.formula("GEX ~ Age + Sex + (1|Subject)")
  tryCatch(expr = {
    fit <- lmer(myformula, data = data)
    calcVarPart(fit)
  },error = function(e) {
    return(NA)})
}

gene_list<-list()
for (gene in rownames(exprs_data_flt)){
  exprs_data_flt_mrgd<-cbind(exprs_data_flt[gene,],pheno_data_short)
  colnames(exprs_data_flt_mrgd)<-c("GEX","Age","Sex","file","Subject","Sample")
  varpart_result<-fit_gene_wise(exprs_data_flt_mrgd)
  gene_list[[gene]]<-varpart_result
}
df_varp<-do.call(rbind,gene_list)
df_varp<-as.data.frame(df_varp)
df_varp$Gene<-rownames(df_varp)


bulkVES <- read.csv(snakemake@input[["vpars"]])
rownames(bulkVES)<-bulkVES$Gene
genes<-intersect(bulkVES$Gene,rownames(df_varp))
bulkVES<-bulkVES[genes,]
df_varp<-df_varp[genes,]
df_varp_mltd<-reshape2::melt(df_varp)
df_varp$Subject.NICA<-bulkVES$Subject.ID
ggplot(df_varp,aes(x=Subject,y=Subject.NICA,label=Gene))+geom_point()+geom_abline(intercept = 0, slope = 1, color="blue", lwd=1) + geom_label_repel(data=df_varp[(df_varp$Subject>0.6)&(df_varp$Subject.NICA>0.6),],max.overlaps = 100)+geom_density_2d(bins=100)+ theme_bw()+stat_cor(method = "spearman")


### Finnish linear model
inf_marker<-c("ADAR","IP6K2","PML","DDX58","IRF7","PSMB8","DDX60","IRF9","RSAD2","EIF2AK2","ISG15","SOCS1","EIF4A3","ISG20","SP100","HERC5","KPNA2","STAT1","IFI27","MT2A","STAT2","IFI35","MX1","TAP1","IFI6","MX2","TRIM21","IFIH1","OAS1","TRIM22","IFIT1","OAS2","TRIM38","IFIT2","OAS3","TRIM5","IFIT3","OASL","UBE2L6","IFITM1","PARP9","USP18","IFITM3","PLSCR1","XAF1")

strd_marker<-c("ABCA2","NEFL","ABCA3","NR3C2","ABHD2","PTAFR","BCL2","PTCH1","CATSPER1","PTGER2","DEFA3","RARG","EEF2","RORA","FBXO32","RORC","KLF9","RXRA","MBD4","SGK1","MBP","SPP1","NCOA4","WNT7A")

cluster_gene_df<-ageTraj #read.csv("/home/lb2336/project/NICA/supp_tables_S6.csv")
c1_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==1),]$Gene
c2_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==2),]$Gene
c3_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==3),]$Gene
c4_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==4),]$Gene
c5_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==5),]$Gene
c6_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==6),]$Gene
c7_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==7),]$Gene
c8_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==8),]$Gene

#create list of cluster genesets
genesetl<-list()
genesetl[["C1: Bcells,\n lymphocyte activation"]]<-c1_markers
genesetl[["C2: tRNA processing,\n hydrogen peroxide,\n hemoglobin,translation"]]<-c2_markers
genesetl[["C3: Steroid hormones,\n leukocyte adhesion"]]<-c3_markers
genesetl[["C3: Steroid hormone-specific"]]<-strd_marker
genesetl[["C4: Ox.reductase activity"]]<-c4_markers
genesetl[["C5: Cell cycle,division,\n interferon response"]]<-c5_markers
genesetl[["C5: Interferon-specific"]]<-inf_marker
genesetl[["C6: Adhesion, proliferation,\n cytoskeleton"]]<-c6_markers
genesetl[["C7: Integrins, hormones"]]<-c7_markers
genesetl[["C8: Ig heavy and light chain,\n B cells, cell cycle"]]<-c8_markers

df_meta<-pheno_data_short
colnames(df_meta)<-c("Age.months","Sex","file","subject","rep")
design <- model.matrix(~Age.months + Sex, df_meta)
fit <- lmFit(exprs_data_flt, design = design)
fit <- eBayes(fit)
toptab<-topTable(fit, coef = "Age.months", number = Inf) %>% rownames_to_column("gene")
fc <- toptab$logFC
names(fc) <- toptab$gene
set.seed(1)
fgseaRes <- fgsea(pathways = genesetl, stats = fc, minSize=15)
fgseaRes$nlog_pval <- -log10(fgseaRes$padj)
fgsea_dat <- fgseaRes
fgsea_dat$cohort<-"Bulk.FINN"
fgsea_dat %>%
  ggplot(aes(x = cohort, y = pathway)) +
  geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_gradient2(low = "blue",mid = "white",high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Age-trajectories - FINN cohort") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
