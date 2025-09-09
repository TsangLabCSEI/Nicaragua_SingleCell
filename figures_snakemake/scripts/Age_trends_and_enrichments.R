#load libraries
library(ggplot2)
library(fgsea)
library(edgeR)
library(reshape2) 
library(broom)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)


#load age-trajectory cluster data
cluster_gene_df<-read.csv("data/input/supp_tables_S6.csv")
c1_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==1),]$Gene
c2_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==2),]$Gene
c3_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==3),]$Gene
c4_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==4),]$Gene
c5_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==5),]$Gene
c6_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==6),]$Gene
c7_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==7),]$Gene
c8_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==8),]$Gene

# bulk RNAseq subject-variance enriched genesets
nrchd_gene_df<-read.csv("data/input/supp_tables_S10.csv") 


#create list of cluster genesets
genesetl<-list()
genesetl[["C1: Bcells,\n lymphocyte activation"]]<-c1_markers
genesetl[["C2: tRNA processing,\n hydrogen peroxide,\n hemoglobin,translation"]]<-c2_markers
genesetl[["C3: Steroid hormones,\n leukocyte adhesion"]]<-c3_markers
genesetl[["C4: Ox.reductase activity"]]<-c4_markers
genesetl[["C5: Cell cycle,division,\n interferon response"]]<-c5_markers
genesetl[["C6: Adhesion, proliferation,\n cytoskeleton"]]<-c6_markers
genesetl[["C7: Integrins, Hormones"]]<-c7_markers
genesetl[["C8: Ig heavy and light chain,\n B cells, cell cycle"]]<-c8_markers

#load variance partitioning results
vpars<-readRDS("data/input/NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_vpars_object.rds") # aggregated
all_tts<-readRDS("data/input/NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_spline.rds") #raw with stats

#load normalized pseudobulk scRNAseq data
dge <- readRDS("data/input/NICA_batch2to4_new_timepoints_cleaned_normalized.rds")


### Computation of age variance gene set enrichments

#set seed
set.seed(1)


### Correlation analysis of bulk cluster 3 gene module with age across celltypes
cluster_tt_means<-list()
for (ct in names(all_tts)[which(!(names(all_tts) %in% c("doublets","ILC")))]){
dge_c<-dge[[ct]]
dge_c$normalizedExpr <- dge_c$normalizedExpr[, which(as.numeric(dge_c$samples$matched.timepoint.age)<=14)]
tt<-all_tts[[ct]]
tt$cell.type<-ct
tt$gene<-row.names(tt)
overlap<-tt[which((tt$P.Value<=1)&(tt$gene%in%c3_markers)&(tt$cell.type==ct)),]$gene
tps<-as.vector(sapply(colnames(dge_c$normalizedExpr),function(x){as.numeric(str_split_i(x,"_",2))}))
inds<-as.vector(sapply(colnames(dge_c$normalizedExpr),function(x){str_split_i(x,"_",1)}))
dge_t<-as.data.frame(t(dge_c$normalizedExpr))
dge_t<-as.data.frame(scale(dge_t,scale = FALSE, center=TRUE))
dge_t$timepoints<-tps
dge_t$individuals<-inds
dge_t<-dge_t[which(!is.na(dge_t$timepoints)),]

dfm<-melt(dge_t,id.vars = "timepoints",measure.vars = overlap)#
df2 <- dfm %>% group_by(timepoints,variable) %>% summarise(across(value,median),.groups = 'drop') %>%as.data.frame()

data<-reshape(df2, idvar = "variable", timevar = "timepoints", direction = "wide")
if(!"value.17" %in% names(data)){
  data$value.17<-NA
}
if(!"value.16" %in% names(data)){
  data$value.16<-NA
}
if(!"value.15" %in% names(data)){
  data$value.15<-NA
}
if(!"value.2" %in% names(data)){
  data$value.2<-NA
}
data2<-data[,seq(2,17)]

data2<-sapply(data2,as.numeric)
colnames(data2)<-as.vector(sapply(colnames(data2),function(x){str_split_i(x,"value.",2)}))
row.names(data2)<-data$variable
cluster_tt_means[[ct]]<-colMeans(data2)
}
data3<-melt(do.call(cbind,cluster_tt_means))
colnames(data3)<-c("timepoints","celltype","value")
ggplot(
  data=data3, #c("B_Mem","B_Naive","CD4_Naive","CD4_Mem","CD8_Mem","CD8_Naive","MAIT","Mono_Classical","gdT_Vd2", "NK_CD16hi"))
  aes(x=timepoints, y=value,color=celltype)) + stat_smooth(geom='line', alpha=0.25, se=FALSE) + stat_smooth(data=data3[which(data3$celltype %in% c("B_Naive","B_Mem","CD4_Naive","CD4_Mem","CD8_Mem","CD8_Naive","MAIT","gdT_Vd1","gdT_Vd2","Mono_Classical","NK_CD16hi")),], geom='line', se=FALSE) + theme_classic()+ theme(text = element_text(size=15)) + scale_x_continuous(limits=c(3, 14),breaks=seq(3,14,2))+ggtitle("C3: Steroid response module across celltypes")+theme(legend.position = "none")
ggsave("data/output/C3_STRD_AgeTrend.pdf")


### Correlation analysis of bulk cluster 5 gene module with age across celltypes
cluster_tt_means<-list()
for (ct in names(all_tts)[which(!(names(all_tts) %in% c("doublets","ILC")))]){
dge_c<-dge[[ct]]
dge_c$normalizedExpr <- as.data.frame(dge_c$normalizedExpr)[, which(as.numeric(dge_c$samples$matched.timepoint.age)<=14)]
tt<-all_tts[[ct]]
tt$cell.type<-ct
tt$gene<-row.names(tt)
overlap<-tt[which((tt$P.Value<1)&(tt$gene%in%c5_markers)&(tt$cell.type==ct)),]$gene
tps<-as.vector(sapply(colnames(dge_c$normalizedExpr),function(x){as.numeric(str_split_i(x,"_",2))}))
inds<-as.vector(sapply(colnames(dge_c$normalizedExpr),function(x){str_split_i(x,"_",1)}))
dge_t<-as.data.frame(t(dge_c$normalizedExpr))
dge_t<-as.data.frame(scale(dge_t,scale = FALSE, center=TRUE))
dge_t$timepoints<-tps
dge_t$individuals<-inds
dge_t<-dge_t[which(!is.na(dge_t$timepoints)),]

dfm<-melt(dge_t,id.vars = "timepoints",measure.vars = overlap)#
df2 <- dfm %>% group_by(timepoints,variable) %>% summarise(across(value,median),.groups = 'drop') %>%as.data.frame()

data<-reshape(df2, idvar = "variable", timevar = "timepoints", direction = "wide")
if(!"value.17" %in% names(data)){
  data$value.17<-NA
}
if(!"value.16" %in% names(data)){
  data$value.16<-NA
}
if(!"value.15" %in% names(data)){
  data$value.15<-NA
}
if(!"value.2" %in% names(data)){
  data$value.2<-NA
}
data2<-data[,seq(2,17)]
data2<-sapply(data2,as.numeric)
colnames(data2)<-as.vector(sapply(colnames(data2),function(x){str_split_i(x,"value.",2)}))
row.names(data2)<-data$variable
cluster_tt_means[[ct]]<-colMeans(data2)
}
data3<-melt(do.call(cbind,cluster_tt_means))
colnames(data3)<-c("timepoints","celltype","value")
ggplot(
  data=data3, 
  aes(x=timepoints, y=value,color=celltype)) + stat_smooth(geom='line', alpha=0.25, se=FALSE, span=1) + stat_smooth(data=data3[which(data3$celltype %in% c("B_Naive","B_Mem","CD4_Naive","CD4_Mem","CD8_Naive","CD8_Mem","cDC","MAIT","gdT_Vd1","gdT_Vd2","Mono_Classical","Mono_NonClassical","NK_CD16hi")),], geom='line', se=FALSE, span=1) + theme_classic()+ theme(text = element_text(size=15)) + scale_x_continuous(limits=c(3, 14),breaks=seq(3,14,2))+ggtitle("C5: IFN response module across celltypes")+theme(legend.position = "none")
ggsave("data/output/C5_INF_AgeTrend.pdf")


# Simplified age model results
fgsea_list<-readRDS("/gpfs/gibbs/pi/csei/users/lb2336/pseudo_bulk/NICA_batch2to4_new_timepoints_limma_fgsea_simplified_age_results.rds")
fgsea_dat <- bind_rows(fgsea_list, .id = "celltype")
fgsea_dat$nlog_pval <- -log10(fgsea_dat$padj)

mean_abs_nes_dat <- fgsea_dat %>%
  group_by(pathway) %>%
  summarise(mean_abs_nes = mean(abs(NES)), mean_pval = mean(padj)) %>%
  mutate(mean_abs_nes_rank_up = rank(mean_abs_nes)) %>%
  mutate(mean_abs_nes_rank_dn = rank(-mean_abs_nes))


keep_pathways <- c(
  "reactome_Costimulation by the CD28 family",                #CD8_Mem / CD4_Mem                  2  
  "btm_S0_T cell surface signature",                          #CD4_Mem                            4
  "btm_M4.5_mitotic cell cycle in stimulated CD4 T cells",    #CD4_Mem
  "btm_M4.1_cell cycle (I)",                                  #CD8_Mem / CD4_Mem / NK_CD16hi
  "HALLMARK_G2M_CHECKPOINT",                                  #CD8_Mem / CD4_Mem / NK_CD16hi
  "btm_M4.0_cell cycle and transcription",                    #CD8_Naive / gdT_Vd1 / NK_CD16hi
  "btm_M16_TLR and inflammatory signaling",                   #CD8_Naive / Mono_Classic           1
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "reactome_Signaling by the B Cell Receptor (BCR)",
  "reactome_Fcgamma receptor (FCGR) dependent phagocytosis", 
  "reactome_Fc epsilon receptor (FCERI) signaling", 
  "reactome_FCERI mediated NF-kB activation"                  #B_Mem / B_Naive                    3
  #M4.2
  #btm_M46
  #"btm_M103_cell cycle (III)",                                #CD8_Mem / CD4_Mem / NK_CD16hi
)

fgsea_dat %>% filter(pathway %in% keep_pathways) %>%
    mutate(pathway = factor(pathway, levels = keep_pathways)) %>%
    ggplot(aes(x = celltype, y = pathway)) +
    geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
    scale_color_viridis_c() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Selected genesets that change with age - celltype specific") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("data/output/sc_age_gene_set_enrichment_selected.pdf")



#Visualize all pathways passing minor filter critera
pbulk <- dge
keep_pathways <- mean_abs_nes_dat %>%
	  	 filter(mean_pval < .25) %>%
	    	 pull(pathway)

ggplot(fgsea_dat %>% filter(pathway %in% keep_pathways), aes(x = celltype, y = pathway)) +
	geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
	scale_color_gradient2(low = "blue",mid = "white",high = "red") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	ggtitle("Genesets that change with age")

ggsave("data/output/limma_fgsea_bubble_celltype_unspecific.pdf", height = 20, width = 14)

#Visualize all pathways passing minor filter critera celltype specific
celltypes <- names(pbulk)
names(celltypes) <- celltypes
n_celltype_signif_dat <- fgsea_dat %>%group_by(pathway) %>% summarise(n_celltypes_signif = sum(padj < .05),var_NES = var(NES))

keep_pathways_celltype_specific_list <- lapply(celltypes, function(nm){
									fgsea_dat %>% filter(celltype == nm, padj < .05) %>%
									left_join(n_celltype_signif_dat) %>%
									filter(n_celltypes_signif < 5, var_NES > 1.15) %>%
									pull(pathway)})

keep_pathways_celltype_specific <- unique(unlist(keep_pathways_celltype_specific_list))

fgsea_dat %>% filter(pathway %in% keep_pathways_celltype_specific) %>%
		mutate(pathway = factor(pathway, levels = keep_pathways_celltype_specific)) %>%
		ggplot(aes(x = celltype, y = pathway)) +
			geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
			scale_color_gradient2(low = "blue",mid = "white",high = "red") +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			ggtitle("Genesets that change with age - celltype specific")

ggsave("data/output/limma_fgsea_bubble_celltype_specific.pdf", height = 20, width = 14)


#Bcells
keep_bcell_pathways <- c(
  "reactome_Signaling by the B Cell Receptor (BCR)",
  "reactome_Fcgamma receptor (FCGR) dependent phagocytosis", 
  "reactome_Fc epsilon receptor (FCERI) signaling", 
  "reactome_FCERI mediated NF−kB activation"
)

keep_bcell_pathways_ix <- which(fgsea_dat$celltype == "B_Mem" & fgsea_dat$pathway %in% keep_bcell_pathways)

keep_bcell_genes <- unique(unlist(fgsea_dat$leadingEdge[keep_bcell_pathways_ix]))

bcell_cpm <- edgeR::cpm(pbulk$B_Mem, log = TRUE)
bcell_plot_dat <- bcell_cpm[keep_bcell_genes, ] %>%
  t %>% scale %>% t %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(key = "Sample", value = "zscore_log_cpm", -gene) %>%
  left_join(pbulk$B_Mem$samples %>%  rownames_to_column("Sample"))


ggplot(bcell_plot_dat %>% filter(!grepl("IG", gene)), aes(x = matched.timepoint.age, y = zscore_log_cpm)) +
  geom_line(stat="smooth", formula = y ~ x,
            size = 1,
            #linetype ="dashed",
            aes(color = gene),
            alpha = 0.5) +
  #theme_bw() +
  theme_classic() +
  geom_smooth(col = "black", lwd = 1.5) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13)) +
  xlab("Age (Years)") +
  ylab("Scaled Expression") +
  scale_color_viridis_d(option = "mako") +
  ggtitle("BCR and FC receptor signaling in Memory B cells") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(3, 14))

ggsave("data/output/sc_FGCR_genes_age_trajectory.pdf")


#Monocyte Antigen processing and presentation
keep_mono_pathways <- c(
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
)

keep_mono_pathways_ix <- which(fgsea_dat$celltype == "Mono_Classical" & fgsea_dat$pathway %in% keep_mono_pathways)

keep_mono_genes <- unique(unlist(fgsea_dat$leadingEdge[keep_mono_pathways_ix]))

mono_cpm <- edgeR::cpm(pbulk$Mono_Classical, log = TRUE)
mono_plot_dat <- mono_cpm[keep_mono_genes, ] %>%
  t %>% scale %>% t %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(key = "Sample", value = "zscore_log_cpm", -gene) %>%
  left_join(pbulk$Mono_Classical$samples %>%  rownames_to_column("Sample"))# %>%

ggplot(mono_plot_dat %>% filter(!grepl("IG", gene)), aes(x = matched.timepoint.age, y = zscore_log_cpm)) +
  #geom_smooth(aes(color = gene), alpha = .1) +
  geom_line(stat="smooth", formula = y ~ x,
          size = 1,
          #linetype ="dashed",
          aes(color = gene),
          alpha = 0.5) +
  #theme_bw() +
  theme_classic() +
  geom_smooth(col = "black", lwd = 1.5) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13)) +
  xlab("Age (Years)") +
  ylab("Scaled Expression") +
  #scale_color_viridis_d(option = "viridis") +
  scale_color_viridis_d(option = "mako") +
  ggtitle("Antigen Presentation in Classical Monocytes") +
  coord_cartesian(ylim = c(-1, .8), xlim = c(3, 14))+theme(legend.position = "none")

ggsave("data/output/sc_AntigenPresentation_genes_age_trajectory.pdf")


#CD8 TLR & inflammatory pathwyasy
keep_ifn_pathways <- unique(fgsea_dat[str_detect(fgsea_dat$pathway,"TLR and"),]$pathway)
keep_mono_pathways_ifn <- which(fgsea_dat$celltype == "CD8_Naive" & fgsea_dat$pathway %in% keep_ifn_pathways)

keep_ifn_genes <- unique(unlist(fgsea_dat$leadingEdge[keep_mono_pathways_ifn]))

mono_ifn_cpm <- edgeR::cpm(pbulk$CD8_Naive, log = TRUE)
mono_ifn_plot_dat <- mono_ifn_cpm[keep_ifn_genes, ] %>%
  t %>% scale %>% t %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(key = "Sample", value = "zscore_log_cpm", -gene) %>%
  left_join(pbulk$CD8_Naive$samples %>%  rownames_to_column("Sample"))# %>%
# group_by(Sample, gene) %>%
# mutate(zscore_log_cpm = scale(log_cpm)) %>%
# ungroup

ggplot(mono_ifn_plot_dat, aes(x = matched.timepoint.age, y = zscore_log_cpm)) +
  geom_line(stat="smooth", formula = y ~ x,
            size = 1,
            #linetype ="dashed",
            aes(color = gene),
            alpha = 0.5) +
  #theme_bw() +
  theme_classic() +
  geom_smooth(col = "black", lwd = 1.5) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13)) +
  xlab("Age (Years)") +
  ylab("Scaled Expression") +
  scale_color_viridis_d(option = "mako") +
  ggtitle("TLR inflammatory signaling genes in CD8 Naive cells") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(3, 14))+theme(legend.position = "none")

ggsave("data/output/sc_TLR_genes_age_trajectory.pdf")


#Cell cycle pathways in CD8_Naive
keep_cc_pathways <- unique(fgsea_dat[str_detect(fgsea_dat$pathway,"btm_M4.0"),]$pathway)
keep_nk_pathways_cc <- which(fgsea_dat$celltype == "CD8_Naive" & fgsea_dat$pathway %in% keep_cc_pathways)

keep_cc_genes <- unique(unlist(fgsea_dat$leadingEdge[keep_nk_pathways_cc]))

nk_cc_cpm <- edgeR::cpm(pbulk$CD8_Naive, log = TRUE)
nk_cc_plot_dat <- nk_cc_cpm[keep_cc_genes, ] %>%
  t %>% scale %>% t %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(key = "Sample", value = "zscore_log_cpm", -gene) %>%
  left_join(pbulk$CD8_Naive$samples %>%  rownames_to_column("Sample"))

ggplot(nk_cc_plot_dat, aes(x = matched.timepoint.age, y = zscore_log_cpm)) +
  geom_line(stat="smooth", formula = y ~ x,
            size = 1,
            #linetype ="dashed",
            aes(color = gene),
            alpha = 0.5) +
  #theme_bw() +
  theme_classic() +
  geom_smooth(col = "black", lwd = 1.5) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13)) +
  xlab("Age (Years)") +
  ylab("Scaled Expression") +
  scale_color_viridis_d(option = "mako") +
  ggtitle("Cell cycle and transcription genes in CD8_Naive") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(3, 14))

ggsave("data/output/sc_CellCycle_genes_age_trajectory.pdf")
