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
cluster_gene_df<-read.csv(snakemake@input[["cluster_gene_df"]])
c1_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==1),]$Gene
c2_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==2),]$Gene
c3_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==3),]$Gene
c4_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==4),]$Gene
c5_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==5),]$Gene
c6_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==6),]$Gene
c7_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==7),]$Gene
c8_markers<-cluster_gene_df[which(cluster_gene_df$Cluster==8),]$Gene
inf_marker<-c("ADAR","IP6K2","PML","DDX58","IRF7","PSMB8","DDX60","IRF9","RSAD2","EIF2AK2","ISG15","SOCS1","EIF4A3","ISG20","SP100","HERC5","KPNA2","STAT1","IFI27","MT2A","STAT2","IFI35","MX1","TAP1","IFI6","MX2","TRIM21","IFIH1","OAS1","TRIM22","IFIT1","OAS2","TRIM38","IFIT2","OAS3","TRIM5","IFIT3","OASL","UBE2L6","IFITM1","PARP9","USP18","IFITM3","PLSCR1","XAF1")
strd_marker<-c("ABCA2","NEFL","ABCA3","NR3C2","ABHD2","PTAFR","BCL2","PTCH1","CATSPER1","PTGER2","DEFA3","RARG","EEF2","RORA","FBXO32","RORC","KLF9","RXRA","MBD4","SGK1","MBP","SPP1","NCOA4","WNT7A")

# bulk RNAseq subject-variance enriched genesets
nrchd_gene_df<-read.csv(snakemake@input[["nrchd_gene_df"]]) 


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
genesetl[["C7: Integrins, Hormones"]]<-c7_markers
genesetl[["C8: Ig heavy and light chain,\n B cells, cell cycle"]]<-c8_markers


#load variance partitioning results
vpars<-readRDS(snakemake@input[["vpars"]])
all_tts<-readRDS(snakemake@input[["all_tts"]])

#load normalized pseudobulk scRNAseq data
dge <- readRDS(snakemake@input[["dge"]]) 


### Computation of age variance gene set enrichments

#set seed
set.seed(1)


### Correlation analysis of bulk cluster 3 gene module with age across celltypes
cluster_tt_means<-list()
for (ct in names(all_tts)[which(!(names(all_tts) %in% c("T_Platelet_bind","ILC")))]){
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

data3<-data3[!is.na(data3$value),]
cor.list<-list()
for (ct in unique(data3$celltype)){
    model<-augment(loess(value ~ timepoints, data = data3[data3$celltype==ct,]),data3[data3$celltype==ct,])
    cor.list[[ct]]<-cor.test(model$timepoints,model$.fitted,method = "spearman")$p.value
    if(cor.list[[ct]]<=0.05){
        print(ct)
    }
}
print(cor.list)

### Correlation analysis of bulk cluster 5 gene module with age across celltypes
cluster_tt_means<-list()
for (ct in names(all_tts)[which(!(names(all_tts) %in% c("T_Platelet_bind","ILC")))]){
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

data3<-data3[!is.na(data3$value),]
cor.list<-list()
for (ct in unique(data3$celltype)){
    model<-augment(loess(value ~ timepoints, data = data3[data3$celltype==ct,]),data3[data3$celltype==ct,])
    cor.list[[ct]]<-cor.test(model$timepoints,model$.fitted,method = "spearman")$p.value
    if(cor.list[[ct]]<=0.05){
        print(ct)
    }
}
print(cor.list)

# Simplified age model results
sc_gene_fit_list<-readRDS(snakemake@input[["fgsea_simpl_fitl"]])
toptab_list <- lapply(sc_gene_fit_list, function(fit){
  topTable(fit, coef = "age", number = Inf) %>%
  rownames_to_column("gene")})

fgsea_list <- lapply(toptab_list, function(toptab){
  fc <- toptab$logFC
  names(fc) <- toptab$gene
  set.seed(1)
  fgseaRes <- fgsea(pathways = genesetl, stats = fc, minSize=15)
})

fgsea_dat <- bind_rows(fgsea_list, .id = "celltype")
fgsea_dat$nlog_pval <- -log10(fgsea_dat$padj)

fgsea_dat %>%
  ggplot(aes(x = celltype, y = pathway)) +
  geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_gradient2(low = "blue",mid = "white",high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Selected genesets that change with age - celltype specific") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("data/output/sc_age_gene_set_enrichment_bulk_trends.pdf",width=10)


fgsea_list<-readRDS(snakemake@input[["fgsea_simpl"]])
fgsea_dat <- bind_rows(fgsea_list, .id = "celltype")
fgsea_dat$nlog_pval <- -log10(fgsea_dat$padj)

mean_abs_nes_dat <- fgsea_dat %>%
  group_by(pathway) %>%
  summarise(mean_abs_nes = mean(abs(NES)), mean_pval = mean(padj)) %>%
  mutate(mean_abs_nes_rank_up = rank(mean_abs_nes)) %>%
  mutate(mean_abs_nes_rank_dn = rank(-mean_abs_nes))


keep_pathways <- c(
  "GO_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE",
  "reactome_Chemokine receptors bind chemokines",
  "GO_COMPLEMENT_ACTIVATION",
  "reactome_Costimulation by the CD28 family",                #CD8_Mem / CD4_Mem                  2  
  "btm_S0_T cell surface signature",                          #CD4_Mem                            4
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
    scale_shape_manual(values = c(1, 16)) +
    scale_color_gradient2(low = "blue",mid = "white",high = "red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Selected genesets that change with age - celltype specific") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("data/output/sc_age_gene_set_enrichment_selected.pdf",width=10)



#Visualize all pathways passing minor filter critera
pbulk <- dge
#keep_pathways <- mean_abs_nes_dat %>%
#	  	 filter(mean_pval < .25) %>%
#	    	 pull(pathway)

#ggplot(fgsea_dat %>% filter(pathway %in% keep_pathways), aes(x = celltype, y = pathway)) +
#	geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
#	scale_shape_manual(values = c(1, 16)) +
#	scale_color_gradient2(low = "blue",mid = "white",high = "red") +
#	theme_bw() +
#	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#	ggtitle("Genesets that change with age")

#ggsave("data/output/limma_fgsea_bubble_celltype_unspecific.pdf", height = 20, width = 16)

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
			scale_shape_manual(values = c(1, 16)) +
			scale_color_gradient2(low = "blue",mid = "white",high = "red") +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			ggtitle("Genesets that change with age - celltype specific")

ggsave("data/output/limma_fgsea_bubble_celltype_specific.pdf", height = 20, width = 16)


keep_pathways <- mean_abs_nes_dat %>%
	         filter(mean_pval < .5) %>%
		 pull(pathway)

keep_pathways_celltype_unspecific_list <- lapply(celltypes, function(nm){ 
				                                        fgsea_dat %>% filter(celltype == nm, padj < .05) %>%
				                                        left_join(n_celltype_signif_dat) %>%
                                                                        filter(n_celltypes_signif >= 5) %>%
                                                                        pull(pathway)})

keep_pathways_celltype_unspecific <- unique(unlist(keep_pathways_celltype_unspecific_list))
				    
ggplot(fgsea_dat %>% filter(pathway %in% keep_pathways) %>% filter(pathway %in% keep_pathways_celltype_unspecific), aes(x = celltype, y = pathway)) +
			geom_point(aes(size = -log10(padj), color = NES, shape = padj < .05)) +
			scale_shape_manual(values = c(1, 16)) +
			scale_color_gradient2(low = "blue",mid = "white",high = "red") +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			ggtitle("Genesets that change with age")


ggsave("data/output/limma_fgsea_bubble_celltype_unspecific.pdf", height = 20, width = 16)



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
bcell_plot_dat$matched.timepoint.age<-as.numeric(bcell_plot_dat$matched.timepoint.age)

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
mono_plot_dat$matched.timepoint.age<-as.numeric(mono_plot_dat$matched.timepoint.age)

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
  coord_cartesian(ylim = c(-1, .8), xlim = c(3, 14))

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
mono_ifn_plot_dat$matched.timepoint.age<-as.numeric(mono_ifn_plot_dat$matched.timepoint.age)


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
  coord_cartesian(ylim = c(-1, 1), xlim = c(3, 14))

ggsave("data/output/sc_TLR_genes_age_trajectory.pdf")


#CD8 complement activation pathway
keep_ifn_pathways <- unique(fgsea_dat[str_detect(fgsea_dat$pathway,"COMPLEMENT_ACTIVATION"),]$pathway)
keep_mono_pathways_ifn <- which(fgsea_dat$celltype == "CD8_TEMRA" & fgsea_dat$pathway %in% keep_ifn_pathways)

keep_ifn_genes <- unique(unlist(fgsea_dat$leadingEdge[keep_mono_pathways_ifn]))

mono_ifn_cpm <- edgeR::cpm(pbulk$CD8_TEMRA, log = TRUE)
mono_ifn_plot_dat <- mono_ifn_cpm[keep_ifn_genes, ] %>%
	  t %>% scale %>% t %>%
	    as.data.frame() %>%
	      rownames_to_column("gene") %>%
	        gather(key = "Sample", value = "zscore_log_cpm", -gene) %>%
		  left_join(pbulk$CD8_TEMRA$samples %>%  rownames_to_column("Sample"))# %>%
	  # group_by(Sample, gene) %>%
	  # mutate(zscore_log_cpm = scale(log_cpm)) %>%
	  # ungroup
	  mono_ifn_plot_dat$matched.timepoint.age<-as.numeric(mono_ifn_plot_dat$matched.timepoint.age)


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
  ggtitle("Complement activation genes in CD8_TEMRA cells") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(3, 14))

ggsave("data/output/sc_CMPL_genes_age_trajectory.pdf")



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
nk_cc_plot_dat$matched.timepoint.age<-as.numeric(nk_cc_plot_dat$matched.timepoint.age)

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
