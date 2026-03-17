# set seed
set.seed(4)

# load libraries
library(ImmuneAgeStability)
library(pheatmap)
library(stringr)
library(patchwork)
library(viridis)
library(circlize)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)


# load data
pbulk_data<-read.csv(snakemake@input[["data"]])

#get background genes
background_genes_list<-list()
for (ct in names(pbulk_data)){
  print(ct)
  emat <- pbulk_data[[ct]]$counts
  
  gene_means <- apply(emat, 1, mean)
  gene_sds <- apply(emat, 1, sd)
  
  gene_vars <- gene_sds
  
  n <- 5
  
  stopifnot(identical(names(gene_means), names(gene_sds)))
  mat <- cbind(scale(gene_sds), scale(gene_means))
  #rownames(mat) <- names(gene_vars)
  
  dist_mat <- as.matrix(dist(mat))
  
  background_genes <- lapply(all_sigs, function(sig_genes){
    print("starting")
    names(sig_genes) <- sig_genes
    topn_closest <- lapply(sig_genes, function(gene){
      if(!gene %in% rownames(dist_mat)){
        return(rep(NA, n))
      }
      names(sort(dist_mat[gene, setdiff(colnames(dist_mat), gene)]))[1:n]
    })
  })
  
  background_unlist <- unique(unlist(background_genes, use.names = F))
  
  background_unlist <- setdiff(background_unlist, unlist(all_sigs_sub))
  
  background_genes_list[[ct]]<-background_unlist
}

background_genes_list$bulk<-ImmuneAgeStability::baseline_signatures$background_bulk

#visualize subject variance partitioning for selected immune response signatures
sigs<-ImmuneAgeStability::baseline_signatures[c("IHM","IFN","ia_bcell","iaa","GPR56highvslow_CD8_EM_sig","DE_CD29hi_CD8_genes")]
names(sigs)<-c("IHM","IFN","IAB B-cell","IAA","CD8 VM GPR56hi","CD8 VM CD29hi")
multisignature_stability_df<-ImmuneAgeStability::compute_multisignature_stability(sigs)

background_stability_df<-ImmuneAgeStability::compute_multisignature_stability(background_genes_list)
background_stability_df<-background_stability_df[background_stability_df$celltype==background_stability_df$signature,]
background_stability_df$signature<-"Background"

stability_df<-rbind(background_stability_df,multisignature_stability_df)
pdf("data/output/signature_stability_heatmap.pdf")
ImmuneAgeStability::plot_signature_stability_heatmap(stability_df)
dev.off()

#visualize subject variance time trends for selected immune response signature

stability_slope_df<-ImmuneAgeStability::compute_signature_stability_slope(signature_name = "IAB-Bcell", signature = ImmuneAgeStability::baseline_signatures$ia_bcell)
pdf("data/output/IAB-Bcell-stability_trend.pdf")
ImmuneAgeStability::plot_signature_stability_bars(stability_slope_df)
dev.off()

stability_slope_df<-ImmuneAgeStability::compute_signature_stability_slope(signature_name = "IFN", signature = ImmuneAgeStability::baseline_signatures$IFN)
pdf("data/output/IFN-stability_trend.pdf")
ImmuneAgeStability::plot_signature_stability_bars(stability_slope_df)
dev.off()

stability_slope_df<-ImmuneAgeStability::compute_signature_stability_slope(signature_name = "CD8 VM GPR56hi", signature = ImmuneAgeStability::baseline_signatures$GPR56highvslow_CD8_EM_sig)
pdf("data/output/CD8_VM_GPR56hi-stability_trend.pdf")
ImmuneAgeStability::plot_signature_stability_bars(stability_slope_df)
dev.off()

background_gene_VES<-list()
for(ct in names(background_genes_list)){
  print(ct)
  background_gene_VES[[ct]]<-compute_signature_stability_slope(signature_name = ct, signature = background_genes_list[[ct]])
}
background_gene_VES_df<-do.call(rbind,background_gene_VES)
background_gene_VES_df<-background_gene_VES_df[background_gene_VES_df$celltype==background_gene_VES_df$signature,]
background_gene_VES_df$signature<-"Background"
select_celltypes<-c("bulk","B_Mem","B_Naive","CD4_Mem","CD4_Naive","CD8_Mem","CD8_Naive","NK_CD16hi","Mono_Classical","Mono_NonClassical","gdT_Vd1","gdT_Vd2","MAIT","pDC")
background_gene_VES_df <- background_gene_VES_df %>% filter(celltype %in% select_celltypes) %>% mutate(celltype = factor(celltype, levels = select_celltypes)) %>% arrange(celltype)
pdf("data/output/Background-stability_trend.pdf")
plot_signature_stability_bars(background_gene_VES_df)
dev.off()


#visualize explicit subject variance time trend of IFN signature for monocytes
stability_rateOfchange_df<-compute_signature_stability_rateofchange(signature_name = "IFN", signature = baseline_signatures$IFN)
pdf("data/output/IFN-stability_trend_explicit_monocytes.pdf")
plot_signature_agetrend(signature_rateOfchange = stability_rateOfchange_df,cell_subsets = c("Mono_Classical","Mono_NonClassical"))
dev.off()

stability_rateOfchange_df<-compute_signature_stability_rateofchange(signature_name = "IFN - genes", signature = age_trajectories$`IFN genes`)
pdf("data/output/IFN-stability_trend_explicit_monocytes_genes.pdf")
plot_signature_gene_agetrend(signature_rateOfchange = stability_rateOfchange_df,cell_subsets = c("Mono_Classical"))
dev.off()
