library(fgsea)
library(dplyr)
library(stringr)

set.seed(1)

vpars<-readRDS(snakemake@input[["vpars"]]) #"data/output/NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_vpar.rds") #snakemake@input[["vpars"]])
nrchd_gene_df<-read.csv(snakemake@input[["nrchd_genes"]]) #"data/input/paper_supplementary_tables/supp_tables_S10.csv") #snakemake@input[["nrchd_genes"]]) #leading edge genes

#aggregate VES
for (ct in names(vpars)){
  vpars[[ct]]$age_variance_explained<-rowSums(vpars[[ct]][names(vpars[[ct]])[startsWith(names(vpars[[ct]]),"X")]])
  vpars[[ct]]$celltype<-ct
  vpars[[ct]]$gene<-sapply(rownames(vpars[[ct]]),function(x){str_split_i(x,"\\.",1)})
}
vpars<-do.call(rbind,vpars)[c("batch","matched.individual","gender","age_variance_explained","Residuals","celltype","gene")]
colnames(vpars)<-c("batch_variance_explained","subject_variance_explained","sex_variance_explained","age_variance_explained","residual_variance_explained","celltype","gene")
saveRDS(vpars,"data/output/NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_vpars_object.rds")

#retrieve VES enriched genes from bulk
genesetl<-list()
for (i in seq(1,dim(nrchd_gene_df)[1])){
  genesetl[[nrchd_gene_df[i,"pathway"]]]<-unlist(str_split(nrchd_gene_df[i,"leadingEdge"]," "))
}

#double enrich VES enriched genes in single cell
enrichment_list<-list()
for (ct in unique(vpars$celltype)){
  vpars_ct<-vpars[which(vpars$celltype==ct),]
  ranks <- vpars_ct$subject_variance_explained
  names(ranks)<-vpars_ct$gene
  ranks <- sort(ranks, decreasing = T)

  enrichment_list[[ct]]<-fgsea(pathways = genesetl, ranks)
  enrichment_list[[ct]]<-enrichment_list[[ct]][enrichment_list[[ct]]$pval<0.05]
  enrichment_list[[ct]]$celltype<-ct
}

enrichment_df<-do.call(rbind,enrichment_list)
enrichment_df$nlog_pval<-(-log10(enrichment_df$pval))

saveRDS(enrichment_df,"data/output/NICA_batch2to4_new_timepoints_fsex_subject_enriched_with_labels.rds")

#renorm vpars by taking out batch variance

sc_varpart <- vpars %>% mutate(var_bio_sum = subject_variance_explained + sex_variance_explained + age_variance_explained + residual_variance_explained) %>% mutate(ves_renorm = subject_variance_explained / var_bio_sum, vex_renorm= sex_variance_explained / var_bio_sum, vea_renorm= age_variance_explained / var_bio_sum, res_renorm=residual_variance_explained / var_bio_sum)
sc_varpart<-sc_varpart[c("ves_renorm","vex_renorm","vea_renorm","res_renorm","celltype","gene")]
colnames(sc_varpart)<-c("subject_variance_explained","sex_variance_explained","age_variance_explained","residual_variance_explained","celltype","gene")
saveRDS(sc_varpart,"data/output/NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_vpars_object_renorm.rds")
