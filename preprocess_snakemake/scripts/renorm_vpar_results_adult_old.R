library(dplyr)

set.seed(4)

nind_list<-readRDS(snakemake@input[["vpars"]])

#renorm vpars by taking out batch variance
nind_list_new<-list()
for (ct in names(nind_list)){
  nind_list_new.idx<-list()
  pseudo.sce.list<-nind_list[[ct]]
  for (idx in seq(1,20)){
    pseudo.sce.list.idx<-pseudo.sce.list[[idx]][c("Batch","Donor_id","Age","Residuals","cell.type")]
    colnames(pseudo.sce.list.idx)<-c("batch","subject_variance_explained","age_variance_explained","residual_variance_explained","celltype")
    sc_varpart <- pseudo.sce.list.idx %>% mutate(var_bio_sum = subject_variance_explained + age_variance_explained + residual_variance_explained) %>% mutate(ves_renorm = subject_variance_explained / var_bio_sum, vea_renorm= age_variance_explained / var_bio_sum, res_renorm=residual_variance_explained / var_bio_sum)
    sc_varpart$gene <- rownames(pseudo.sce.list.idx)
    sc_varpart<-sc_varpart[c("ves_renorm","vea_renorm","res_renorm","celltype","gene")]
    colnames(sc_varpart)<-c("subject_variance_explained","age_variance_explained","residual_variance_explained","celltype","gene")
    nind_list_new.idx[[idx]]<-sc_varpart
  }
  nind_list_new[[ct]]<-nind_list_new.idx
}


saveRDS(nind_list_new,"data/output/age_subject_variance_VESlist_adult_old_aging_cohort_renorm.RDS")
