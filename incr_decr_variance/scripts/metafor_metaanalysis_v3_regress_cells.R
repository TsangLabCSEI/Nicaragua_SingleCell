library(tidyverse)
library(metafor)
library(SummarizedExperiment)

options(stringsAsFactors = FALSE)

source("scripts/util/jaccard_matrix.R")

CELLTYPE <- snakemake@params[["celltype"]]

#IN_DIR <- "data/analysis_out/variancePartition/age_sliding_window_bootstrap/varpart_objects/"
#IN_DIR <- snakemake@input[[1]]
#files <- list.files(IN_DIR, full.names = TRUE)
#print(str(snakemake@input))
#print(snakemake@input)
files <- unlist(snakemake@input)


SUMM_OUT_PATH <- snakemake@output[[1]]
DAT_FULL_OUT_PATH <- snakemake@output[[2]]

print(1)
# META_DIR <- "data/analysis_out/variancePartition/age_sliding_window/metadata_subsets/"
# 
# meta_files <- list.files(META_DIR, full.names = TRUE)
# 
# meta_objs <- lapply(meta_files, read_tsv)
vp_files <- grep("varpart.rds", files, value = T)
meta_files <- grep("vst.rds", files, value = T)

#META_DIR <- unique(dirname(files))
#META_DIR <- gsub("varpart_objects", "metadata_subsets", META_DIR)
#stopifnot(length(META_DIR) == 1)
#meta_files <- list.files(META_DIR, full.names = TRUE)
#meta_files <- grep("meta", meta_files, value = T)
#print(str(meta_files))

print("36")
dgelist_objs <- lapply(meta_files, readRDS)
print("37")
meta_objs <- lapply(dgelist_objs, `[[`, "samples")
print(str(meta_objs))

names(meta_objs) <- sapply(strsplit(basename(meta_files), "\\."), `[[`, 1)
names(meta_objs) <- gsub("age_", "", names(meta_objs))
names(meta_objs) <- gsub("_vst", "", names(meta_objs))

meta_obj_ages <- sapply(strsplit(names(meta_objs), "to"), `[[`, 1)
meta_obj_order <- order(as.numeric(meta_obj_ages))

meta_objs <- meta_objs[meta_obj_order]

n_samp_dat <- sapply(meta_objs, nrow) %>% enframe(name = "window", value = "n_samp")

n_subj_dat <- sapply(meta_objs, function(x){length(unique(x$matched.individual))}) %>% 
        enframe(name = "window", value = "n_subj")


subj_list <- lapply(meta_objs, `[[`, "matched.individual")
samp_list <- lapply(meta_objs, rownames)

subj_jacc_mat <- jaccardMat(subj_list)
samp_jacc_mat <- jaccardMat(samp_list)

#I will use this combined jaccard similarity as an estimate of the correlation that should exist
combined_jacc_mat <- (subj_jacc_mat + samp_jacc_mat) / 2

pdf(snakemake@output[[3]])
pheatmap::pheatmap(combined_jacc_mat, cluster_rows = F, cluster_cols = F)
dev.off()
write_tsv(data.frame(combined_jacc_mat), snakemake@output[[4]])

vp_objs <- lapply(vp_files, readRDS)
id <- sapply(strsplit(basename(vp_files), "_"), function(x) paste(x[2:3], collapse = "_"))
names(vp_objs) <- id

dat <- lapply(vp_objs, function(vp){
  data.frame(vp) %>%
  rownames_to_column("gene")
}) %>% bind_rows(.id = "id")

# dat <- dat %>%
#         separate(id, into =c("window", "boot_iter"), sep = "_")

id_split <- strsplit(dat$id, "_")

dat$window <- sapply(id_split, `[[`, 1)
dat$boot_iter <- sapply(id_split, `[[`, 2)

dat <- dat %>%
        mutate(window_number = as.numeric(sapply(strsplit(window, "to"), `[[`, 1)))

head(dat %>% filter(gene =="DDX11L1"))
#head(dat)

dat_renorm <- with(dat, data.frame(subj_div_all = Subject.ID,
                        subj_div_subj_age_gender_resid = Subject.ID / (Subject.ID + Age.months + sex.numeric + Residuals), 
                        subj_div_subj_resid = Subject.ID / (Subject.ID + Residuals),
                        age_div_subj_age_gender_resid = Age.months / (Subject.ID + Age.months + sex.numeric + Residuals),
                        age_div_age_resid = Age.months / (Age.months + Residuals),
                        sex_div_subj_age_sex_resid =  sex.numeric / (Subject.ID + Age.months + sex.numeric + Residuals),
                        sex_div_sex_resid =  sex.numeric / (sex.numeric + Residuals),
                        window = window, 
                        window_number = window_number,
                        gene = gene
                        ))
head(dat_renorm %>% filter(gene =="DDX11L1"))

keepcols <- setdiff(colnames(dat_renorm), c("window", "window_number", "gene"))
names(keepcols) <- keepcols

summ_list <- lapply(keepcols, function(nm){
  print(nm)
  dat_tmp <- data.frame(x = dat_renorm[[nm]], window_number = dat_renorm$window_number, window = dat_renorm$window, gene = dat_renorm$gene)
  dat_tmp %>% group_by(window, gene) %>%
        summarise(med_vexp = median(x), mean_vexp = mean(x),
                  var_vexp = var(x), window_number = unique(window_number)
        )

})
names(summ_list) <- keepcols

summ_list <- lapply(summ_list, function(x){
  x %>% left_join(n_subj_dat)
})


#SUMM_OUT_PATH <- paste0("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/bootstrap_summary_dat_list_",CELLTYPE,".rds")
dir.create(dirname(SUMM_OUT_PATH))
saveRDS(summ_list, SUMM_OUT_PATH)
#summ_list <- readRDS(SUMM_OUT_PATH)


res_list <- lapply(summ_list["subj_div_subj_age_gender_resid"], function(summ){
  print("starting next")
  genes <- unique(summ$gene)
  names(genes) <- genes

  res_list <- lapply(genes, function(gene_name){
    gene_index <- which(genes == gene_name)
    if(gene_index %% 500 == 0){
      print(gene_index)
    }
    summ_single <- summ %>% filter(gene == gene_name) %>%
            arrange(window_number) %>%
            mutate(yi = mean_vexp) %>%
            mutate(window_number_demeaned = window_number - median(.$window_number))
    
    print("Variance mat, all variance put in matrix form to be scaled by jaccmat")
    print(str(sqrt(outer(summ_single$var_vexp, summ_single$var_vexp))))
    print("combined_jacc_mat")
    print(str(combined_jacc_mat))
    V <- combined_jacc_mat * sqrt(outer(summ_single$var_vexp, summ_single$var_vexp))
    
    tryCatch({
    mod <- rma.mv(yi ~ window_number_demeaned + n_subj, V, data = summ_single) 
  
    data.frame(term = rownames(mod$beta), beta = mod$beta[,1], se = mod$se, zval = mod$zval, pval = mod$pval)
  
    }, error = function(e) NULL)
  })
  res_dat_slopeint <- bind_rows(res_list, .id = "gene")

})

#DAT_FULL_OUT_PATH <- paste0("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/results_dat_list_with_intercept_",CELLTYPE,".rds")

saveRDS(res_list, DAT_FULL_OUT_PATH)
