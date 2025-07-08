library(tidyverse)
library(dplyr)
library(metafor)
library(SummarizedExperiment)

options(stringsAsFactors = FALSE)

source("scripts/util/jaccard_matrix.R")

CELLTYPE <- snakemake@params[["celltype"]]

IN_DIR <- "data/analysis_out/variancePartition/age_sliding_window_bootstrap/varpart_objects/"

files <- list.files(IN_DIR, full.names = TRUE, pattern=CELLTYPE)

# META_DIR <- "data/analysis_out/variancePartition/age_sliding_window/metadata_subsets/"
# 
# meta_files <- list.files(META_DIR, full.names = TRUE)
# 
# meta_objs <- lapply(meta_files, read_tsv)

META_DIR <- paste0("data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells_",CELLTYPE)
meta_files <- list.files(META_DIR, full.names = TRUE)

vst_objs <- lapply(meta_files, readRDS)
meta_objs <- lapply(vst_objs, function(x){x$samples})

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


vp_objs <- lapply(files, readRDS)
id <- sapply(strsplit(basename(files), "_"), function(x) paste(x[2:5], collapse = "_"))
names(vp_objs) <- id

dat <- lapply(vp_objs, function(vp){
  data.frame(vp) %>%
  rownames_to_column("gene")
}) %>% bind_rows(.id = "id")

# dat <- dat %>%
#         separate(id, into =c("window", "boot_iter"), sep = "_")

id_split <- strsplit(dat$id, "_")
dat$window <- sapply(id_split, `[[`, 1)

id_split <- strsplit(dat$id, "boot-")
dat$boot_iter <- sapply(id_split, `[[`, 2)

dat <- dat %>%
        mutate(window_number = as.numeric(sapply(strsplit(window, "to"), `[[`, 1)))

head(dat %>% dplyr::filter(gene =="DDX11L1"))
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
head(dat_renorm %>% dplyr::filter(gene =="DDX11L1"))

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


SUMM_OUT_PATH <- paste0("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/bootstrap_summary_dat_list_",CELLTYPE,".rds")
dir.create(dirname(SUMM_OUT_PATH), recursive = TRUE)
saveRDS(summ_list, SUMM_OUT_PATH)
#summ_list <- readRDS(SUMM_OUT_PATH)


#summ_single <- summ %>% filter(gene == "STAT1") %>%
#        arrange(window_number) %>%
#        mutate(yi = mean_vexp) %>%
#        mutate(window_number_demeaned = window_number - median(.$window_number))
#
#
#V <- combined_jacc_mat * sqrt(outer(summ_single$var_vexp, summ_single$var_vexp))
#
##just checking that I made my covariance matrix correctly
#cov2cor(V) - combined_jacc_mat < 1e-15
#
#mod <- rma.mv(yi ~ window_number_demeaned, V, data = summ_single) 

res_list <- lapply(summ_list["subj_div_subj_age_gender_resid"], function(summ){
  print("starting next")
  genes <- unique(summ$gene)
  names(genes) <- genes

  res_list <- lapply(genes, function(gene_name){
    gene_index <- which(genes == gene_name)
    if(gene_index %% 500 == 0){
      print(gene_index)
    }
    summ_single <- summ %>% dplyr::filter(gene == gene_name) %>%
            arrange(window_number) %>%
            mutate(yi = mean_vexp) %>%
            mutate(window_number_demeaned = window_number - median(.$window_number))
    
    V <- combined_jacc_mat * as.vector(sqrt(outer(summ_single$var_vexp, summ_single$var_vexp)))
    
    tryCatch({
    mod <- rma.mv(yi ~ window_number_demeaned + n_subj, V, data = summ_single) 
  
    data.frame(term = rownames(mod$beta), beta = mod$beta[,1], se = mod$se, zval = mod$zval, pval = mod$pval)
  
    }, error = function(e) NULL)
  })
  res_dat_slopeint <- bind_rows(res_list, .id = "gene")

})

DAT_FULL_OUT_PATH <- paste0("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/results_dat_list_with_intercept_",CELLTYPE,".rds")

saveRDS(res_list, DAT_FULL_OUT_PATH)

#res_list$sex_div_subj_age_sex_resid %>% 
#        filter(term == "window_number_demeaned") %>%
#        mutate(padj = p.adjust(pval)) %>%
#        arrange(pval) %>%
#        head(20)
#
#
#res_list$sex_div_sex_resid %>% 
#        filter(term == "window_number_demeaned") %>%
#        mutate(padj = p.adjust(pval)) %>%
#        arrange(pval) %>%
#        head(20)
#
#res_list$age_div_subj_age_gender_resid %>% 
#        filter(term == "window_number_demeaned") %>%
#        mutate(padj = p.adjust(pval)) %>%
#        arrange(pval) %>%
#        head(20)
#
#res_list$subj_div_subj_age_gender_resid %>% 
#        filter(term == "window_number_demeaned") %>%
#        mutate(padj = p.adjust(pval)) %>%
#        arrange(pval) %>%
#        head(20)

