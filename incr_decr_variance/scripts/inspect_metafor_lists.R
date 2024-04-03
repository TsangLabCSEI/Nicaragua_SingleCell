library(tidyverse)

in_dir <- "../data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/"
files <- list.files(in_dir)

files <- grep("results_dat", files, value = T)

dat_list <- lapply(files, function(f){ readRDS(file.path(in_dir, f))})
ids <- sapply(strsplit(files, "__"), `[[`, 2)
ids <- gsub(".rds", "", ids)
dat_list <- lapply(dat_list, `[[`, 1)
names(dat_list) <- ids

combined_dat <- bind_rows(dat_list, .id = "celltype")

combined_dat <- combined_dat %>%
        filter(term == "window_number_demeaned") %>%
        group_by(celltype) %>%
        mutate(padj = p.adjust(pval, method = "fdr"))

n_signif_dat <- combined_dat %>% 
        summarise(
                  n_p05_signif = sum(pval < .05),
                  n_p20_signif = sum(pval < .2),
                  n_fdr05_signif = sum(padj < .05),

        )
n_signif_dat

n_signif_dat %>% write_csv("../data/analysis_out/n_signif_dat_sliding_window.csv")

combined_dat_filtered_signif <- combined_dat %>% filter(pval < .05)
