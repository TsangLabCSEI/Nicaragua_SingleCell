library(SummarizedExperiment)
library(limma)

#change input file based on celltype
vst_IN_PATH <- snakemake@input[[1]] 

#change OUT paths
OUT_DIR <- dirname(snakemake@output[[1]]) 
dir.create(OUT_DIR, recursive = TRUE)

vst <- readRDS(vst_IN_PATH)

#vst <- vst[, vst$Sample.type == "case"]
meta <- as.data.frame(vst$samples)

min_age_mo <- min(meta$matched.timepoint.age)
max_age_mo <- max(meta$matched.timepoint.age)


window_size <- 3 # 3 years
age_windows_start <- seq(3, 14-window_size)
age_window_end <- age_windows_start + window_size

window_dat <- data.frame(window_name = paste(age_windows_start, age_window_end, sep = "to"),
                         start_age = age_windows_start,
                         end_age = age_window_end
)

vst_window_list <- lapply(1:nrow(window_dat), function(i){
  min_age <- window_dat$start_age[i]
  max_age <- window_dat$end_age[i]
  vst[, meta$matched.timepoint.age >= min_age & meta$matched.timepoint.age <= max_age]
})

lapply(vst_window_list, dim)
lapply(vst_window_list, function(x){table(table(x$samples$matched.individual))})

for(i in seq_along(vst_window_list)){
  window_name <- window_dat$window_name[[i]]
  out_path <- file.path(OUT_DIR, paste0("age_", window_name, "_vst.rds"))
  print(out_path)
  obj <- vst_window_list[[i]]
  saveRDS(obj, out_path)
}
