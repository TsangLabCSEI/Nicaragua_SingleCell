# Nicaragua_SingleCell
This repository contains notebooks to reproduce the single cell analysis of the Nicaragua project - revealing baseline immune signatures in blood transcriptomes of children.

![](analysis_overview.png)

## Data Analysis Instructions for Reproducibility

### Figure 4B,C - Single Cell UMAP and Age Range of Individuals
-> **Data:** Download single cell RDS object from Zenodo: "NICAall_combined_manual_labeled_cleaned_final.rds" \
-> **Code to run:** downstream_analysis/Individuals_and_UMAP.Rmd 

### Pre-processing: Pseudobulk Generation from Single Cell Data
-> **Data:** Download single cell RDS object from Zenodo: "NICAall_combined_manual_labeled_cleaned_final.rds" \
-> **Code to run:** pseudobulk_processing/export_pseudobulk.R

### Pre-processing: Pseudobulk Normalization
-> **Data:** Download pseubulk RDS object from Zenodo: "NICA_batch2to4_new_timepoints_cleaned.rds" \
-> **Code to run:** pseudobulk_processing/normalize_and_batch_correct_pseudobulk.R

### Pre-processing: Age and Subject Variance Partitioning
-> **Data:** Download normalized pseudobulk RDS object from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized.rds" \
-> **Code to run:** pseudobulk_processing/fit_linear_models_pseudobulk.Rmd

### Figure 4D - Age Trajectories
-> **Data:** bulk age trajectory gene sets from Zenodo: "supp_tables_S6.csv" \
-> **Data:** bulk subject variance gene sets from Zenodo: "supp_tables_S10.csv" \
-> **Data:** raw single cell variance partition results from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_spline.rds" \
-> **Data:** aggregated single cell variance partition results from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_vpars_object.rds" \
-> **Data:** Download normalized pseudobulk RDS object from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized.rds" \
-> **Data:** Download simplified age variance partition results: "/Users/leon/Documents/NICA/sc_analysis/data/NICA_batch2to4_new_timepoints_limma_fgsea_simplified_age_results.rds" \
-> **Code to run:** downstream_analysis/Age_trends_and_enrichments.Rmd


