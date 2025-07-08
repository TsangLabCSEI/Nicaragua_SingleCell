# Nicaragua_SingleCell
This repository contains notebooks to reproduce the single cell analysis of the Nicaragua project - revealing baseline immune signatures in blood transcriptomes of children.

![](analysis_overview.png)

## Data Availability

Processed data files to reproduce the analysis are available on zenodo under: \
Longitudinal blood transcriptomic profiling of Nicaraguan children reveals trajectories of pediatric immune development and individuality \
https://zenodo.org/api/records/15151094

## Data Analysis Instructions for Reproducibility

### Pre-processing: Pseudobulk Generation from Single Cell Data
-> **Data:** Download single cell RDS object from Zenodo: "NICAall_combined_manual_labeled_cleaned_final.rds" \
-> **Code to run:** pseudobulk_processing/export_pseudobulk.R

### Pre-processing: Pseudobulk Normalization
-> **Data:** Download pseubulk RDS object from Zenodo: "NICA_batch2to4_new_timepoints_cleaned.rds" \
-> **Code to run:** pseudobulk_processing/normalize_and_batch_correct_pseudobulk.R

### Pre-processing: Age and Subject Variance Partitioning
-> **Data:** Download normalized pseudobulk RDS object from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized.rds" \
-> **Code to run:** pseudobulk_processing/fit_linear_models_pseudobulk.Rmd

### Pre-processing: Increasing and Decreasing Subject Individuality
-> **Data:** Download celltype divided pseudobulk RDS objects from Zenodo: "pseudobulk_per_celltype.zip" \
-> **Data:** Download singularity image from Zenodo: "intrinsicness_0.1.sif" \
-> **Prepare:** Extract celltype specific pseudobulk folder into incr_decr_variance/data/processed/ \
-> **Code to run:** incr_decr_variance/Snakefile \
-> **Example execution:**
```snakemake -j 16 --configfile config.yaml --use-singularity --executor slurm```

### Figure 4B,C - Single Cell UMAP and Age Range of Individuals
-> **Data:** Download single cell RDS object from Zenodo: "NICAall_combined_manual_labeled_cleaned_final.rds" \
-> **Code to run:** downstream_analysis/Individuals_and_UMAP.Rmd

### Figure 4D - Age Trajectories
-> **Data:** bulk age trajectory gene sets from Zenodo: "supp_tables_S6.csv" \
-> **Data:** bulk subject variance gene sets from Zenodo: "supp_tables_S10.csv" \
-> **Data:** raw single cell variance partition results from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_spline.rds" \
-> **Data:** aggregated single cell variance partition results from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized_fsex_vpars_object.rds" \
-> **Data:** Download normalized pseudobulk RDS object from Zenodo: "NICA_batch2to4_new_timepoints_cleaned_normalized.rds" \
-> **Data:** Download simplified age variance partition results: "/Users/leon/Documents/NICA/sc_analysis/data/NICA_batch2to4_new_timepoints_limma_fgsea_simplified_age_results.rds" \
-> **Code to run:** downstream_analysis/Age_trends_and_enrichments.Rmd

### Raw data processing: Cellranger, Normalization, Demuxlet
-> **Data:** Due to privacy concerns we did not make the raw data available \
-> **Data:** HTO mapping of timepoints: "HTO_matching_table_new_cleaned.xlsx" \
-> **Code to run:** Configurations how cellranger was run can be found here - singlecell_processing/example_cellranger.config \
-> **Code to run:** Example code how the cellranger output was normalized and annotated can be found here - singlecell_processing/sc_processing.R
