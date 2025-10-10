# Nicaragua_SingleCell
This repository contains notebooks to reproduce the single cell analysis of the Nicaragua project - revealing baseline immune signatures in blood transcriptomes of children.

## Data Availability

Processed data files to reproduce the analysis are available on zenodo under: \
Longitudinal blood transcriptomic profiling of Nicaraguan children reveals trajectories of pediatric immune development and individuality \
https://zenodo.org/api/records/15151094

## Data Analysis Instructions for Reproducibility

Below are pointers to scripts to reproduce figures and respective pre-processing steps. Each analysis below the following captions can be reproduced independently based on provided intermediate results in the zenodo repository. Raw data processing can not be reproduced as we are unable to share this data, nevertheless we share the scripts how the analysis was carried out.

### Environment and Container
To reproduce the analysis conda environment can be created from one of the yml files in the environment folder or the following singularity container available on https://cloud.sylabs.io/ can be used: \
```library://leon.bichmann/nicaragua/nicaragua:latest```

In the case of choosing one of the provided conda environments, also install the local Rpackage ImmuneAgeStability: \
```install.packages("Nicaragua_SingleCell/ImmuneAgeStability", repos = NULL, type="source")``` \

### ImmuneAgeStability Rpackage
This Rpackage provides a way to reproduce some of the analysis carried out in our study as described in the scripts below. In addition it also enables users to assess any of their own gene expression signatures for aspects of temporal stabilty during childhood. 

## All-in-one snakemake workflows 

### Snakemake workflow to reproduce all figures at once
If you are familiar with snakemake and simply want to recreate all figures in one step, just run the follwing: \
-> **GoTo:** figures_snakemake \
-> **Data:** download all data from the zenodo repository into figures_sakemake/data/input/ \
-> **Execution:**  \
```snakemake --use-singularity```

### Snakemake workflow to reproduce all preprocessing steps at once
If you are familiar with snakemake and simply want to recreate all preprocessing in one step apart from the individuality, just run the follwing: \
-> **GoTo:** preprocess_snakemake \
-> **Data:** download all data from the zenodo repository into preprocess_sakemake/data/input/ \
-> **Execution:**  \
```snakemake --use-singularity```

### Snakemake workflow to compute sliding window increasing and decreasing VES
To re-compute the individuality scores - run the following separate snakemake workflow: \
-> **Data:** Download celltype divided pseudobulk RDS objects from Zenodo: "pseudobulk_per_celltype.zip" \
-> **Data:** Download singularity image from Zenodo: "intrinsicness_0.1.sif" \
-> **Prepare:** Extract celltype specific pseudobulk folder into incr_decr_variance/data/processed/ \
-> **Code to run:** incr_decr_variance/Snakefile \
-> **Example execution:** \
```snakemake -j 16 --configfile config.yaml --use-singularity --executor slurm```

## Additional details on pre-processing not covered in snakemake workflows

### Pre-processing: Pseudobulk Generation in Adult Old Cohort (Terekhova et al, 2023)
We separated this processing step out since the data was very large and needed high memory to process. It can be done following these steps: \
-> **Data:** Download public data set from synapse.org (syn49637038): "all_pbmcs_rna_harmony.h5ad" \
-> **Data:** Download public data set from synapse.org (syn49637038): "pbmc_gex_raw_with_var_obs.h5ad" \
-> **Data:** Download bootstrapped donor list from Zenodo: "age_subject_variance_donors_adult_old_aging_cohort.RDS" \
-> **Code to run:** adult_old_cohort/export_pseudobulk_adult_old.R

### Raw data processing: Cellranger, Normalization, Demuxlet
Due to privacy concerns we did not make the raw data available but we outlined the steps how it was done here: \
-> **Data:** HTO mapping of timepoints: "HTO_matching_table_new_cleaned.xlsx" \
-> **Code to run:** Configurations how cellranger was run can be found here - singlecell_processing/example_cellranger.config \
-> **Code to run:** Example code how the cellranger output was normalized and annotated can be found here - singlecell_processing/sc_processing.R
