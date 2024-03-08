singularity: "file:///gpfs/gibbs/pi/csei/users/lb2336/incr_decr_svariance/singularity/intrinsicness_0.1.sif" 

##### alignment and quantification
#### bash ../data/raw/nf-core_rnaseq/nf_call
#### bash ../data/raw/nf-core_rnaseq_pilot/nf_call
#Command: can put this in the workflow 
#nextflow run nf-core/rnaseq -resume -profile conda --reads '/hpcdata/sg/sg_data/illumina_NCI_runs/Nicaragua/{JohnTsang_CS026456_4libpools_121319_HNVMLDSXX,JohnTsang_CS026580_4libpools_01132020_HW3HJDSXX}/*/*_R{1,2}_001.fastq.gz' --genome GRCh37 --saveAlignedIntermediates

###
##### variant calling
####bash ../variant_calling/nf-core_rnavar/nf_call
##### Imputation prep
####data/analysis_out/variants/imputation/prep_for_imputation/
#### eqtl analysis
#### processing
#### post imputation
###Rscript variants/michigan_imputation/0_merge_vcfs/merge_chromosomes.sh
###Rscript variants/michigan_imputation/convert_vcf.R

##### seq2hla
####bash seq2hla/run_seq2hla/sm_call

#inputs

#If starting from scratch, make sure that have 
#data/raw/nf-core_rnaseq/results/featureCounts/merged_gene_counts.txt
#data/raw/nf-core_rnaseq_pilot/results/featureCounts/merged_gene_counts.txt

## metadata - this can be skipped and assumed to be given
# output: data/metadata/meta_all.rds
# might want to update this to have version without extraneous things that people don't need
#Rscript processing/metadata/compile_metadata_pilot.R
#Rscript processing/metadata/compile_metadata_pilot.R
#Rscript processing/metadata/combine_pilot_full_add_clinical.R

CELLTYPES = ["B_Mem","B_Naive","CD4_Mem","CD4_Naive","CD8_Mem","CD8_Naive","cDC","gdT_Vd1","gdT_Vd2","MAIT","Mono_Classical","Mono_NonClassical","NK_CD16hi","NK_CD56hi","pDC"]


rule all:    
    input:
        expand("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/fgsea/slope_{celltype}.tsv", celltype=CELLTYPES),
        expand("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/bootstrap_summary_dat_list_{celltype}.rds", celltype=CELLTYPES)

## processsing -------------------------------------------------------

#sliding window analysis --------------------------------------

WINDOWS = []
AGE_GROUPS = [1,2,3]

for i in range(1,13):
    WINDOWS.append(str(i) + "to" + str(i + 3))

BOOT_ITER = range(0,100)

rule make_window_subsets_rm_cells:
    input:
        "data/processed/NICA_batch1to4_new_normalized_{celltype}.rds" ### not sure what is this
    output:
        expand("data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells_{celltype}/age_{window}_vst.rds", window = WINDOWS)
    script:
        "scripts/subset_data_ctrl_cells.R"

rule run_varpart_window_jacknife:
    input:
        "data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells_{celltype}/age_{window}_vst.rds"
    output:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/varpart_objects/window_{window}_{celltype}_{boot}_varpart.rds",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/metadata_subsets/window_{window}_{celltype}_{boot}_meta.tsv",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/subject_boots/window_{window}_{celltype}_{boot}_subjects.txt"
    params:
        boot = lambda wildcards: wildcards.boot
    script:
        "scripts/run_varpart_snakemake.R"

rule run_varpart_age_group_jacknife_regress_cells:
    input:
        "data/analysis_out/variancePartition/age_group_regress_cells/subsetted_vst_regress_cells_{celltype}/age_{window}_vst.rds"
    output:
        "data/analysis_out/variancePartition/age_group_bootstrap_regress_cells/varpart_objects/window_{window}_{celltype}_{boot}_varpart.rds",
        "data/analysis_out/variancePartition/age_group_bootstrap_regress_cells/metadata_subsets/window_{window}_{celltype}_{boot}_meta.tsv",
        "data/analysis_out/variancePartition/age_group_bootstrap_regress_cells/subject_boots/window_{window}_{celltype}_{boot}_subjects.txt"
    params:
        boot = lambda wildcards: wildcards.boot
    script:
        "scripts/run_varpart_snakemake_regress_cells.R"

rule metafor_metaanalysis_v3_regress_cells:
    input:
        expand("data/analysis_out/variancePartition/age_sliding_window_bootstrap/varpart_objects/window_{window}_{celltype}_{boot}_varpart.rds", window=WINDOWS, boot=BOOT_ITER),
    output:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/bootstrap_summary_dat_list_{celltype}.rds",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/results_dat_list_with_intercept_{celltype}.rds"
    script:
        "scripts/metafor_metaanalysis_v3_regress_cells.R"

rule sliding_window_fgsea_regress_cells:
    input:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/results_dat_list_with_intercept_{celltype}.rds",
        COMBINED_GENESETS_IN_PATH="data/combined_gene_sets.RDS"
    output:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/fgsea/slope_{celltype}.tsv"
    script:
        "scripts/run_fgsea_v3_regress_cells.R"
