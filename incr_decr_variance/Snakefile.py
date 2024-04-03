singularity: "file:///hpcdata/sg/sg_data/PROJECTS/Nicaragua_Project/analysis/singularity/intrinsicness_0.1.sif" 

CELLTYPES = ["B_Mem","B_Naive","CD4_Mem","CD4_Naive","CD8_Mem","CD8_Naive","cDC","gdT_Vd1","gdT_Vd2","MAIT","Mono_Classical","Mono_NonClassical","NK_CD16hi","NK_CD56hi","pDC"]

WINDOWS = []

for i in range(3,12):
    WINDOWS.append(str(i) + "to" + str(i + 3))

BOOT_ITER = range(0,20)

rule all:    
    input:
        expand("data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/fgsea/{celltype}/slope.tsv", celltype=CELLTYPES),
        #expand("data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/bootstrap_summary_dat_list_{celltype}.rds", celltype=CELLTYPES)
        #expand("data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells/{celltype}/age_{window}_vst.rds", window = WINDOWS, celltype = CELLTYPES)

## processsing -------------------------------------------------------

#sliding window analysis --------------------------------------

rule make_window_subsets_rm_cells:
    input:
        "data/pbulk_objects/{celltype}.rds"
    output:
        expand("data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells/{{celltype}}/age_{window}_vst.rds", window = WINDOWS)
    script:
        "scripts/subset_data_ctrl_cells.R"

rule run_varpart_window_jacknife:
    input:
        expand("data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells/{{celltype}}/age_{window}_vst.rds", window = WINDOWS)
    output:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/varpart_objects/{celltype}/window_{window}_{boot}_varpart.rds",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/metadata_subsets/{celltype}/window_{window}_{boot}_meta.tsv",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/subject_boots/{celltype}/window_{window}_{boot}_subjects.txt"
    #params:
    #    celltype = lambda wildcards: wildcards.celltype,
    #    window = lambda wildcards: wildcards.window,
    #    boot = lambda wildcards: wildcards.boot
    script:
        "scripts/run_varpart_snakemake.R"

rule metafor_metaanalysis_v3_regress_cells:
    input:
        expand("data/analysis_out/variancePartition/sliding_age_window/subsetted_vst_ctrl_cells/{{celltype}}/age_{window}_vst.rds", window = WINDOWS),
        expand("data/analysis_out/variancePartition/age_sliding_window_bootstrap/varpart_objects/{{celltype}}/window_{window}_{boot}_varpart.rds", window=WINDOWS, boot=BOOT_ITER),
    output:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/bootstrap_summary_dat_list__{celltype}.rds",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/results_dat_list_with_intercept__{celltype}.rds",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/combined_jaccmat__{celltype}.pdf",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/combined_jaccmat__{celltype}.tsv"
    params:
        celltype = lambda wildcards: wildcards.celltype,
    script:
        "scripts/metafor_metaanalysis_v3_regress_cells.R"

rule sliding_window_fgsea_regress_cells:
    input:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/results_dat_list_with_intercept__{celltype}.rds",
        COMBINED_GENESETS_IN_PATH="data/combined_gene_sets.RDS"
    output:
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/fgsea/{celltype}/slope.tsv",
        "data/analysis_out/variancePartition/age_sliding_window_bootstrap/gls_model_v3_include_nsubj/fgsea/{celltype}/slope.pdf"
    params:
        celltype = lambda wildcards: wildcards.celltype,
    script:
        "scripts/run_fgsea_v3_regress_cells.R"
