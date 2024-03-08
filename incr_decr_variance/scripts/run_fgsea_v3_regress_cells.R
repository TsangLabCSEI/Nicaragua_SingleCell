library(tidyverse)
library(readr)
library(fgsea)
library(data.table)
library(cowplot)


#Set paths ---------------------------------------------------------------
#Inputs
RES.IN.PATH <- "data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/results_dat_list_with_intercept.rds"
COMBINED.GENESETS.IN.PATH <- "ref/genesets/processed/combined_gene_sets.RDS"

#Out path
TSV.OUT.DIR <- "data/analysis_out/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/fgsea/"
dir.create(TSV.OUT.DIR, recursive = TRUE)
FIG.OUT.DIR <- "plots/variancePartition/age_sliding_window_bootstrap_regress_cells/gls_model_v3_include_nsubj/fgsea/"
dir.create(FIG.OUT.DIR, recursive = TRUE)

#Load data ---------------------------------------------------------------
res <- readRDS(RES.IN.PATH)

slope_dat <- res$subj_div_subj_age_gender_resid %>% filter(term == "window_number_demeaned")
int_dat <- res$subj_div_subj_age_gender_resid %>% filter(term == "intrcpt")

genesetLL <- readRDS(COMBINED.GENESETS.IN.PATH)

##Change rownames to gene names -------------------------------------------
#varpart <- as.data.frame(varpart)

#concatenate geneset list into single list
names(genesetLL$reactome) <- paste0("reactome_", names(genesetLL$reactome))
names(genesetLL$btms) <- paste0("btm_", names(genesetLL$btms))

geneset.list <- Reduce(c, genesetLL)

#Run enrichment -----------------------------------------------------------
universe <- rownames(slope_dat$gene)

dat_list <- list(slope= slope_dat, intercept = int_dat)

set.seed(1)
for(column in names(dat_list)){
        tstat <- dat_list[[column]]$zval
        names(tstat) <- dat_list[[column]]$gene

        fgseaRes <- fgsea(pathways = geneset.list, 
                          stats = tstat,
                          minSize=15,
                          maxSize=500,
                          nperm=100000)

        topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=10), pathway]
        topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=10), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

        plot.path <- paste0(FIG.OUT.DIR, column, ".pdf")
        pdf(plot.path, height =5, width = 12)
        plotGseaTable(geneset.list[topPathways], tstat, fgseaRes, 
                      gseaParam = 0.5)
        dev.off()

        #Save ----------------------------------------------------------------------
        tsv.out.path <- paste0(TSV.OUT.DIR, column, ".tsv" )
        #fwrite(fgseaRes[padj < .05], file=tsv.out.path, sep="\t", sep2=c("", " ", ""), 
        fwrite(fgseaRes, file=tsv.out.path, sep="\t", sep2=c("", " ", ""), 
               nThread = 1)

}


library(tidyverse)
fgseaRes %>% 
        arrange(NES) %>%
        select(-leadingEdge) %>%
        head(20)

fgseaRes %>% 
        arrange(-NES) %>%
        select(-leadingEdge) %>%
        head(20)
