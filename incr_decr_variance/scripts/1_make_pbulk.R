
pbulk <- readRDS("../../../NICA_batch1to4_pseudobulk_normalized.rds")

outdir <- "../data/pbulk_objects/"
dir.create(outdir)

for(nm in names(pbulk)){
  f <- paste0(outdir, nm, ".rds")
  saveRDS(pbulk[[nm]], f)
}


