# set seed
set.seed(4)

# load libraries
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(stringr)

# load single cell seurat object
dataset.query <- readRDS(snakemake@input[["data"]]) 
prots<-rownames(dataset.query)[!(str_detect(rownames(dataset.query),"iso")|str_detect(rownames(dataset.query),"Iso"))]


# vizualize UMAP
Idents(dataset.query)<-dataset.query$manual.cid
dataset.query <- FindVariableFeatures(dataset.query, selection.method = "vst", nfeatures = 2000)
dataset.query <- ScaleData(dataset.query)#, vars.to.regress = "batch")
dataset.query <- RunPCA(dataset.query, features = VariableFeatures(object = dataset.query))
#dataset.query.harmonized <- dataset.query
dataset.query.harmonized <- RunHarmony(dataset.query, "batch")
dataset.query.harmonized <- FindNeighbors(dataset.query.harmonized, dims = 1:15, reduction="harmony", features=VariableFeatures(dataset.query.harmonized))
dataset.query.harmonized <- FindClusters(dataset.query.harmonized, resolution = 1, features=VariableFeatures(dataset.query.harmonized))
dataset.query.harmonized <- RunUMAP(dataset.query.harmonized, dims=1:15, reduction = "harmony") #batch3, n.neighbors = 5L, min.dist = 0.3 #batch4 dims=1:15 #all dims=1:15
pdf("data/output/UMAP.pdf")
DimPlot(dataset.query.harmonized, reduction = "umap", group.by = "manual.cid", label = TRUE, repel = TRUE, raster=FALSE, label.size = 4, shuffle = TRUE)
dev.off()


# aggregate celltypes, individuals, gender and timepoints
nica_list <- list()
for (i in unique(dataset.query$matched.individual)){
  print(i)
  dataset.i<-dataset.query[,which(dataset.query$matched.individual==i)]
  for (t in unique(dataset.i$matched.timepoint.age[which(!is.na(dataset.i$matched.timepoint.age))])){
    dataset.t<-dataset.i[,which(dataset.i$matched.timepoint.age==t)]
    b<-unique(dataset.t$batch)
    g<-unique(dataset.t$gender)
    nica_list[[paste0(i,"_",t,"_",b,"_",g)]]<-table(dataset.t$manual.cid)
  }
}

nica.celltypes<-as.data.frame(do.call(cbind, nica_list))
nica.celltypes<-data.frame(t(nica.celltypes))
nica.celltypes.sums<-rowSums(nica.celltypes,na.rm = TRUE)
for (t in seq(1,dim(nica.celltypes)[1])) {
  nica.celltypes[t,]<-nica.celltypes[t,]/nica.celltypes.sums[t]
}
nica.celltypes$sums<-nica.celltypes.sums
nica.celltypes$individual<-as.vector(t(data.frame(sapply(rownames(nica.celltypes),function(x){str_split(x,"_")[[1]]})))[,1])
nica.celltypes$timepoint<-as.vector(t(data.frame(sapply(rownames(nica.celltypes),function(x){str_split(x,"_")[[1]]})))[,2])
nica.celltypes$timepoint<-as.numeric(nica.celltypes$timepoint)#*12
nica.celltypes$batch<-as.vector(t(data.frame(sapply(rownames(nica.celltypes),function(x){str_split(x,"_")[[1]]})))[,3])
nica.celltypes$gender<-as.vector(t(data.frame(sapply(rownames(nica.celltypes),function(x){str_split(x,"_")[[1]]})))[,4])


#visualize distribution of timepoints per individual
age_df<-nica.celltypes[c("individual","timepoint","gender")]
age_df<-age_df %>% group_by(individual) %>% mutate(min=min(timepoint),max=max(timepoint)) %>% ungroup()
age_df<-age_df[order(age_df$timepoint,age_df$max-age_df$min, decreasing = F),]
age_df$idx<-as.numeric(factor(age_df$individual,levels=unique(age_df$individual)))
age_df %>%
    ggplot(aes(x=timepoint, y = idx, col=gender, group=idx)) +
    geom_point() + 
    geom_line() + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
    scale_x_continuous(breaks = seq(1,20,1)) +
    xlab("Followed age range")
ggsave("data/output/individuals.pdf")

