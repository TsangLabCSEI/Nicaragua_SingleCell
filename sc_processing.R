#Load libraries
library(Seurat)
library(dsb)
library(scCustomize)


#Define protein markers and standards
markers<-read.csv("CITE_protein_markerset.csv")
iso<-c("RatIgG1kiso","RatIgG2akiso","RatIgG2bkIso","ArmenianHamsterIgGiso","IgG2bkiso","IgG2akiso")
markers<-c(markers,iso)


#Iterate over sequencing runs and store demultiplexed data in lists
seurat_object_list <- list()
all_demuxlet_lists <- list()

for (i in c(1,2,3,4,5,6,7,8)){

	sample<-paste0("nica0",i)
	sample10x<-paste0("/cellranger_bash/NICARAGUA_SCMPX0",i,"/outs/per_sample_outs/NICARAGUA_SCMPX0",i,"/count/sample_feature_bc_matrix/")
	print(paste0("... now processing sample ",sample))

	#Read demuxlet output
	dmx<-read.table(paste0("/demuxlet/new.result.",sample,".grch38.best"),header = 1)

	#Filter demuxlet output
	tbmx<-table(dmx[which(dmx$PRB.SNG1>0.8),]$BEST)
	bcds<-dmx[which(dmx$PRB.SNG1>0.8),]$BARCODE
	fdmx<-dmx[which(dmx$PRB.SNG1>0.8),]

	#Store demuxlet output for each patient
	demuxlet_list<-list()

	#control
	demuxlet_list[["CHI14"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_CHI-014"),]$BARCODE
	#batch1
	demuxlet_list[["5763"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_5763_11D"),]$BARCODE
	demuxlet_list[["4416"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_4416_10C"),]$BARCODE
	demuxlet_list[["4916"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_4916_10C"),]$BARCODE
	demuxlet_list[["6936"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_6936_11V"),]$BARCODE
	demuxlet_list[["6959"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_6959_11V"),]$BARCODE
	#batch2
	#demuxlet_list[["4416"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_4416_10C"),]$BARCODE
	#demuxlet_list[["4514"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_4514_10C"),]$BARCODE
	#demuxlet_list[["5689"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_5689_11D"),]$BARCODE
	#demuxlet_list[["5868"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_5868_10D"),]$BARCODE
	#demuxlet_list[["6462"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_6462_10D"),]$BARCODE
	#demuxlet_list[["6853"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_6853_12D"),]$BARCODE
	#batch3
	#demuxlet_list[["4763"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_4763_10C"),]$BARCODE
	#demuxlet_list[["4767"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_4767_10C"),]$BARCODE
	#demuxlet_list[["5862"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_5862_10C"),]$BARCODE
	#demuxlet_list[["6895"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_6895_12V"),]$BARCODE
	#batch4
	#demuxlet_list[["7260"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7260_11c"),]$BARCODE
	#demuxlet_list[["7393"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7393_11c"),]$BARCODE
	#demuxlet_list[["7323"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7323_12d"),]$BARCODE
	#demuxlet_list[["7524"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7524_13c"),]$BARCODE
	#demuxlet_list[["7586"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7586_12c"),]$BARCODE
	#demuxlet_list[["7726"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7726_13c"),]$BARCODE
	#demuxlet_list[["7851"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_7851_12c"),]$BARCODE
	#demuxlet_list[["8429"]]<-fdmx[which(fdmx$BEST=="SNG-NICA_8429_12c"),]$BARCODE

	all_demuxlet_lists[[sample]]<-demuxlet_list

	#Read cellranger mappings
	df<-Read10X(sample10x)
	#batch1
	HTOdf<-df$`Antibody Capture`[c(seq(1,10),14,15),]
	df$`Antibody Capture`<-df$`Antibody Capture`[which(rownames(df$`Antibody Capture`) %in% markers),]
	#batch2
	#HTOdf<-df$`Antibody Capture`[seq(1,10),]
	#batch3
	#HTOdf<-df$`Antibody Capture`[c(seq(1,10),12,13,14),]
	#batch4
	#HTOdf<-df$`Antibody Capture`[seq(1,7),]

	#create and filter seurat object for barcodes
	seurat_object = CreateSeuratObject(counts = df$`Gene Expression`) #[['HTO']]
	seurat_object[['CITE']] = CreateAssayObject(counts = df$`Antibody Capture`)
	seurat_object[['HTO']] = CreateAssayObject(counts = HTOdf)

	seurat_object_list[[sample]]<-seurat_object 
} 

#Merge seurat objects from all lanes
seurat_object_all <- Merge_Seurat_List(list_seurat = seurat_object_list, add.cell.ids = names(seurat_object_list))

#Merge demuxlet barcodes per patient per sample
dmxlt<-as.data.frame(do.call(cbind, all_demuxlet_lists))

#Demultiplex based on HTOs
seurat_object_all <- NormalizeData(seurat_object_all, assay = "HTO", normalization.method = "CLR")
seurat_object_all <- HTODemux(seurat_object_all, assay = "HTO", positive.quantile = 0.99)
Idents(seurat_object_all) <- "HTO_maxID"

#Filter droplets for singles and negatives
seurat_object_allpos<-seurat_object_all[,which(seurat_object_all$HTO_classification.global=="Singlet")]
seurat_object_allneg<-seurat_object_all[,which(seurat_object_all$HTO_classification.global=="Negative")]


#DSB normalization of protein counts
adt_norm = DSBNormalizeProtein(
 cell_protein_matrix = as.matrix(GetAssayData(seurat_object_allpos[["CITE"]], slot = "counts")),
 empty_drop_matrix = as.matrix(GetAssayData(seurat_object_allneg[["CITE"]], slot = "counts")),
 denoise.counts = TRUE,
 use.isotype.control = TRUE,
 isotype.control.name.vec = as.vector(sapply(iso, function(x){paste0("PROT-",x)})))

seurat_object_allpos[["CITE"]] <- SetAssayData(seurat_object_allpos[["CITE"]], slot = "data", new.data = adt_norm)


#Calculcate mitochondrial transcript percentage
seurat_object_allpos[["percent.mt"]] <- PercentageFeatureSet(seurat_object_allpos, pattern = "^MT-")

#Log Normalization of RNA counts
seurat_object_allpos <- NormalizeData(seurat_object_allpos, normalization.method = "LogNormalize", scale.factor = 10000)

#QC Filter based on RNA features, count and percent mt
seurat_object.filtRNA <- subset(seurat_object_allpos, subset = nFeature_RNA > 600 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 15000)

#Prepare dimensionality reductions for seurat object using variable feature detection and PCA 
seurat_object.filtRNA <- FindVariableFeatures(seurat_object.filtRNA, selection.method = "vst", nfeatures = 2000)
seurat_object.filtRNA <- ScaleData(seurat_object.filtRNA, features = all.genes)
seurat_object.filtRNA <- RunPCA(seurat_object.filtRNA, features = VariableFeatures(object = seurat_object.filtRNA))
seurat_object.filtRNA <- FindNeighbors(seurat_object.filtRNA, dims=NULL, reduction=NULL)
seurat_object.filtRNA <- FindClusters(seurat_object.filtRNA, resolution = 1)
seurat_object.filtRNA.filtdmx <- RunUMAP(seurat_object.filtRNA.filtdmx, dims = 1:10)


