library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BSgenome.Drerio.UCSC.danRer11)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
library(AcidGenomes)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
myRDS <- args[1]

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]
print(mysample)

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"


all.genes <- rownames(myObject)
myObject<- NormalizeData(myObject)
myObject <- FindVariableFeatures(myObject, selection.method = "vst", nfeatures = 2000)
myObject<- ScaleData(myObject, features = all.genes)
myObject <- RunPCA(myObject, features = VariableFeatures(object = myObject))
myObject <- FindNeighbors(myObject, dims = 1:15)
myObject <- FindClusters(myObject, resolution = 3)
myObject <- RunUMAP(myObject, metric = "correlation",  n.components = 2L, dims = 1:15)

myObject@meta.data[,'sample']<-apply(as.matrix(rownames(myObject@meta.data)), 1, function(X1){X2<-strsplit(X1,'-')[[1]][2];})
myObject@meta.data[,'sample'] <- recode(myObject@meta.data[,'sample'], "1" = "S1", "2" = "S2")


figure_name <- ""
figure_name <- paste(mysample, "_UMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)

