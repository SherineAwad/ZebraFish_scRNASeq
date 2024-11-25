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
mysample <- args[1]

#Run on merged 
myRDS <- paste(mysample, "_filtered.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"

myObject <- NormalizeData(myObject)
myObject <- FindVariableFeatures(myObject)
myObject <- ScaleData(myObject)
myObject <- RunPCA(myObject)
myObject <- FindNeighbors(myObject, dims = 1:8)
myObject <- FindClusters(myObject, resolution = 2.0)
myObject <- RunUMAP( myObject, dims = 1:8)

myObject[["RNA"]] <- JoinLayers(myObject[["RNA"]])

figure_name <- ""
figure_name <- paste(mysample, "_UMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)

