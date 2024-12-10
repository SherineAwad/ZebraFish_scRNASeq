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

cell = args[2] 

DefaultAssay(myObject) <- "RNA"
cell_values <- c(cell)
mySubset <- subset(myObject, idents = cell_values, invert = FALSE)
head(mySubset)

mySubset <- FindNeighbors(mySubset, dims = 1:15)
mySubset <- FindClusters(mySubset, resolution = 3.0)
mySubset <- RunUMAP(object = mySubset, dims = 1:15)

table(mySubset@active.ident)
table(mySubset@meta.data[,'sample'])

sep = paste("_", cell, sep="") 
figure_name <- ""
figure_name <- paste(mysample, "_UMAP.pdf", sep=sep)
pdf(file =figure_name, width =12)
DimPlot(mySubset, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(mySubset, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
DimPlot(mySubset, reduction = "umap", split.by="sample", repel = TRUE) +ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, ".rds", sep=sep)
saveRDS(mySubset, file = myRDS)





