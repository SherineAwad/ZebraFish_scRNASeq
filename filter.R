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

nCount1 = args[2] 
nCount2 = args[3]
nfeatures1 = args[4]
nfeatures2 = args[5]
mt = args[6] 

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]
print(mysample)


myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"


myObject <- subset(myObject, subset = nFeature_RNA > nfeatures1 & nFeature_RNA < nfeatures2 & nCount_RNA > nCount1 & nCount_RNA < nCount2 & percent.mt < mt)
figure_name <- paste(mysample, "_QC_afterFiltervlnplot.pdf", sep="")
pdf(file=figure_name, width=12)
VlnPlot(object = myObject, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),pt.size=0.1, ncol = 6, group.by ="orig.ident")
dev.off()
myRDS <- paste(mysample, "_filtered.rds", sep="") 
saveRDS(myObject, file = myRDS)



