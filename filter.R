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

mysample = args[1] 
#filter first round 
nRNA1 =  35000
nRNA2 =  250
features = 250  
mt = 25

nRNA1
nRNA2
features 
mt

myRDS = paste(mysample, "_preprocessed.rds", sep="") 
myRDS
myObject <- readRDS(myRDS)

myObject <- subset(x = myObject, nCount_RNA < nRNA1 & nCount_RNA > nRNA2 & nFeature_RNA > features  & percent.mt < mt)
figure_name <- paste(mysample, "_QC_afterFiltervlnplot.pdf", sep="")
pdf(file=figure_name, width=12)
VlnPlot(object = myObject, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),pt.size=0.1, ncol = 6)
dev.off()
myRDS <- paste(mysample, "_filtered.rds", sep="") 
saveRDS(myObject, file = myRDS)



