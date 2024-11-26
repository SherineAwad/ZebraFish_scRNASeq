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

myRDS <- paste(mysample, "_analysed.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

#head(myObject) 
#levels(myObject)
#head(Idents(myObject) 


myObject <- RenameIdents(
  object = myObject,
  "1" = "MG", 
  "2" = "MG", 
  "3" = "MG", 
  "4" = "MG",
  "5" = "MG", 
  "7" = "MG", 
  "8" = "MG", 
  "9" = "MG", 
  "15" = "MG", 
  "16" = "MG", 
  "21" = "MG", 
  "43" = "MG", 
  "48" = "MG", 
  "53" = "MG", 
  
  "13" = "Rod", 
  "14" = "Rod", 
  "18" = "Rod",
  "39" = "Rod", 

  "20" = "Cone",
  "41" = "Cone", 

  "12" = "BC",
  "28" = "BC",
  "40" = "BC", 
  "44" = "BC", 
  "45" = "BC", 
  "55" = "BC",  
  
  "11" = "AC", 
  "22" = "AC",
  "34" = "AC", 
  "35" = "AC", 
  "49"= "AC", 
   
  "17" = "RGC", 

  "24" = "HC",
  "26" = "HC",
  "29" = "HC",
  "30" = "HC",
  "36" = "HC",
  "42" = "HC",
  "46" = "HC",
  "47" = "HC",

  "0" = "Progenitors", 
  "6" = "Progenitors", 
  "10" = "Progenitors", 
  "19" = "Progenitors", 
  "23" = "Progenitors", 
  "25" = "Progenitors",
  "27" = "Progenitors",
  "31" = "Progenitors",
  "32" = "Progenitors",
  "37" = "Progenitors",
  "52" = "Progentiors")

figure_name <- ""
figure_name <- paste(mysample, "_AnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = myRDS)

