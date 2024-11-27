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
  "6" = "MG", 
  "8" = "MG", 
  "10" = "MG", 
  "21" = "MG", 
  "32" = "MG", 
  
  "0" = "Rod", 
  "5" = "Rod", 
  "9" =  "Rod",
  "23" = "Rod",
  "30" = "Rod",
  "31" = "Rod", 
  "37" = "Rod", 
  
  "13" = "Cone", 
  "33" = "Cone",
  "7" = "Cone", 
  "17" = "Cone", 
  "34" = "Cone", 
  "44" = "Cone",  

  "12" = "BC", 
  "18"  = "BC",
  "29"  = "BC",
  "38"  = "BC",
  "40"  = "BC",
  "49"  = "BC",

  "14" = "RGC", 
  "45" = "RGC", 
 
  "24" = "AC", 
  "27" = "AC",
  "28" = "AC",
  "35" = "AC",
  "39" = "AC",
  "41" = "AC",
  "42" = "AC",
  "43" = "AC",
  "47" = "AC",
  "50" = "AC",
  "51" = "AC",

  "48" = "Microglia",
  "11" = "Olignocytes", 

  "15" = "Progenitors", 
  "16" = "Progenitors",
  "19" = "Progenitors",
  "20" = "Progenitors",
  "22" = "Progenitors",
  "25" = "Progenitors", 
  "26" = "Progenitors",
  "36" = "Progenitors") 

figure_name <- ""
figure_name <- paste(mysample, "_AnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = myRDS)

