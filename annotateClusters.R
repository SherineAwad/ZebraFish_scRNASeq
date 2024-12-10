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


myObject <- RenameIdents(
  object = myObject,
  "0" = "MG", 
  "2" = "MG", 
  "4" = "MG", 
  "5" = "MG", 
  "6" = "MG", 
  "8" = "MG",
  "9" = "MG", 
  "11" = "MG", 
  "15" = "MG", 
  "22" = "MG", 
  "25" = "MG", 
  "28" = "MG", 
  "32" = "MG", 
  "35" = "MG", 
  "46" = "MG", 
  "57" = "MG", 

  "1" = "Rod", 
  "3" = "Rod", 
  "7" = "Rod", 
  "30" = "Rod",  
  "38" = "Rod",  
  "49" = "Rod",  
  "53" = "Rod",  
  "10" = "Rod", 

  "21" = "Cone", 
  "37" = "Cone",
  "55" = "Cone",


  "12" = "HC", 
  "29" = "HC",
  "39" = "HC", 
  "40" = "HC", 
  "43" = "HC", 
  "54" = "HC", 
  "59" = "HC", 

  "47" = "BC", 
  "14"  = "BC",
  "13"  = "BC",
  "14"  = "BC", 
  "36"  = "BC", 
  "44"  = "BC", 
  "61"  = "BC", 



  "34" = "RGC", 
  "23" = "RGC", 

  "26" = "AC", 
  "27" = "AC",
  "31" = "AC", 
  "45" = "AC",
  "50" = "AC",
  "51" = "AC",
  "52" = "AC", 
  "56" = "AC", 
  "60" = "AC",
  "63" = "AC", 
  "58" = "AC", 
  "33" = "AC", 
  "62" = "AC", 


  "64" = "Microglia",

  "16" = "Progenitors", 
  "17" = "Progenitors",
  "18" = "Progenitors",
  "19" = "Progenitors", 
  "20" = "Progenitors", 
  "24" = "Progenitors", 
  "41" = "Progenitors", 
  "42" = "Progenitors", 
  "48" = "Progenitors") 


figure_name <- ""
figure_name <- paste(mysample, "_AnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = myRDS)

