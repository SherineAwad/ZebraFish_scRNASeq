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

#head(myObject) 
#levels(myObject)
#head(Idents(myObject) 


myObject <- RenameIdents(
  object = myObject,
  "0" = "Rod", 
  "3" = "Rod", 
  "9" =  "Rod",
  "25" = "Rod",
  "31" = "Rod",
  "34" = "Rod", 
  "41" = "Rod", 

  "13" = "BC", 
  "17"  = "BC",
  "30"  = "BC",
  "36"  = "BC",
  "48"  = "BC",

  "44" = "RGC", 
  "26" = "RGC", 
  "38" = "RGC", 
 
  "23" = "AC", 
  "28" = "AC",
  "29" = "AC",
  "35" = "AC",
  "37" = "AC",
  "40" = "AC",
  "43" = "AC",
  "45" = "AC",
  "46" = "AC", 
  "49" = "AC", 
  "50" = "AC",

  "12" = "Progenitors", 
  "18" = "Progenitors",
  "19" = "Progenitors",
  "20" = "Progenitors",
  "21" = "Progenitors",
  "24" = "Progenitors", 
  "27" = "Progenitors",
  "1" = "MG", 
  "2" = "MG", 
  "4"  = "MG",
  "5"  = "MG",
  "6" = "MG",
  "8" = "MG",
  "10" = "MG",
  "11" = "MG",
  "15"  = "MG",
  "22"  = "MG",
  "47" = "MG",
  
  "7" = "HC", 
  "16" = "HC",
  "33" = "HC",
  "42" = "HC",


  "14" = "Cones",
  "32"  = "Cones",
  "39"  = "Cones") 

figure_name <- ""
figure_name <- paste(mysample, "_SAnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_Sannotated.rds", sep="")
saveRDS(myObject, file = myRDS)

