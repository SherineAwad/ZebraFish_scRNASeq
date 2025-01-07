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
  "6" = "Rod", 
  "11" =  "Rod",
  "24" = "Rod",
  "48" = "Rod", 
  "29" = "Rod", 
  "49" = "Rod",

  "13" = "BC", 
  "15"  = "BC",
  "31"  = "BC",
  "38"  = "BC",
  "39"  = "BC",

  "23" = "RGC", 
  "25" = "RGC", 
 
  "18" = "AC", 
  "26" = "AC",
  "34" = "AC",
  "35" = "AC",
  "37" = "AC",
  "40" = "AC",
  "41" = "AC",
  "42" = "AC",
  "43" = "AC", 
  "45" = "AC",
  "47" = "AC",

  "9" = "Progenitors", 
  "12" = "Progenitors",
  "14" = "Progenitors",
  "17" = "Progenitors",
  "21" = "Progenitors",
  "22" = "Progenitors", 
  "30" = "Progenitors",
  
  "7" = "HC", 
  "19" = "HC",
  "36" = "HC",

  "1" = "MG", 
  "2"  = "MG",
  "3"  = "MG",
  "4"  = "MG",
  "5"  = "MG",
  "8"  = "MG",
  "10"  = "MG",
  "16"  = "MG",
  "27"  = "MG",
  "28"  = "MG",
  "32"  = "MG",

  "20" = "Cones",
  "33"  = "Cones",
  "44"  = "Cones",
  "46" = "Cones") 

figure_name <- ""
figure_name <- paste(mysample, "_SAnnotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_Sannotated.rds", sep="")
saveRDS(myObject, file = myRDS)

