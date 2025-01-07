library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
args <- commandArgs(trailingOnly = TRUE)

myRDS <- args[1] 
mysample1 <- args[2] 
mysample2 <- args[3] 
mysample3 <- args[4] 
mysample4 <- args[5] 
mysample5 <- args[6] 
mysample6 <- args[7] 
mysample7 <- args[8] 
mysample8 <- args[9] 

myObject1 <- paste(mysample1,"_filtered.rds", sep="")
myObject2 <- paste(mysample2,"_filtered.rds", sep="") 
myObject3 <- paste(mysample3,"_filtered.rds", sep="")
myObject4 <- paste(mysample4,"_filtered.rds", sep="")
myObject5 <- paste(mysample5,"_filtered.rds", sep="")
myObject6 <- paste(mysample6,"_filtered.rds", sep="")
myObject7 <- paste(mysample7,"_filtered.rds", sep="")
myObject8 <- paste(mysample8,"_filtered.rds", sep="")

myRDS <- paste(myRDS, ".rds", sep="")

myObject1
myObject2
myObject3
myObject4
myObject5
myObject6
myObject7
myObject8

myRDS 

myObject1 <- readRDS(myObject1) 
myObject2 <- readRDS(myObject2) 
myObject3 <- readRDS(myObject3) 
myObject4 <- readRDS(myObject4) 
myObject5 <- readRDS(myObject5) 
myObject6 <- readRDS(myObject6) 
myObject7 <- readRDS(myObject7)
myObject8 <- readRDS(myObject8)

myObject1@meta.data[,'sample']<-apply(as.matrix(rownames(myObject1@meta.data)), 1, function(X1){X2<-strsplit(X1,'-')[[1]][2];})
myObject1@meta.data[,'sample'] <- recode(myObject1@meta.data[,'sample'], "1" = "S1", "2" = "S2")


myObject <- merge(myObject1, y = c(myObject2, myObject3, myObject4, myObject5, myObject6, myObject7, myObject8), add.cell.ids = c("mysample1", "mysample2", "mysample3","mysample4", "mysample5", "mysample6", "mysample7", "mysample8"), project = "mergedZebrafish")
saveRDS(myObject, file = myRDS)






