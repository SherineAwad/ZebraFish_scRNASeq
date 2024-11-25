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
myh5 <- paste(mysample,"_filtered_feature_bc_matrix.h5", sep="")
myRDS <- paste(mysample, ".rds", sep="")

mysample
myh5
myRDS

## Prepare genome annotation
annotation <- makeGRangesFromGff("Danio_rerio.GRCz11.105.gtf", level = c("genes", "transcripts"))
#Need to check why the annotation doesn't work 
#seqlevelsStyle(annotation) <- "UCSC"
#unique(seqnames(annotation)) 

###load sample, the RNA and ATAC data
counts <- Read10X_h5(myh5)
# Check counts
lapply(counts, dim)
lapply(counts, class)

# create a Seurat object containing the RNA adata
myObject <- CreateSeuratObject(counts = counts,assay = "RNA", project=mysample)

myObject[["percent.mt"]] <- PercentageFeatureSet(myObject, pattern = "^mt-", assay="RNA") 
### Check cell quality
figure_name <- ""
figure_name <- paste(mysample, "_QC_beforeFiltervlnplot.pdf", sep="")
pdf(file=figure_name, width=12)
VlnPlot(object = myObject, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),pt.size=0.1, ncol = 6)
dev.off()
myRDS <- paste(mysample, "_preprocessed.rds", sep="")
saveRDS(myObject, file = myRDS)





