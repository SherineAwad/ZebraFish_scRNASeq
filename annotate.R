library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

myRDS <- args[1]
annotation_file <- args[2] 

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]
print(mysample)

# Read Seurat object
myObject <- readRDS(myRDS)

# Read annotations from file
annotations <- read.csv(annotation_file, header = FALSE, col.names = c("cluster", "label"))

# Convert to named vector
annotation_vector <- setNames(as.character(annotations$label), as.character(annotations$cluster))

# Apply dynamic renaming
myObject <- RenameIdents(myObject, annotation_vector)


figure_name <- ""
figure_name <- paste(mysample, "_annotatedUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

# Save updated object
saveRDS(myObject, file = paste0(mysample, "_annotated.rds"))


