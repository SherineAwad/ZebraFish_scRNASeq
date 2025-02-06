library(Seurat)
library(dplyr)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

myRDS <- args[1]   # First argument: input Seurat object file
annotation_file <- args[2]  # Second argument: file with annotations

# Read Seurat object
myObject <- readRDS(myRDS)

# Read annotations from file
annotations <- read.csv(annotation_file, header = FALSE, col.names = c("cluster", "label"))

# Convert to named vector
annotation_vector <- setNames(as.character(annotations$label), as.character(annotations$cluster))

# Apply dynamic renaming
myObject <- RenameIdents(myObject, annotation_vector)

# Save updated object
saveRDS(myObject, file = paste0("annotated_", myRDS))



