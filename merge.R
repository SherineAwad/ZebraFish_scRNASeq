library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

set.seed(1234)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the output RDS filename
myRDS <- paste0(args[1], ".rds")

# The rest are sample names
samples <- args[-1]  # Exclude the first argument (output filename)

# Load RDS files dynamically
sample_objects <- list()

for (sample in samples) {
  sample_file <- paste0(sample, "_filtered.rds")
  sample_objects[[sample]] <- readRDS(sample_file)
  
  # Update metadata
  sample_objects[[sample]]@meta.data[,'sample'] <- apply(
    as.matrix(rownames(sample_objects[[sample]]@meta.data)), 1, function(X1) {
      X2 <- strsplit(X1, '-')[[1]][2]
    }
  )
}

# Merge objects dynamically
myObject <- Reduce(function(x, y) merge(x, y), sample_objects)

# Save the merged object
saveRDS(myObject, file = myRDS)


