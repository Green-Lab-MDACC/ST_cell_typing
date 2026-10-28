rm(list = ls())

library(tidyverse)
library(Seurat)
library(patchwork)
library(ggpubr)

# seurat5


# input arguments
sobj_name <- "../../example_data/sobj_NESC.rds"
assay_name <- "Nanostring"


# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^3)

# load
sobj <- readRDS(sobj_name)

# Norm, embed, cluster -------
sobj <- SCTransform(sobj, assay = assay_name, clip.range = c(-10, 10), variable.features.n = 3000, verbose = FALSE)

npcs <- 50
sobj <- RunPCA(sobj, npcs = npcs, verbose = FALSE)

# variance explained
var.features <- VariableFeatures(sobj)
total_variance <- sum(matrixStats::rowVars(sobj@assays$SCT@scale.data[var.features, ]))
eigValues <- (sobj@reductions$pca@stdev)^2  ## EigenValues
varExplained <- eigValues / total_variance
cumVarExplained <- cumsum(varExplained)

# elbow plot
tb <- tibble(PC = 1:length(cumVarExplained),
             percentVarExplained = cumVarExplained*100)
pdf("smi_data_prep_pca_elbow.pdf", height = 4, width = 5)
ElbowPlot(sobj, ndims = npcs)
ggline(tb, "PC", "percentVarExplained", ylab="% Variance Explained")
dev.off()

# cluster
n_dims <- 50
k.param <- 50

sobj <- RunUMAP(sobj, dims = 1:n_dims, verbose = FALSE)
sobj <- FindNeighbors(sobj, k.param = k.param, dims = 1:n_dims, verbose = FALSE)

sobj <- FindClusters(sobj, resolution = 0.1, verbose = FALSE)

pdf("smi_data_prep_umap_seurat_clusters.pdf", height = 5, width = 6)
DimPlot(sobj, group.by = "seurat_clusters")
dev.off()

# save
sobj %>%
  saveRDS("smi_data_prep.rds")


sessionInfo()
