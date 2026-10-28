rm(list = ls())

library(tidyverse)
library(Seurat)
library(randomForest)
library(harmony)

# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^3)


sel_lin <- "NHC"


# celltype colors used in LymphoMAP
celltype_cols <- tribble(
  ~celltype, ~color,
  "B", "#fdbf6f",
  "CD4T", "#33a02c",
  "CD8T", "#e31a1c",
  "NK", "#a6cee3",
  "Mac", "#1f78b4",
  "Mo", "#ff7f00",
  "DC", "#b2df8a",
  "pDC", "#fb9a99",
  "NESC", "#b15928",
  "BEC", "#6a3d9a",
  "LEC", "#ffee00",
  "Epi", "#0000FF",
  "NA", "#7f7f7f" # "grey50"
) %>%
  deframe()
sorted_celltypes <- names(celltype_cols)


# cells from scRNALymphoMAP
sobj <- readRDS(paste0("../../LymphoMAP_snRNA_ref/Seurat_objects_randcells/scRNALymphoMAP_ref_lineage_", 
                       sel_lin, 
                       ".rds"))
levels(sobj$celltype) <- ifelse(levels(sobj$celltype) == "BEC/HEV",
                                "BEC",
                                levels(sobj$celltype))
sel_celltypes <- sobj$celltype %>% unique() %>% sort()
sobj_ref <- sobj

# load data
sobj <- readRDS("smi_data_prep.rds")
sobj_query <- sobj

# common genes
genes <- intersect(rownames(sobj_ref),
                   rownames(sobj_query))

# Create object
sobj <- CreateSeuratObject(
  counts = sobj_ref@assays$originalexp$counts[genes, ],
  meta.data = sobj_ref@meta.data,
  project = "ref"#,
  # min.cells = 3,    # drop genes seen in <3 cells
  # min.features = 200
)
sobj$dataset <- "ref"
sobj_ref <- sobj

sobj <- CreateSeuratObject(
  counts = sobj_query@assays$Nanostring$counts[genes, ],
  meta.data = sobj_query@meta.data,
  project = "query"#,
  # min.cells = 3,    # drop genes seen in <3 cells
  # min.features = 200
)
sobj$dataset <- "query"
sobj_query <- sobj

sobj_list <- as.list(c(sobj_ref, sobj_query))


# Perform integration
sobj_list <- lapply(X = sobj_list, FUN = SCTransform, assay = "RNA", clip.range = c(-10, 10), variable.features.n = 3000, verbose = FALSE)
# features <- SelectIntegrationFeatures(object.list = sobj_list, nfeatures = 3000)


# features
tb <- read_tsv("../../LymphoMAP_snRNA_ref/deg/deg_celltypeInLineage.tsv")
tb <- tb %>%
  mutate(celltype = str_replace(celltype, "BEC/HEV", "BEC"))

stopifnot(length(setdiff(sel_celltypes, tb$celltype)) == 0)
tb <- tb %>%
  filter(celltype %in% sel_celltypes)
features <- tb %>%
  group_by(celltype) %>%
  arrange(pvals_adj, desc(logfoldchanges)) %>%
  do(head(., n = 100)) %>%
  .$gene %>%
  unique()

stopifnot(length(setdiff(features, rownames(sobj))) == 0)


sobj_list <- PrepSCTIntegration(object.list = sobj_list, anchor.features = features)
sobj_list <- lapply(X = sobj_list, FUN = RunPCA, features = features)

n_dims <- 50
k.anchor <- 20 # 50 tested: not improved
sobj.anchors <- FindIntegrationAnchors(object.list = sobj_list, normalization.method = "SCT",
                                       reference = 1, # index of reference object in object.list so that the reference SCT model is kept
                                       anchor.features = features, dims = 1:n_dims, reduction = "cca", k.anchor = k.anchor)
sobj <- IntegrateData(anchorset = sobj.anchors, normalization.method = "SCT", dims = 1:n_dims)
DefaultAssay(sobj) <- "integrated"
sobj <- RunPCA(sobj, npcs = n_dims, verbose = FALSE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:n_dims)

pdf("scRNALymphoMAP_coembed_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "celltype", cols = celltype_cols)
DimPlot(sobj[, sobj$dataset == "ref"], group.by = "celltype", cols = celltype_cols) +
  ggtitle("ref")
DimPlot(sobj, group.by = "dataset")
dev.off()

# save
sobj %>%
  saveRDS("scRNALymphoMAP_coembed.rds")


sessionInfo()
