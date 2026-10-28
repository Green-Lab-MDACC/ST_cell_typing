rm(list = ls())

library(tidyverse)
library(Seurat)
library(randomForest)

# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^3)


rfPred <- function(x_embed, y, query_cells = NULL) {
  # labels
  stopifnot(all.equal(names(y), rownames(x_embed)))
  y[y == "NA"] <- NA_character_
  
  classes <- y[!is.na(y)] %>% sort() %>% unique()
  
  # split
  train_cells <- names(y[!is.na(y)])
  pred_cells <- names(y[is.na(y)])
  
  train_x <- x_embed[train_cells, ]
  train_y <- y[train_cells]
  train_y <- factor(train_y, levels = classes)
  
  # train randomForest
  set.seed(42)
  rf.model <- randomForest(x = train_x, y = train_y)
  print(rf.model)
  
  # pred
  if (is.null(query_cells)) {
    pred_x <- x_embed[pred_cells, ]
  } else {
    # final pred of all query cells
    stopifnot(length(setdiff(query_cells, rownames(x_embed))) == 0)
    pred_x <- x_embed[query_cells, ]
  }
  
  pred_labels  <- predict(rf.model, pred_x)
  pred_probs <- predict(rf.model, pred_x, type = "prob")
  pred_class_prob <- pred_probs[cbind(1:nrow(pred_probs), match(pred_labels, colnames(pred_probs)))]
  
  tb <- enframe(pred_labels,
                name = "cell",
                value = "celltype") %>%
    mutate(celltype = as.character(celltype)) %>%
    mutate(celltype_prob = pred_class_prob)
  
  # ref + query
  tb <- enframe(train_y,
                name = "cell",
                value = "celltype") %>%
    mutate(celltype = as.character(celltype)) %>%
    mutate(celltype_prob = 1) %>%
    filter(!(cell %in% tb$cell)) %>%
    bind_rows(tb)
  
  tb <- tibble(cell = rownames(x_embed)) %>%
    left_join(tb, by = "cell")
  
  return(tb)
}


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
sobj <- readRDS("../../LymphoMAP_snRNA_ref/Seurat_objects_randcells/scRNALymphoMAP_ref_all_lineage.rds")
levels(sobj$celltype) <- ifelse(levels(sobj$celltype) == "BEC/HEV",
                                "BEC",
                                levels(sobj$celltype))
sobj_ref <- sobj

# load data
sobj <- readRDS("../step1_st_prepLabel/smi_data_prep.rds")

# remove Epi
tb <- read_tsv("../step1_st_prepLabel/cellRefine_knnLabel_rfPred.tsv")
tb <- tb %>%
  filter(celltype != "Epi")
sobj <- sobj[, colnames(sobj) %in% tb$cell]

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
tb <- read_tsv("../../LymphoMAP_snRNA_ref/deg/deg_celltype.tsv")

stopifnot(all.equal(tb$celltype %>% unique() %>% sort(),
                    sobj_ref$celltype %>% as.character() %>% unique() %>% sort()))

features <- tb %>%
  filter(pvals_adj < 0.05 & logfoldchanges > 0.58) %>%
  group_by(celltype) %>%
  arrange(pvals_adj, desc(logfoldchanges)) %>%
  do(head(., n = 100)) %>%
  .$gene %>%
  unique() %>%
  sort()

stopifnot(length(setdiff(features, rownames(sobj))) == 0)


sobj_list <- PrepSCTIntegration(object.list = sobj_list, anchor.features = features)
sobj_list <- lapply(X = sobj_list, FUN = RunPCA, features = features)

n_dims <- 50
k.anchor <- 20 # 50 tested: not improved
sobj.anchors <- FindIntegrationAnchors(object.list = sobj_list, normalization.method = "SCT",
                                       anchor.features = features, dims = 1:n_dims, reduction = "rpca", k.anchor = k.anchor)
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


# augment high-conf cellRefineKnnLabel
tb <- read_tsv("../step1_st_prepLabel/cellRefine_knnLabel.tsv")
tb <- tb %>%
  mutate(cellRefineKnnLabel = ifelse(celltype_fract > 0.75, celltype, NA_character_)) %>%
  select(cell, cellRefineKnnLabel)
message("fract of high-conf labels")
sum(!is.na(tb$cellRefineKnnLabel))/nrow(tb)

tb <- tibble(cell = colnames(sobj)) %>%
  left_join(tb, by = "cell")

sobj$cellRefineKnnLabel <- tb$cellRefineKnnLabel

pdf("scRNALymphoMAP_coembed_cellRefineKnnLabel_umap.pdf", height = 5, width = 7)
DimPlot(sobj[, sobj$dataset == "query"], group.by = "cellRefineKnnLabel", cols = celltype_cols)
dev.off()

# labels
sobj$celltype <- ifelse(!is.na(sobj$celltype),
                        sobj$celltype,
                        sobj$cellRefineKnnLabel)

# merge Mac/Mo before rfPred
sobj$celltype <- ifelse(!is.na(sobj$celltype) & sobj$celltype == "Mo",
                        "Mac",
                        sobj$celltype)

tb <- tibble(cell = colnames(sobj),
             celltype = sobj$celltype)
tb_knnLabel <- tb


# rfPred -----
x_embed <- sobj@reductions$pca@cell.embeddings

# labels
tb <- tb_knnLabel
stopifnot(all.equal(tb$cell, colnames(sobj)))

y <- tb %>% deframe()
y[y == "NA"] <- NA_character_

# pred
tb <- rfPred(x_embed, 
             y, 
             query_cells = NULL)

tb <- tibble(cell = rownames(sobj@meta.data)) %>%
  left_join(tb, by = "cell")

# write query
tb %>%
  filter(cell %in% colnames(sobj)[sobj$dataset == "query"]) %>%
  write_tsv("scRNALymphoMAP_coembed_label_query.tsv")

stopifnot(sum(is.na(tb$celltype)) == 0)

# plot
sobj$celltype <- tb$celltype

# reformat celltype
celltypes <- sobj$celltype %>% as.character()
celltypes <- ifelse(is.na(celltypes), "NA", celltypes) # NA is not appropriately visualized in ImageDimPlot
celltype_cats <- celltypes %>% unique() %>% sort()
celltype_cats <- c(celltype_cats[celltype_cats != "NA"], "NA")
celltypes <- factor(celltypes, levels = celltype_cats)
sobj$celltype <- celltypes

# # use Seurat default color (ggplot2 default color)
# celltype_cols <- structure(c(scales::hue_pal()(length(celltype_cats) - 1), "grey50"),
#                            names = c(celltype_cats))

pdf("scRNALymphoMAP_coembed_label_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "celltype", cols = celltype_cols)
DimPlot(sobj[, sobj$dataset == "ref"], group.by = "celltype", cols = celltype_cols) +
  ggtitle("ref")
DimPlot(sobj[, sobj$dataset == "query"], group.by = "celltype", cols = celltype_cols) +
  ggtitle("query")
dev.off()


# query
sobj <- sobj[, sobj$dataset == "query"]

# deg
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("scRNALymphoMAP_coembed_label_query_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("scRNALymphoMAP_coembed_label_query_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


sessionInfo()
