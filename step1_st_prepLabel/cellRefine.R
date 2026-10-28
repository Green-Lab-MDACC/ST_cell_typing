rm(list = ls())

library(tidyverse)
library(Seurat)
library(randomForest)

options(future.globals.maxSize = 8000 * 1024^3)


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


# embed
run_embed <- function(sobj, features = NULL, assay = "Nanostring", return.model = FALSE) {
  sobj <- SCTransform(sobj, assay = assay, clip.range = c(-10, 10), variable.features.n = 3000, verbose = FALSE)
  
  npcs <- 50
  sobj <- RunPCA(sobj, npcs = npcs, features = features, verbose = FALSE)
  
  n_dims <- 50
  k.param <- 50
  
  sobj <- RunUMAP(sobj, dims = 1:n_dims, 
                  return.model = return.model,
                  verbose = FALSE)
  sobj <- FindNeighbors(sobj, k.param = k.param, dims = 1:n_dims, verbose = FALSE)
  
  return(sobj)
}


# load data
sobj <- readRDS("smi_data_prep.rds")

# labels
tb <- read_tsv("init_label.tsv")
df <- tibble(cell = colnames(sobj)) %>%
  left_join(tb, by = "cell") %>%
  select(-cell) %>%
  as.data.frame()
sobj@meta.data <- sobj@meta.data[, !(colnames(sobj@meta.data) %in% colnames(df))] %>%
  cbind(df)

sobj$celltype <- sobj$labels

# reformat celltype
celltypes <- sobj$celltype %>% as.character()
celltypes <- ifelse(is.na(celltypes), "NA", celltypes) # NA is not appropriately visualized in ImageDimPlot
celltype_cats <- celltypes %>% unique() %>% sort()
celltype_cats <- c(celltype_cats[celltype_cats != "NA"], "NA")
celltypes <- factor(celltypes, levels = celltype_cats)
sobj$celltype <- celltypes

# embed
if ("SCT_nn" %in% Graphs(sobj)) {
  message("SCT_nn exists, skip run_embed")
} else {
  sobj <- run_embed(sobj, features = NULL, assay = "Nanostring")
}

sobj$labels <- sobj$celltype
pdf("cellRefine_initLabel_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "labels", cols = celltype_cols)
dev.off()


# knnLabel restriction -----
library(igraph)

# knn with fraction vote
fract_vote <- function(cells, labels, fract_cutoff) {
  freq_table <- labels[cells] %>% table()
  fracts <- freq_table/sum(freq_table)
  if (max(fracts) > fract_cutoff) {
    return(names(freq_table)[which.max(freq_table)])
  } else {
    return(NA_character_)
  }
}

fract_score <- function(cells, labels) {
  freq_table <- labels[cells] %>% table()
  fracts <- freq_table/sum(freq_table)
  return(max(fracts))
}

knnLabel <- function(sobj, fract_cutoff = 0.5) {
  # labels
  stopifnot("labels" %in% colnames(sobj@meta.data))
  labels <- sobj$labels
  
  stopifnot(!("NA" %in% labels[!is.na(labels)]))
  labels[is.na(labels)] <- "NA"
  
  # K nearest neighbors graph from Seurat
  mt_nn <- sobj@graphs$SCT_nn #%>% as.matrix()
  
  k_cells <- sum(mt_nn[1, ] > 0)
  stopifnot(all(rowSums(mt_nn > 0) == k_cells))
  
  # check if cell names are consistent
  stopifnot(all.equal(rownames(mt_nn), colnames(mt_nn)))
  stopifnot(all.equal(names(labels), rownames(mt_nn)))
  
  # # get cells of nn
  # pos_names <- function(x) {
  #   names(x)[x > 0]
  # }
  # cells_nn <- apply(mt_nn, 1, pos_names) # this requires high memory
  
  # get cells of nn
  
  g <- graph_from_adjacency_matrix(mt_nn, weighted = TRUE)
  edge_matrix <- as_edgelist(g)
  
  # cells_nn <- edge_matrix %>%
  #   as_tibble() %>%
  #   spread(key = V1, value = V2) %>% # this is slow
  #   as.matrix()
  # cells_nn <- cells_nn[, rownames(mt_nn)]
  
  tb <- edge_matrix %>%
    as_tibble() %>%
    group_by(V1) %>%
    summarise(V2 = list(V2)) %>%
    unnest_wider(V2, names_sep = "_")
  cells_nn <- as.matrix(tb[, -1])
  rownames(cells_nn) <- tb[[1]]
  cells_nn <- t(cells_nn)
  cells_nn <- cells_nn[, rownames(mt_nn)]
  
  knn_labels <- apply(cells_nn, 2, fract_vote, labels = labels, fract_cutoff = fract_cutoff)
  knn_labels[knn_labels == "NA"] <- NA_character_
  
  # knn_labels[is.na(knn_labels)] <- labels[is.na(knn_labels)]
  # use fract-voted labels to reduce false positives from SingleR that doesn't consider unknown cell types
  
  knn_scores <- apply(cells_nn, 2, fract_score, labels = labels)
  stopifnot(all.equal(names(knn_scores), names(knn_labels)))
  
  tb <- tibble(cell = names(knn_labels),
               celltype = knn_labels,
               celltype_fract = knn_scores)
  return(tb)
}

tb <- knnLabel(sobj)
tb %>%
  write_tsv("cellRefine_knnLabel.tsv")

tb_anno <- tb

# umap
df <- tibble(cell = colnames(sobj)) %>%
  left_join(tb, by = "cell") %>%
  select(-cell) %>%
  as.data.frame()
sobj@meta.data <- sobj@meta.data[, !(colnames(sobj@meta.data) %in% colnames(df))] %>%
  cbind(df)

fract_cutoff <- 0.5
sobj$celltype[sobj$celltype_fract < fract_cutoff] <- NA_character_

pdf("cellRefine_knnLabel_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "celltype", cols = celltype_cols) +
  ggtitle(paste0("fract_cutoff ", fract_cutoff, " NA rate ", sum(is.na(sobj$celltype))/ncol(sobj)))
dev.off()


# rfPred -----
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


x_embed <- sobj@reductions$pca@cell.embeddings

# labels
tb <- tb_anno

y <- tb %>% deframe()
y[y == "NA"] <- NA_character_

# pred
tb <- rfPred(x_embed, 
             y,
             query_cells = NULL)

tb <- tibble(cell = rownames(sobj@meta.data)) %>%
  left_join(tb, by = "cell")
tb %>%
  write_tsv("cellRefine_knnLabel_rfPred.tsv")

tb_anno <- tb

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

pdf("cellRefine_knnLabel_rfPred_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "celltype", cols = celltype_cols)
dev.off()


# deg
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("cellRefine_knnLabel_rfPred_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("cellRefine_knnLabel_rfPred_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


sessionInfo()
