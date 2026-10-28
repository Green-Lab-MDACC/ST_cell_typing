rm(list = ls())

library(tidyverse)
library(Seurat)
library(randomForest)
library(harmony)

# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^3)


sel_ct <- "NESC"


# load
sobj <- readRDS("scRNALymphoMAP_coembed.rds")
DefaultAssay(sobj) <- "integrated"

sobj_all <- sobj

# ref
sobj <- sobj_all[, sobj_all$dataset == "ref"]

n_dims <- 50
sobj <- RunPCA(sobj, npcs = n_dims, verbose = FALSE)

sel_celltypes <- sobj$celltype %>% unique() %>% sort()
sobj$celltype <- sobj$celltype %>%
  factor(levels = sel_celltypes)

# use Seurat default color (ggplot2 default color)
celltype_cols <- structure(scales::hue_pal()(length(sel_celltypes)),
                           names = c(sel_celltypes))

# harmony
sobj$batch_id <- sobj$sample_id
stopifnot(sum(is.na(sobj$batch_id)) == 0)

n_dims_umap_clust <- 50
minor_celltypes <- c("NK",
                     "Mo",
                     "DC",
                     "pDC",
                     "LEC")
if (sel_ct %in% minor_celltypes) {
  n_dims_umap_clust <- 15
}

k.param <- 50
if (length(unique(sobj$batch_id)) > 1) {
  sobj <- sobj %>%
    RunHarmony("batch_id",
               dims.use = 1:n_dims_umap_clust,
               max.iter.harmony = 100)
  sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:n_dims_umap_clust, return.model = TRUE, verbose = FALSE)
  sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:n_dims_umap_clust, k.param = k.param, verbose = FALSE)
} else {
  sobj <- RunUMAP(sobj, dims = 1:n_dims_umap_clust, return.model = TRUE, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:n_dims_umap_clust, k.param = k.param, verbose = FALSE)
}

# sobj <- FindClusters(sobj, resolution = 0.1, verbose = FALSE)

# save
sobj %>%
  saveRDS("scRNALymphoMAP_coembed_integRefHarmony.rds")


pdf("scRNALymphoMAP_coembed_integRefHarmony_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "celltype", cols = celltype_cols)
dev.off()

pdf("scRNALymphoMAP_coembed_integRefHarmony_umap_sample.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "sample_id") + NoLegend()
dev.off()


sobj_ref <- sobj


# labelTransfer
sobj_query <- sobj_all[, sobj_all$dataset == "query"]
# sobj_list <- as.list(c(sobj_ref, sobj_query))

n_dims <- 50
sobj_anchors <- FindTransferAnchors(reference = sobj_ref,
                                    query = sobj_query,
                                    normalization.method = "SCT",
                                    recompute.residuals = FALSE,
                                    reference.assay = "integrated",
                                    query.assay = "integrated",
                                    scale = FALSE,
                                    dims = 1:n_dims,
                                    reference.reduction = "pca")
sobj_query <- MapQuery(anchorset = sobj_anchors,
                       reference = sobj_ref,
                       query = sobj_query,
                       refdata = list(celltype = "celltype"),
                       reference.reduction = "pca",
                       reduction.model = "umap")

pdf("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_sep_umap.pdf", height = 5, width = 14)
p1 <- DimPlot(sobj_ref,
              reduction = "umap",
              group.by = "celltype",
              cols = celltype_cols,
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
  # NoLegend() + 
  ggtitle("Reference annotations")
p2 <- DimPlot(sobj_query,
              reduction = "ref.umap",
              group.by = "predicted.celltype",
              cols = celltype_cols,
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
  # NoLegend() + 
  ggtitle("Query transferred labels")
p1 + p2
dev.off()


# transferred umap
mt <- sobj_query@reductions$ref.umap@cell.embeddings

# write
tibble(cell = rownames(mt)) %>%
  bind_cols(as_tibble(mt)) %>%
  write_tsv("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_umap.tsv")


# query
tb <- tibble(cell = colnames(sobj_query),
             celltype = sobj_query$predicted.celltype,
             celltype_score = sobj_query$predicted.celltype.score)
tb %>%
  write_tsv("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_celltype.tsv")


# markers ------
# ref
sobj <- sobj_ref

# deg
DefaultAssay(sobj) <- "integrated"
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_ref_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_ref_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


# SCT
DefaultAssay(sobj) <- "SCT"

# deg
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_ref_sct_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_ref_sct_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


# query
sobj <- sobj_query
sobj$celltype <- sobj$predicted.celltype

stopifnot(all.equal(
  sobj$celltype %>% unique() %>% as.character() %>% sort(),
  sobj_ref$celltype %>% unique() %>% as.character() %>% sort()
))
sobj$celltype <- factor(sobj$celltype,
                        levels = levels(sobj_ref$celltype))

# deg
DefaultAssay(sobj) <- "integrated"
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_query_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_query_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


# SCT
DefaultAssay(sobj) <- "SCT"

# deg
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_query_sct_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("scRNALymphoMAP_coembed_integRefHarmony_labelTransfer_query_sct_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


sessionInfo()
