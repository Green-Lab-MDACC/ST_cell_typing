rm(list = ls())

library(tidyverse)
library(Seurat)


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
  "Epi", "#fb8072",
  "NA", "#7f7f7f" # "grey50"
) %>%
  deframe()
sorted_celltypes <- names(celltype_cols)


# load
tb <- read_tsv("singleR_refLM.tsv")

# merge Mac/Mo
for (x in c("labels", "pruned_labels")) {
  tb[[x]] <- ifelse(!is.na(tb[[x]]) & tb[[x]] == "Mo",
                    "Mac",
                    tb[[x]])
}

tb_ct <- tb

# epi
tb <- read_tsv("singleR_refBEL.tsv")
tb <- tb %>%
  filter(labels %in% c("Epi"))
if (nrow(tb) > 0) {
  tb_ct <- tb_ct %>%
    filter(!(cell %in% tb$cell)) %>%
    bind_rows(tb)
}

# write
tb_ct %>%
  write_tsv("init_label.tsv")

# umap
sobj <- readRDS("smi_data_prep.rds")

tb <- tb_ct %>%
  dplyr::rename(celltype = labels) %>%
  mutate(celltype = replace_na(celltype, "NA"))

df <- tibble(cell = colnames(sobj)) %>%
  left_join(tb, by = "cell") %>%
  select(-cell) %>%
  as.data.frame()
sobj@meta.data <- sobj@meta.data[, !(colnames(sobj@meta.data) %in% colnames(df))] %>%
  cbind(df)

stopifnot(length(setdiff(sobj$celltype, sorted_celltypes)) == 0)
sorted_celltypes <- sorted_celltypes[sorted_celltypes %in% sobj$celltype]
sobj$celltype <- factor(sobj$celltype, levels = sorted_celltypes)

pdf("init_label_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "celltype", cols = celltype_cols)
dev.off()

# deg
Idents(sobj) <- "celltype"
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE)
tb <- sobj.markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58)

tb %>%
  write_tsv("init_label_deg.tsv")

markers <- tb %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(p_val, desc(avg_log2FC)) %>%
  do(head(., n = 10)) %>%
  .$gene %>%
  unique()

pdf("init_label_deg_dot.pdf", height = 7, width = 30)
DotPlot(sobj, #group.by = "seurat_clusters",
        features = markers,
        cols = "RdBu") + # cols = "RdYlBu"
  RotatedAxis()
dev.off()


sessionInfo()
