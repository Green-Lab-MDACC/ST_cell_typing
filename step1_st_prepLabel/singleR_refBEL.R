rm(list = ls())

library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(SingleR)


# load data
sobj <- readRDS("smi_data_prep.rds")
sce <- as.SingleCellExperiment(sobj)

# sce <- readRDS(sce_path)
sce <- scuttle::logNormCounts(sce) 

# run singleR to get BlueprintEncode grouped lineages
ref <- readRDS("../../SingleR_ref/BlueprintEncode_grouped_lineage_celldex/BlueprintEncodeData_group_250610.rds")
pred <- SingleR(test = sce, 
                ref = ref, 
                labels = ref$label.main)
pred$labels %>% table() %>% sort(decreasing = TRUE) %>% print()
summary(is.na(pred$pruned.labels)) %>% print()

tb <- tibble(cell = rownames(pred),
             spearman_corr_score = apply(pred$scores, 1, max),
             labels = pred$labels,
             pruned_labels = pred$pruned.labels) 
tb %>%
  write_tsv("singleR_refBEL.tsv")

# plot
stopifnot(all.equal(tb$cell, colnames(sobj)))
sobj$labels <- tb$labels

pdf("singleR_refBEL_umap.pdf", height = 5, width = 7)
DimPlot(sobj, group.by = "labels")
dev.off()


sessionInfo()
