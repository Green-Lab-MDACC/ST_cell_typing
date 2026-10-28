# scRNA-seq Anchored Cell Typing for Spatial Transcriptomics

A scalable, hierarchical reference-anchoring framework designed for high-resolution cell-type annotation of high-plex targeted spatial transcriptomics data (e.g., **CosMx Spatial Molecular Imager [SMI]**) using single-nucleus or single-cell RNA sequencing (**sn/scRNA-seq**) references.

---

## Overview

Compared to sc/snRNA-seq, which profiles whole cells or nuclei via untargeted sequencing with molecular amplification, targeted in situ imaging platforms like CosMx SMI detect fewer transcripts per cell and are more susceptible to lateral transcript misassignment caused by cell segmentation errors in dense tissue microenvironments. These technical factors compress cell populations in gene expression space, making de novo unsupervised clustering and cell-type annotation challenging.

To address this, this framework performs **hierarchical, supervised reference co-embedding and transfer mapping**, leveraging matched or atlas-scale sn/scRNA-seq data with high-resolution cell annotations (e.g., large B-cell lymphoma snRNAseq atlas; [PMID: 40920660](https://pubmed.ncbi.nlm.nih.gov/40920660/)).

```
         Step 1: Initial Cell Label Preparation (SingleR + KNN)
                      Raw Spatial / SMI Data
                                │
                   ┌────────────┴────────────┐
                   ▼                         ▼
         Blueprint/ENCODE SingleR      Custom snRNA-seq SingleR
                   │                         │
         Filter Out-of-Reference             │
         Lineages (e.g., Epithelial)         │
                   │                         │
                   └────────────┬────────────┘
                                ▼
         Refine Labels with KNN (k=50, ≥75% agreement)
                                │
                    High-Confidence SMI Labels
                                │
                                ▼
         Step 2: Lineage Separation (RPCA + Random Forest)
         Features: Top 100 DEGs per cell type 
         Embedding: RPCA
         Model: Random Forest 
  ──────┬──────────────┬───────────────┬──────────────┬────────
        ▼              ▼               ▼              ▼
     B Cells        T / NK          Myeloid          NHCs
        │              │               │              │
        │              └───────┬───────┴──────────────┘
        │                      ▼
        │        Step 3: Cell Type Classification (CCA + MapQuery)
        │        • T/NK: CD4 T, CD8 T, NK
        │        • Myeloid: Macrophages, Monocytes, DCs, pDCs
        │        • NHC: NESCs, BECs, LECs
        │                      │
        │                      ▼
        │        Step 4: Cell Subset & State Classification (CCA + MapQuery)
        │        • Finer-resolution subsets (e.g., CD4 T C0–C7)
        ▼                      ▼
 ─────────────────────────────────────────────────────────────
              Final High-Resolution Spatial Atlas
```

---

## Key Features

- **Out-of-Reference Lineage Filtering:** Identifies and separates tissue-specific cell types absent from the single-cell reference (e.g., epithelial elements) prior to integration.
- **Hierarchical Gating Strategy:** Progresses from major lineages to cell types, and finally to fine-grained cell states/subsets, optimizing the feature space at each tier.
- **Strategic Integration Anchoring:**
  - **Reciprocal PCA (RPCA):** Preserves biological heterogeneity and malignant-specific expression profiles during major lineage separation.
  - **Canonical Correlation Analysis (CCA):** Enforces alignment across technologies for shared non-malignant microenvironment populations.
- **Semi-Supervised Label Propagation:** Augments reference training data with high-confidence spatial cells verified via local neighborhood consensus ($k$-NN), capturing assay-specific technical signals.

---

## Workflow Breakdown

### Step 1: Out-of-Reference Lineage Identification & Initial Label Assignment
1. **Unconstrained Query Profiling:** Predict initial cell identities across all spatial cells with [`SingleR`] (`v2.0.0`) using the broad human **Blueprint/ENCODE** reference ([PMID: 30643263](https://pubmed.ncbi.nlm.nih.gov/30643263/)).
2. **Exclude Out-of-Reference Classes:** Remove lineages absent in the reference snRNA-seq (e.g., tonsillar epithelial cells).
3. **Custom Reference Mapping:** Generate pseudobulk profiles from reference snRNA-seq populations ($N \ge 100$ cells) into a custom `celldex` reference to predict spatial cell identities.
4. **$k$-Nearest Neighbor (KNN) Smoothing:** Refine predicted labels in PCA space ($k = 50$). Cells where $\ge 75%$ of neighbors share the identical predicted label are designated as **high-confidence SMI anchors**.

### Step 2: Major Lineage Separation
1. **Feature Space Optimization:** Subset the feature space to the top 100 differentially expressed genes (DEGs; adjusted $P < 0.05$, $\text{FC} > 1.5$) per cell type from the snRNA-seq reference.
2. **RPCA Co-Embedding:** Normalize both datasets using `SCTransform`, prepare integration with `PrepSCTIntegration()`, and compute the top 50 PCs. Identify mutual nearest neighbors via `FindIntegrationAnchors(reduction = "rpca")` and align data via `IntegrateData()`.
3. **Random Forest Classification:** Train a Random Forest classifier ($n_{\text{tree}} = 500$) using the top 50 integrated PCs on the combined reference pool (snRNA-seq ground truth + high-confidence SMI anchors).
4. **Hierarchical Splitting:** Predict lineages for all spatial cells and partition data into four primary lineages:
   - **B cells**
   - **T / NK cells**
   - **Myeloid cells**
   - **Non-hematopoietic cells (NHCs)**

### Step 3: Cell Type Classification Within Lineages
For non-B cell lineages (T/NK, Myeloid, NHC):
1. **Lineage-Specific DEGs:** Isolate top 100 DEGs distinguishing intra-lineage cell types:
   - *T/NK:* CD4 T cells, CD8 T cells, NK cells.
   - *Myeloid:* Macrophages, Monocytes, Dendritic Cells (DCs), Plasmacytoid DCs (pDCs).
   - *NHC:* Non-endothelial stromal cells (NESCs), Blood endothelial cells (BECs), Lymphatic endothelial cells (LECs).
2. **CCA Integration:** Perform `FindIntegrationAnchors(reduction = "cca")` specifying the snRNA-seq dataset as reference to map microenvironment populations across technologies.
3. **Batch Correction & Label Transfer:** Run `Harmony` on reference cells to correct sample-level variance, establish a reference UMAP, and transfer discrete cell-type labels and coordinate projections to spatial cells using `FindTransferAnchors()` and `MapQuery()`.

### Step 4: Fine-Grained Subset & State Classification
1. **Subpopulation Dissection:** Recapitulate Step 3 at deep resolution (e.g., CD4 T cell subsets C0–C7).
2. **State-Level Marker Sets:** Re-compute top 100 DEGs differentiating intra-cell-type phenotypes.
3. **CCA Mapping & Transfer:** Execute reference-constrained CCA co-embedding and project final subset calls and latent coordinates onto query cells.

---

## Scalability via Chunked Partitioning

To handle spatial datasets scaling to millions of cells without exceeding system RAM in the Seurat workflow:
- Cells are shuffled randomly at each tier and partitioned into parallel blocks of 200,000 cells.
- The pipeline was validated across chunk ranges of 100,000 to 500,000 cells, showing consistent annotation concordance across independently processed partitions.
- Annotated partitions are concatenated prior to the next downstream hierarchical gate.

---

## Dependencies & Environment

Key dependencies include:

- [Seurat](https://satijalab.org/seurat/) (v5)
- [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html) ($\ge 2.0.0$)
- [celldex](https://bioconductor.org/packages/release/data/experiment/html/celldex.html)
- [randomForest](https://cran.r-project.org/web/packages/randomForest/)
- [harmony](https://github.com/immunogenomics/harmony)
