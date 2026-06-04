# Fetal Brain SVZ Single-Cell Atlas

Single-cell RNA-seq atlas of the fetal macaque subventricular zone (SVZ) at late gestation. Nine samples (5 control, 4 saline-treated) processed through ambient RNA correction, doublet detection, Harmony batch integration, and iterative cell type annotation. Trajectory analysis performed with scVelo and CellRank 2.

## Samples

| ID | Group | Gestational age |
|----|-------|----------------|
| Ctrl13, Ctrl24, Ctrl25, Ctrl26, Ctrl27 | Control | 144–158 days |
| Sal7, Sal9, Sal10, Sal13 | Saline | 128–135 days |

## Cell types (v8 annotations)

| Lineage | Cell types |
|---------|-----------|
| Radial glia | vRG, oRG-G |
| Progenitors | APC, UP/I, IPC |
| Astrocytes | AS0, AS1, AS2 |
| Oligodendrocyte lineage | OPC, COP1, COP2, OL1, OL2 |
| Neurons | ExN, ImmN, CGE, MGE |
| Microglia | MG0, MG1, BAM |
| Non-neural | ChP, Ep, EC, FIB, PC |

## Scripts

### R

| Script | Purpose |
|--------|---------|
| `fb_v4.R` | **Main pipeline.** Per-sample QC (emptyDrops → SoupX → scDblFinder), merge, SCTransform normalization, Harmony batch correction, label transfer from reference atlas, iterative cluster annotation (v1–v8), lineage subset re-clustering, scVelo metadata export. |
| `make_qc_tables.R` | Regenerates QC summary tables from raw Cell Ranger output + the final Seurat object: pipeline cell counts per sample and cell-type composition of the final atlas. |
| `atlas_palette.R` | Shared color palette (`ATLAS_PALETTE`, `SUBSET_PALETTE`). Defines `get_atlas_colors()` and `get_subset_colors()` helpers. Identical colors are used in all R and Python figures. |
| `plot_umap_clusters.R` | Wrapper around `DimPlot_scCustom` + ggrepel. Generates labeled UMAP cluster plots with shadow text at cluster centroids. |
| `plot_umap_genes.R` | Overlays multi-gene expression on a UMAP — expressing cells plotted over a grey background, each gene a distinct color. |
| `plot_harmony_integration.R` | Publication-ready Harmony QC figure: PCA/UMAP pre- and post-correction colored by sample and gestational age, plus iLISI violin plots computed in pure R (no external package). |

### Python

| Script | Purpose |
|--------|---------|
| `run_scvelo_brain_full_3.py` | RNA velocity (stochastic mode) on the full atlas. Loads per-sample loom files, attaches Seurat UMAP and metadata, writes plots and an h5ad checkpoint for CellRank. |
| `run_scvelo_brain_subsets_5.py` | Same as above but run independently on each lineage subset (astrocyte, OL, microglia, neuron, excitatory neuron) using subset-specific re-embedded UMAPs. |
| `run_cellrank_brain_5.py` | CellRank 2 fate analysis on scVelo h5ad outputs. Runs velocity-directed PAGA, identifies terminal/initial states (GPCCA estimator), computes per-cell fate probabilities, identifies lineage driver genes, and plots gene expression trends along pseudotime. |
| `plot_mg_specific_genes.py` | Standalone script to generate pseudotime expression trend plots for specific microglia marker genes (IL1B, CD86) for MG0/MG1 clusters. Injects genes missing from the HVG-filtered h5ad directly from loom files. |

## Data

`fb_seurat_FINAL_v2.RDS` — Final annotated Seurat object (v8 cell type labels in `cell_type_v8`).

## Dependencies

**R:** Seurat, DropletUtils, SoupX, scDblFinder, SingleCellExperiment, harmony, SCTransform, scCustomize, SeuratExtend, dittoSeq, sccomp, FNN, ggrepel, shadowtext, patchwork, tidyverse

**Python:** scvelo, cellrank (≥2.x), scanpy, anndata, scipy, numpy, pandas, matplotlib

## Usage

1. Set `DATA_ROOT` (Cell Ranger output directory) and `ATLAS_REF` (reference Seurat object path) at the top of `fb_v4.R`.
2. Run `fb_v4.R` section by section. The final object is saved as `fb_seurat_FINAL_v2.RDS`.
3. Export per-subset metadata CSVs using `save_subset_for_scvelo()`, then run `run_scvelo_brain_full_3.py` and `run_scvelo_brain_subsets_5.py`.
4. Run `run_cellrank_brain_5.py` using the h5ad checkpoints from step 3.
5. Generate QC tables with `make_qc_tables.R` and integration figures with `plot_harmony_integration.R`.
