# Set your working directory to the root of this project before running.

library(Seurat)
library(DropletUtils)
library(scDblFinder)
library(SingleCellExperiment)
library(tidyverse)

# Read a 10X filtered matrix from either a directory or an .h5 file.
# Uses DropletUtils::read10xCounts() for h5 so hdf5r is not required.
read_filtered_counts <- function(filtered_path, h5_path) {
  mtx_file <- file.path(filtered_path, "matrix.mtx.gz")
  if (file.exists(mtx_file)) {
    return(Seurat::Read10X(filtered_path))
  } else if (file.exists(h5_path)) {
    sce <- DropletUtils::read10xCounts(h5_path)
    m   <- SingleCellExperiment::counts(sce)
    colnames(m) <- sce$Barcode
    return(m)
  }
  stop("No filtered matrix found at ", filtered_path, " or ", h5_path)
}

# --- CONFIG ---
# DATA_ROOT: parent directory containing one subdirectory per sample.
#   Each sample subdirectory must contain an outs/ folder produced by Cell Ranger.
DATA_ROOT <- '/path/to/cellranger/outputs'
# --- END CONFIG ---

data_dirs <- list(
  Ctrl13 = file.path(DATA_ROOT, 'CTRL13', 'outs'),
  Ctrl24 = file.path(DATA_ROOT, 'CTRL24', 'outs'),
  Ctrl25 = file.path(DATA_ROOT, 'CTRL25', 'outs'),
  Ctrl26 = file.path(DATA_ROOT, 'CTRL26', 'outs'),
  Ctrl27 = file.path(DATA_ROOT, 'CTRL27', 'outs'),
  Sal7   = file.path(DATA_ROOT, 'SAL7',   'outs'),
  Sal9   = file.path(DATA_ROOT, 'SAL9',   'outs'),
  Sal10  = file.path(DATA_ROOT, 'SAL10',  'outs'),
  Sal13  = file.path(DATA_ROOT, 'SAL13',  'outs')
)

# =============================================================================
# Stage 1 & 2: CellRanger filtered counts + emptyDrops per sample
# Also builds per-sample Seurat objects (basic QC only, no SoupX) for merging.
# SoupX is skipped — it corrects counts but never removes cells.
# =============================================================================

sample_stats <- list()
seurat_list  <- list()

for (sample_name in names(data_dirs)) {
  cat("Processing", sample_name, "...\n")
  sample_path   <- data_dirs[[sample_name]]
  raw_path      <- file.path(sample_path, "raw_feature_bc_matrix")
  filtered_path <- file.path(sample_path, "filtered_feature_bc_matrix")
  h5_path       <- file.path(sample_path, "filtered_feature_bc_matrix.h5")
  mtx_file      <- file.path(filtered_path, "matrix.mtx.gz")

  tryCatch({
    raw_counts <- Seurat::Read10X(raw_path)
    n_raw <- ncol(raw_counts)

    filtered_counts <- read_filtered_counts(filtered_path, h5_path)
    n_cellranger <- ncol(filtered_counts)

    # emptyDrops
    set.seed(42)
    e.out   <- DropletUtils::emptyDrops(raw_counts)
    is_cell <- which(!is.na(e.out$FDR) & e.out$FDR < 0.01)
    if (length(is_cell) >= 100) {
      cell_counts <- raw_counts[, is_cell]
    } else {
      cat("  emptyDrops called too few cells — using CellRanger filtered matrix\n")
      cell_counts <- filtered_counts
    }
    n_emptyd <- ncol(cell_counts)

    # Basic per-sample QC (nFeature > 200, MT < 5%) — mirrors main pipeline
    fb <- CreateSeuratObject(counts = cell_counts, project = sample_name, min.cells = 1)
    fb$sample <- sample_name
    fb[["percent.mt"]] <- PercentageFeatureSet(fb, pattern = "^MT-")
    fb <- subset(fb, subset = nFeature_RNA > 200 & percent.mt < 5)

    sample_stats[[sample_name]] <- list(
      n_raw        = n_raw,
      n_cellranger = n_cellranger,
      n_emptyd     = n_emptyd,
      n_basic_qc   = ncol(fb)
    )
    seurat_list[[sample_name]] <- fb
    cat("  CellRanger:", n_cellranger,
        "| emptyDrops:", n_emptyd,
        "| after basic QC:", ncol(fb), "\n")

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    sample_stats[[sample_name]] <<- list(n_cellranger = NA, n_emptyd = NA, n_basic_qc = NA)
  })
}

# =============================================================================
# Stage 3: Merge and run scDblFinder to count doublets per sample
# Uses strict QC thresholds from the main pipeline (nFeature 200–7500, MT < 2.5%).
# =============================================================================

sample_ids <- names(seurat_list)
fb_merged  <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[2:length(seurat_list)],
  add.cell.ids = sample_ids,
  project      = "all_samples"
)
fb_merged <- JoinLayers(fb_merged)

cat("\nRunning scDblFinder on merged object (", ncol(fb_merged), "cells) ...\n")
set.seed(42)
sce        <- as.SingleCellExperiment(fb_merged)
sce        <- scDblFinder(sce, samples = "orig.ident")
fb_merged$scDblFinder_class <- sce$scDblFinder.class
rm(sce)

# Doublets per sample
doublet_counts <- table(fb_merged$sample[fb_merged$scDblFinder_class == "doublet"])
rm(fb_merged)

# =============================================================================
# Stage 4: Final annotated object
# =============================================================================

cat("\nLoading fb_seurat_FINAL_v2.RDS ...\n")
fb_seurat        <- readRDS("fb_seurat_FINAL_v2.RDS")
final_counts     <- table(fb_seurat$sample)
cell_type_counts <- table(cell_type_v8 = fb_seurat$cell_type_v8)

# =============================================================================
# Table 1: Pipeline cell counts per sample
# =============================================================================

samples <- names(data_dirs)

pipeline_table <- tibble(
  Sample            = samples,
  Raw_Barcodes      = sapply(samples, function(s) sample_stats[[s]]$n_raw),
  CellRanger_Called = sapply(samples, function(s) sample_stats[[s]]$n_cellranger),
  emptyDrops_Called = sapply(samples, function(s) sample_stats[[s]]$n_emptyd),
  After_Basic_QC    = sapply(samples, function(s) sample_stats[[s]]$n_basic_qc),
  Doublets_Removed  = as.integer(doublet_counts[samples]),
  Final             = as.integer(final_counts[samples])
) %>%
  mutate(Pct_Retained = round(Final / emptyDrops_Called * 100, 1))

totals <- tibble(
  Sample            = "Total",
  Raw_Barcodes      = sum(pipeline_table$Raw_Barcodes,      na.rm = TRUE),
  CellRanger_Called = sum(pipeline_table$CellRanger_Called, na.rm = TRUE),
  emptyDrops_Called = sum(pipeline_table$emptyDrops_Called, na.rm = TRUE),
  After_Basic_QC    = sum(pipeline_table$After_Basic_QC,    na.rm = TRUE),
  Doublets_Removed  = sum(pipeline_table$Doublets_Removed,  na.rm = TRUE),
  Final             = sum(pipeline_table$Final,              na.rm = TRUE),
  Pct_Retained      = round(sum(pipeline_table$Final,        na.rm = TRUE) /
                             sum(pipeline_table$emptyDrops_Called, na.rm = TRUE) * 100, 1)
)

pipeline_table <- bind_rows(pipeline_table, totals)

write_csv(pipeline_table, "qc_pipeline_table.csv")
cat("\nPipeline table saved to qc_pipeline_table.csv\n")
print(pipeline_table, width = Inf)

# =============================================================================
# Table 2: Cell counts per cell type in the final object
# =============================================================================

cell_type_table <- as_tibble(cell_type_counts) %>%
  rename(Cell_Type = cell_type_v8, Count = n) %>%
  mutate(Percent = round(Count / sum(Count) * 100, 2)) %>%
  arrange(desc(Count)) %>%
  bind_rows(tibble(Cell_Type = "Total",
                   Count     = sum(.$Count),
                   Percent   = 100.00))

write_csv(cell_type_table, "cell_type_counts_final.csv")
cat("\nCell type table saved to cell_type_counts_final.csv\n")
print(cell_type_table, n = Inf)
