# =============================================================================
# 1. SETUP
# =============================================================================

# Set your working directory to the root of this project before running.

library(Seurat)
library(DropletUtils)
library(SoupX)
library(scDblFinder)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(scCustomize)
library(harmony)
library(future)
library(SeuratExtend)
library(SeuratDisk)
library(dittoSeq)
library(FNN)
library(ggrepel)
library(shadowtext)
library(gatepoints)
library(shiny)
library(plotly)
library(sccomp)
library(patchwork)

source('atlas_palette.R')
source('plot_umap_clusters.R')
source('plot_umap_genes.R')
source('plot_harmony_integration.R')

gradient <- c("blue", "white", "red")

set.seed(42)
options(future.globals.maxSize = 1e20)
fb_seurat <- readRDS("fb_seurat_FINAL_v2.RDS")

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

# ---- save UMAP function ----
save_umap <- function(plot, filename, width = 15) {
  plot <- plot +
    coord_fixed(ratio = 0.75) +
    theme(aspect.ratio = 0.75)

  ggsave(
    filename,
    plot = plot,
    width = width,
    height = width * 0.75,
    device = "pdf"
  )
}

# Function to save subset metadata and UMAP coordinates
save_subset_for_scvelo <- function(seurat_obj, name) {
  metadata <- seurat_obj@meta.data
  metadata$barcode <- rownames(metadata)
  
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  
  # Fix: rename to lowercase to match Python scripts
  colnames(umap_coords)[colnames(umap_coords) == "UMAP_1"] <- "umap_1"
  colnames(umap_coords)[colnames(umap_coords) == "UMAP_2"] <- "umap_2"
  
  umap_coords$barcode <- rownames(umap_coords)
  
  combined <- merge(metadata, umap_coords, by = "barcode")
  
  write.csv(combined, paste0(name, "_metadata.csv"), row.names = FALSE)
  cat("Saved metadata for", name, "with", nrow(combined), "cells\n")
}

# Interactive Plotly lasso for manual cell gating.
# Saves selected barcodes to gated_cells.rds; load and relabel after running.
run_lasso <- function(seurat_obj, reduction = "umap") {
  umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction))
  colnames(umap_coords) <- tolower(colnames(umap_coords))
  cluster_labels <- as.character(Idents(seurat_obj))

  ui <- fluidPage(
    plotlyOutput("umap", height = "600px"),
    verbatimTextOutput("selected"),
    actionButton("save", "Save Selected Cells")
  )

  server <- function(input, output, session) {
    output$umap <- renderPlotly({
      plot_ly(umap_coords, x = ~umap_1, y = ~umap_2,
              type = "scattergl", mode = "markers",
              marker = list(size = 3),
              text = cluster_labels,
              hovertemplate = "Cluster: %{text}<extra></extra>",
              source = "umap") %>%
        layout(dragmode = "lasso")
    })

    selected_cells <- reactive({
      event_data("plotly_selected", source = "umap")
    })

    output$selected <- renderPrint({
      req(selected_cells())
      paste("Selected:", nrow(selected_cells()), "cells")
    })

    observeEvent(input$save, {
      req(selected_cells())
      cells <- rownames(umap_coords)[selected_cells()$pointNumber + 1]
      saveRDS(cells, "gated_cells.rds")
      showNotification("Cells saved to gated_cells.rds")
    })
  }

  shinyApp(ui, server)
}

# =============================================================================
# 3. DATA LOADING
# =============================================================================

# --- CONFIG ---
# DATA_ROOT: parent directory containing one subdirectory per sample.
#   Each sample subdirectory must contain an outs/ folder produced by Cell Ranger.
#   Example structure:  DATA_ROOT/CTRL13/outs/  DATA_ROOT/SAL7/outs/  ...
# ATLAS_REF: path to the reference Seurat object used for cell type label transfer.
DATA_ROOT <- '/path/to/cellranger/outputs'
ATLAS_REF <- '/path/to/reference_atlas.RDS'
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
# 4. PER-SAMPLE QC & AMBIENT RNA CORRECTION
# =============================================================================

# --- Process each sample individually ---
fb_seurat_list <- list()
for (sample_name in names(data_dirs)) {
  cat("Processing sample:", sample_name, "\n")

  sample_path <- data_dirs[[sample_name]]

  # Check if the required directories exist
  raw_path <- file.path(sample_path, "raw_feature_bc_matrix")
  filtered_path <- file.path(sample_path, "filtered_feature_bc_matrix")

  if (!dir.exists(raw_path) || !dir.exists(filtered_path)) {
    cat("Warning: Required matrix paths not found for", sample_name, "- skipping\n\n")
    next
  }

  # Try to load and process the data with error handling
  tryCatch({
    # Load raw and filtered data
    # Load raw counts (always from directory)
    cat("  Loading counts...\n")
    raw_counts <- Seurat::Read10X(raw_path)

    # Load filtered counts — use .h5 if matrix file is missing
    filtered_mtx <- file.path(filtered_path, "matrix.mtx.gz")
    h5_path <- file.path(sample_path, "filtered_feature_bc_matrix.h5")

    if (!file.exists(filtered_mtx) && file.exists(h5_path)) {
      cat("  matrix.mtx.gz not found — loading filtered counts from .h5 file\n")
      filtered_counts <- Seurat::Read10X_h5(h5_path)
    } else {
      filtered_counts <- Seurat::Read10X(filtered_path)
    }

    # --- NEW: emptyDrops filtering ---
    cat("  Running emptyDrops...\n")
    set.seed(42)
    e.out <- DropletUtils::emptyDrops(raw_counts)

    # Check how many cells are called
    is_cell <- which(e.out$FDR < 0.01)
    cat("  emptyDrops called", length(is_cell), "cells (vs", ncol(filtered_counts), "from CellRanger)\n")

    if (length(is_cell) < 100) {
      cat("  Too few cells from emptyDrops, falling back to CellRanger filtered matrix\n")
    } else {
      filtered_counts <- raw_counts[, is_cell]  # replace with emptyDrops calls
    }
    # --- end emptyDrops section ---


    # Check if data is valid
    if (ncol(filtered_counts) == 0 || nrow(filtered_counts) == 0) {
      cat("Warning: No cells or features in filtered matrix for", sample_name, "- skipping\n\n")
      next
    }

    cat("  Sample has", ncol(filtered_counts), "cells and", nrow(filtered_counts), "features\n")

    # Decide whether to use SoupX based on cell count
    use_soupx <- ncol(filtered_counts) >= 100  # Only use SoupX if >= 100 cells

    if (!use_soupx) {
      cat("  Too few cells for SoupX (", ncol(filtered_counts), " cells). Skipping SoupX correction.\n")
    }

    soupx_success <- FALSE
    corrected_counts <- NULL

    if (use_soupx) {
      # Create preliminary Seurat object for clustering
      cat("  Creating preliminary Seurat object...\n")
      fb_temp <- CreateSeuratObject(counts = filtered_counts, min.cells = 1, min.features = 1)
      fb_temp <- NormalizeData(fb_temp, verbose = FALSE)
      fb_temp <- FindVariableFeatures(fb_temp, verbose = FALSE)
      fb_temp <- ScaleData(fb_temp, verbose = FALSE)
      fb_temp <- RunPCA(fb_temp, verbose = FALSE)
      fb_temp <- FindNeighbors(fb_temp, dims = 1:10, verbose = FALSE)
      fb_temp <- FindClusters(fb_temp, resolution = 0.5, verbose = FALSE)
      prelim_clusters <- fb_temp$seurat_clusters

      # --- Try SoupX correction ---
      cat("  Attempting SoupX correction...\n")

      tryCatch({
        sc <- SoupChannel(raw_counts, filtered_counts)
        sc <- setClusters(sc, prelim_clusters)
        # Use more lenient parameters
        tryCatch({
          sc <- autoEstCont(sc, tfidfMin = 0.5, soupQuantile = 0.9)
        }, error = function(e) {
          cat("  Default SoupX params failed, retrying with lenient params...\n")
          sc <<- autoEstCont(sc, tfidfMin = 0.1, soupQuantile = 0.7)
        })
        corrected_counts <- adjustCounts(sc, roundToInt = TRUE)
        soupx_success <- TRUE
        cat("  SoupX correction successful!\n")
      }, error = function(e) {
        cat("  SoupX failed:", conditionMessage(e), "\n")
        cat("  Proceeding without SoupX correction...\n")
      })
    }

    # --- Create Seurat object (with or without SoupX correction) ---
    if (soupx_success) {
      fb <- CreateSeuratObject(counts = corrected_counts, project = sample_name)
    } else {
      fb <- CreateSeuratObject(counts = filtered_counts, project = sample_name)
    }

    fb$sample <- sample_name
    fb$soupx_corrected <- soupx_success  # Track which samples were corrected

    # --- QC for the current sample ---
    fb[["percent.mt"]] <- PercentageFeatureSet(fb, pattern = "^MT-")

    # Check if any cells pass QC
    cells_before <- ncol(fb)
    fb <- subset(fb, subset = nFeature_RNA > 200 & percent.mt < 5)
    cells_after <- ncol(fb)

    if (cells_after == 0) {
      cat("Warning: No cells passed QC filters for", sample_name, "- skipping\n\n")
      next
    }

    # Add to the list
    fb_seurat_list[[sample_name]] <- fb

    cat("Completed:", sample_name, "with", cells_after, "cells (", cells_before - cells_after, "filtered out)",
        ifelse(soupx_success, "(SoupX corrected)", "(no SoupX)"), "\n\n")

  }, error = function(e) {
    cat("ERROR processing", sample_name, ":\n")
    cat("  ", conditionMessage(e), "\n\n")
  })
}

# Check how many samples were successfully processed
cat("\n========================================\n")
cat("Successfully processed", length(fb_seurat_list), "out of", length(data_dirs), "samples\n")

# Summary of SoupX correction status
if (length(fb_seurat_list) > 0) {
  soupx_status <- sapply(fb_seurat_list, function(x) x$soupx_corrected[1])
  cat("SoupX corrected:", sum(soupx_status), "samples\n")
  cat("Not corrected:", sum(!soupx_status), "samples\n")

  # Show cell counts per sample
  cat("\nCell counts per sample:\n")
  for (sname in names(fb_seurat_list)) {
    cat("  ", sname, ":", ncol(fb_seurat_list[[sname]]), "cells\n")
  }
}

sample_ids <- names(fb_seurat_list)

# =============================================================================
# 5. MERGE, DOUBLET DETECTION & CELL QUALITY FILTERING
# =============================================================================

fb_merged <- merge(x = fb_seurat_list[[1]],
                         y = fb_seurat_list[2:length(fb_seurat_list)],
                         add.cell.ids = sample_ids,
                         project = "all_samples")

#DOUBLET DETECTION WITH scDblFinder
fb_merged <- JoinLayers(fb_merged)
# scDblFinder requires a SingleCellExperiment object
sce <- as.SingleCellExperiment(fb_merged)

sce <- scDblFinder(sce, samples = "orig.ident")

# Add doublet results back to Seurat metadata
fb_merged$scDblFinder_class <- sce$scDblFinder.class
fb_merged$scDblFinder_score <- sce$scDblFinder.score

# Get MT percentages
fb_merged <- PercentageFeatureSet(fb_merged, pattern = "^MT", col.name = "percent.mt")

VlnPlot(fb_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "sample")

ggsave("qc_plots_pre_filt.pdf", width = 30, height = 10)

plot1 <- FeatureScatter(fb_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fb_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Subset based on manual thresholds and doublet results
# Note: thresholds should be adjusted based on your specific dataset QC plots
fb_merged_filt <- subset(fb_merged,
                       subset = nFeature_RNA > 200 &
                                nFeature_RNA < 7500 &
                                percent.mt < 2.5 &
                                scDblFinder_class == "singlet")

VlnPlot(fb_merged_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "sample")

ggsave("qc_plots.pdf", width = 30, height = 10)

plot1 <- FeatureScatter(fb_merged_filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fb_merged_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# =============================================================================
# 6. NORMALIZATION, DIMENSIONALITY REDUCTION & HARMONY INTEGRATION
# =============================================================================

# --- SCTransform normalization ---
fb_merged_filt <- SCTransform(fb_merged_filt, vars.to.regress = "percent.mt")

# --- Perform PCA on the merged object ---
fb_merged_filt <- RunPCA(fb_merged_filt, assay = "SCT", npcs = 50)

DimPlot(fb_merged_filt, reduction = "pca", group.by = "sample")
ggsave("pca_pre_harmony_sample.pdf", width = 15, height = 10)

# --- Generate UMAP before Harmony integration ---
fb_merged_filt <- RunUMAP(fb_merged_filt, reduction = "pca", dims = 1:30,
                             reduction.name = "umap_preharmony",
                             reduction.key = "UMAPpreharmony_")
p <- DimPlot(fb_merged_filt, reduction = "umap_preharmony", group.by = "sample", shuffle = TRUE)
save_umap(p, "umap_preharmony.pdf")

saveRDS(fb_merged_filt, file = "fb_merged_pre_harmony.RDS")

# Run Harmony correcting by sample
fb_merged_filt <- RunHarmony(
  object = fb_merged_filt,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = 1:30,
  theta = 2,
  max.iter.harmony = 20,
  verbose = TRUE
)

# Generate UMAP on harmony-corrected embeddings
fb_merged_filt <- RunUMAP(
  fb_merged_filt,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3
)

# Clustering
fb_merged_filt <- FindNeighbors(fb_merged_filt, reduction = "harmony", dims = 1:30)
fb_merged_filt <- FindClusters(fb_merged_filt, resolution = 0.5)

# Visualize results
p1 <- DimPlot(fb_merged_filt, reduction = "umap", group.by = "sample") +
  ggtitle("By Sample - After Harmony")

save_umap(p1, "post_harmony_umap.pdf")

ggsave("pca_post_harmony_sample.pdf", width = 15, height = 10)

# =============================================================================
# 7. CELL TYPE LABEL TRANSFER
# =============================================================================

# Load your reference atlas (path set in CONFIG block above)
atlas <- readRDS(ATLAS_REF)


atlas <- SCTransform(atlas, vst.flavor = "v2", verbose = TRUE)


# Find transfer anchors between your query and reference
transfer_anchors <- FindTransferAnchors(
  reference = atlas,
  query = fb_merged_filt,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)

# Transfer the cell type labels
predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = atlas$cell_type_v3,
  dims = 1:30
)

# Add predictions to your object
fb_merged_filt <- AddMetaData(fb_merged_filt, metadata = predictions)

# Rename to match your previous naming convention
fb_merged_filt$predicted.celltype <- fb_merged_filt$predicted.id


# Visualize
  ggtitle("Transferred Cell Types - After Harmony")
ggsave("umap_post_harmony_xfer.pdf", width = 15, height = 10)

fb_merged_filt$prediction.score.max <- predictions$prediction.score.max

# Check prediction scores
VlnPlot(fb_merged_filt, features = "prediction.score.max",
        group.by = "predicted.celltype", pt.size = 0) +
  ggtitle("Prediction Confidence by Cell Type")
ggsave("pred_conf_cell_type.pdf", width = 15, height = 5)


p1 <- DimPlot_scCustom(fb_merged_filt, reduction = "umap", group.by = "predicted.celltype",
        label = TRUE, repel = TRUE, label.box = TRUE) +
  ggtitle("Transferred Cell Types - After Harmony")
p2 <- DimPlot_scCustom(fb_merged_filt, reduction = "umap", group.by = "seurat_clusters",
        label = TRUE, repel = TRUE, label.box = TRUE)
ggsave("umaps_post_harmony_and_xfer.pdf", width = 25, height = 10)

# =============================================================================
# 8. INITIAL CLUSTER ANNOTATION (v1–v2)
# =============================================================================

Idents(fb_merged_filt) <- "sample"
fb_seurat <- RenameIdents(fb_merged_filt, c("Ctrl13" = "158.4", "Ctrl24" = "147.4", "Ctrl25" = "144.3", "Ctrl26" = "157.4", 'Ctrl27' = "148.4", "Sal7" = "135.4", "Sal9" = "128.4", "Sal10" = "134.4", "Sal13" = "132.3"))
fb_seurat$gest_age <- Idents(fb_seurat)

saveRDS(fb_seurat, file = "fb_seurat_v1.RDS")

astro_plot <- DotPlot_scCustom(fb_seurat, features = c("GFAP", "NFIA", "SLC1A2", "SLC1A3", "SPARCL1",
                                                          "TIMP3", "APOE", "S100B", "VIM", "MT2A"),
                               group.by = "seurat_clusters", flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

astro_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
ggsave("astro_plot.pdf", plot = astro_plot, width = 10, height = 5)

olig_plot <- DotPlot_scCustom(fb_seurat, features = c("SOX10", "PDGFRA", "OLIG2", "PCDH17", "ENSMMUG00000056728", "PCDH15", "MMP16",
                                                      "CA10", "BCAS1", "ENPP6", "MAL", "MOG", "PLP1", "MBP"), group.by = "seurat_clusters",
                              flip_axes = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

olig_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
ggsave("oligo_plot.pdf", plot = olig_plot, width = 10, height = 5)

ipc_plot <- DotPlot_scCustom(fb_seurat,
                             features = c("DCX", "NEUROD1", "EOMES", "TBR1", "SOX5", "POU3F2", "NR2F1", 'NRP1', "CRYM", "SLC17A6",
                                          "NEUROD2", "TLE4", "SATB2", "UNC5D", "GRIA2", "CUX1", "CUX2", "RORB", "GRIN2B", "SLC17A7",
                                          "CAMK2A", "NRGN"
                                          ),
                             group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

ipc_plot_filtered <- ipc_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
ipc_plot_filtered
ggsave("ipc_plot.pdf", plot = ipc_plot_filtered, width = 10, height = 5)

astro_plot <- DotPlot_scCustom(fb_seurat, features = c("MT2A", "VIM", "S100B", "APOE", "TIMP3", "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP"),
                               group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

astro_plot_filtered <- astro_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
astro_plot_filtered
ggsave("astro_plot_2.pdf", plot = astro_plot_filtered, width = 10, height = 5)

miscel_plot <- DotPlot_scCustom(fb_seurat,
                                features = c("MKI67", "PECAM1", "VWF", "CD34", "CD248", "MYOF", "ABCC9", "GJA4", "HBA", "HBM",
                                             "HBG1", "HBG2", "HBF", "ALAS2", "ENSMMUG00000044429", "ENSMMUG00000041831", "HBD", "HBG1", "HBG2", "FECH",
                                             "BLVRB", "SLC25A37", "GATA1", "KLF1", "ANK1", "HBA2",  "CD3G", "CD3D", "CD69", "ICOS",
                                             "CD4", "CD8A", "FOXJ1", "SOX2", "CD133", "S100B", "TPPP3", "CCDC153", "GFAP", "CD68",
                                             "CD163", "AIF1", "C1QB", "C1QC",  "SALL1", "CX3CR1", "P2RY12", "MRC1", "CSF1R",  "GPR34",
                                             "CD80", "CD86", "IL1B"
                                             ),
                                group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

miscel_plot_filtered <- miscel_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
miscel_plot_filtered
ggsave("misc_plot.pdf", plot = miscel_plot_filtered, width = 10, height = 10)

immune_plot <- DotPlot_scCustom(fb_seurat,
                                features = c("PECAM1", "VWF", "HBA2", "CD3G", "CD3D", "CD69", "ICOS", "CD4", "CD8A", "CD163",
                                             "CD68", "AIF1", "TPI1", "P2RY12", "IRF8", "IL1B", "CD80", "CD86", "IL1", "ITGAM",
                                             "CD40", "IL6", "MAMUDRB1", "LDA", "C1QB", "C1QC", "SALL1", "CX3CR1", "MRC1", "CSF1R",
                                             "TMEM119", "LYVE1", "GPR34", "TREM2", "APOE", "CTSD"),
                                group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

immune_plot_filtered <- immune_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
immune_plot_filtered
ggsave("imm_plot.pdf", plot = immune_plot_filtered, width = 10, height = 10)

micro_plot <- DotPlot_scCustom(fb_seurat,
                               features = c("CD68", "CD163", "AIF1", "C1QB", "C1QC", "SALL1", "CX3CR1", "P2RY12", "MRC1", "CSF1R",
                                            "GPR34", "CD80", "CD86", "IL1B", "CTSD"),
                               group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

micro_plot_filtered <- micro_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
micro_plot_filtered
ggsave("micro_plot.pdf", plot = micro_plot_filtered, width = 10, height = 5)

cge_mge_plot <- DotPlot_scCustom(fb_seurat,
                                 features = c("MKI67", "CDK1", "TOP2A", "PCNA", "ASPM", "ANLN", "HOPX", "FAM107A", "TNC", "HES5",
                                              "NR2F1", "NR2F2", "DLX1", "DLX5", "GAD1", "GAD2", "TUBB3", "SP8", "NRP2", "PROX1",
                                              "CALB2", "KCNC1", "CHRNA2", "VIP", "NKX2-1", "SOX2", "DLX2", "SOX6", "SP9", "NPY",
                                              "LHX6",  "ERBB4", "GRIK3", "KCNC2", "SST", "PVALB"),
                                 group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) +
  labs(x = "Gene", y = "Cluster")

cge_mge_plot_filtered <- cge_mge_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
cge_mge_plot_filtered
ggsave("cge_mge_plot.pdf", plot = cge_mge_plot_filtered, width = 10, height = 10)

cluster_labels <- c(
  "7"  = "OL1",
  "5"  = "COP1",
  "0"  = "OPC",
  "13" = "vRG",
  "6"  = "AS0",
  "11" = "AS1",
  "3"  = "AS2",
  "14" = "oRG-G",
  "2"  = "CGE0",
  "18" = "PC",
  "22" = "FIB",
  "9"  = "EC",
  "10" = "Ep",
  "16" = "MG0",
  "17" = "IMM",
  "23" = "BAM"
)

Idents(fb_seurat) <- "seurat_clusters"
fb_seurat <- RenameIdents(fb_seurat, cluster_labels)
fb_seurat$cell_type_v2 <- Idents(fb_seurat)
saveRDS(fb_seurat, file = "fb_seurat_v2.RDS")


# =============================================================================
# 9. CLUSTER REFINEMENT & ANNOTATION (v3–v8)
# =============================================================================

# ---- Subcluster cluster 8 ----

Idents(fb_seurat) <- "seurat_clusters"
sub8 <- subset(fb_seurat, idents = "8")

sub8 <- SCTransform(sub8, vars.to.regress = "percent.mt")
sub8 <- RunPCA(sub8, npcs = 30)


sub8 <- RunHarmony(
  sub8,
  group.by.vars = "sample",
  reduction.use = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)

sub8 <- FindNeighbors(sub8, reduction = "harmony", dims = 1:20)
sub8 <- FindClusters(sub8, resolution = 0.3)
sub8 <- RunUMAP(sub8, reduction = "harmony", dims = 1:20)

p1 <- DimPlot_scCustom(sub8, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE) +
  ggtitle("Cluster 8 subclusters — Resolution 0.3")
save_umap(p1, "sub8_umap.pdf")

# ---- Subcluster cluster 12 ----

Idents(fb_seurat) <- "seurat_clusters"
sub12 <- subset(fb_seurat, idents = "12")

sub12 <- SCTransform(sub12, vars.to.regress = "percent.mt")
sub12 <- RunPCA(sub12, npcs = 30)


sub12 <- RunHarmony(
  sub12,
  group.by.vars = "sample",
  reduction.use = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)

sub12 <- FindNeighbors(sub12, reduction = "harmony", dims = 1:20)
sub12 <- FindClusters(sub12, resolution = 0.3)
sub12 <- RunUMAP(sub12, reduction = "harmony", dims = 1:20)

p1 <- DimPlot_scCustom(sub12, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE) +
  ggtitle("Cluster 12 subclusters — Resolution 0.3")
save_umap(p1, "sub12_umap.pdf")

# ---- Subcluster cluster 15 ----

Idents(fb_seurat) <- "seurat_clusters"
sub15 <- subset(fb_seurat, idents = "15")

sub15 <- SCTransform(sub15, vars.to.regress = "percent.mt")
sub15 <- RunPCA(sub15, npcs = 30)


sub15 <- RunHarmony(
  sub15,
  group.by.vars = "sample",
  reduction.use = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)

sub15 <- FindNeighbors(sub15, reduction = "harmony", dims = 1:20)
sub15 <- FindClusters(sub15, resolution = 0.3)
sub15 <- RunUMAP(sub15, reduction = "harmony", dims = 1:20)

p1 <- DimPlot_scCustom(sub15, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE) +
  ggtitle("Cluster 15 subclusters — Resolution 0.3")
save_umap(p1, "sub15_umap.pdf")


sub8 <- PrepSCTFindMarkers(sub8)

all_markers <- FindAllMarkers(object = sub8) %>%
  Add_Pct_Diff() 

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE, rank_by = "avg_log2FC")
plots <- Clustered_DotPlot(seurat_object = sub8, features = top_markers, flip = TRUE, x_lab_rotate = 90, elbow_kmax = 30)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = sub8, features = top_markers, flip = TRUE, x_lab_rotate = 90, k = 11)

top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20, data_frame = TRUE,
                                      rank_by = "avg_log2FC")


write.csv(top_20_markers, file = "top_20_markers_immune.csv")

pdf("sub8_mark.pdf", width = 20, height = 8)
dot_plot <- Clustered_DotPlot(seurat_object = sub8, features = top_markers, flip = TRUE, x_lab_rotate = 90, k = 11)
print(dot_plot)
dev.off()

Idents(fb_seurat) <- 'cell_type_v2'

fb_seurat <- FindSubCluster(
  fb_seurat,
  cluster = "8",
  graph.name = "SCT_snn",
  resolution = 0.05
)


Idents(fb_seurat) <- 'sub.cluster'

markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "8_1",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)

top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  pull(gene)

markers <- rownames_to_column(markers, var = "gene")
write_csv(markers, file = "unknown_0_markers.csv")


fb_seurat <- RenameIdents(fb_seurat, "8_1" = "ImmN", "8_0" = "UP/I")
fb_seurat$cell_type_v3 <- Idents(fb_seurat)


# Get barcodes for just cluster 1 from sub15
sub15_cluster1_cells <- colnames(sub15)[Idents(sub15) == "1"]

# Start from the latest label version
fb_seurat$cell_type_v4 <- as.character(fb_seurat$cell_type_v3)

# Overwrite only cluster 1 cells from sub15 with a new label
# Replace "Sub15_1" with whatever identity you assign it
fb_seurat$cell_type_v4[sub15_cluster1_cells] <- "Sub15_1"

fb_seurat$cell_type_v4 <- as.factor(fb_seurat$cell_type_v4)

# Visualize on main UMAP
Idents(fb_seurat) <- "cell_type_v4"
DimPlot_scCustom(fb_seurat, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = TRUE, group.by = "cell_type_v4")


run_lasso(fb_seurat)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
fb_seurat$cell_type_v3 <- as.character(fb_seurat$cell_type_v3)
fb_seurat$cell_type_v3[gated_cells] <- "COP2"

# Convert back to factor if needed
fb_seurat$cell_type_v3 <- as.factor(fb_seurat$cell_type_v3)

Idents(fb_seurat) <- "cell_type_v3"

fb_seurat <- RenameIdents(fb_seurat, "15" = "OL2")
fb_seurat$cell_type_v4 <- Idents(fb_seurat)

markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "20",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)


top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_unk, file = "20_markers.csv")


Idents(fb_seurat) <- "cell_type_v4"
fb_seurat <- RenameIdents(fb_seurat, "20" = "oRG")
fb_seurat$cell_type_v4 <- Idents(fb_seurat)

fb_seurat <- FindSubCluster(
  fb_seurat,
  cluster = "12",
  graph.name = "SCT_snn",
  resolution = 0.1
)

plot1 <- DimPlot_scCustom(fb_seurat, reduction = "umap", label = TRUE, repel = TRUE, aspect_ratio = 0.75, group.by = "sub.cluster")

ggsave("sub_main.pdf", plot = plot1, width = 15, height = 10)


Idents(fb_seurat) <- "sub.cluster"
fb_seurat <- FindSubCluster(
  fb_seurat,
  cluster = "12_0",
  graph.name = "SCT_snn",
  resolution = 0.1
)

cluster_labels <- c(
  "12_0_0"  = "ExN",
  "12_1"  = "ExN",
  "12_0_1"  = "CGE1",
  "12_2" = "IPC"
)

Idents(fb_seurat) <- "sub.cluster"
fb_seurat <- RenameIdents(fb_seurat, cluster_labels)
fb_seurat$cell_type_v5 <- Idents(fb_seurat)

plot1 <- DimPlot_scCustom(fb_seurat, reduction = "umap", label = TRUE, repel = TRUE, aspect_ratio = 0.75, group.by = "cell_type_v5")

ggsave("sub_main.pdf", plot = plot1, width = 15, height = 10)

FeaturePlot_scCustom(fb_seurat, features = "GFAP", reduction = "umap")
ggsave("gfap.pdf", width = 15, height = 10)

plot1 <- DimPlot_scCustom(fb_seurat, reduction = "umap", label = TRUE, repel = TRUE, aspect_ratio = 0.75, group.by = "sample")

ggsave("sample.pdf", plot = plot1, width = 15, height = 10)

FeaturePlot_scCustom(fb_seurat, features = "AIF1", reduction = "umap")
ggsave("AIF1.pdf", width = 15, height = 10)

saveRDS(fb_seurat, file = "fb_seurat_v3.RDS")

run_lasso(fb_seurat)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
fb_seurat$cell_type_v5 <- as.character(fb_seurat$cell_type_v5)
fb_seurat$cell_type_v5[gated_cells] <- "Ctrl27_MG"

# Convert back to factor if needed
fb_seurat$cell_type_v5 <- as.factor(fb_seurat$cell_type_v5)

Idents(fb_seurat) <- "cell_type_v5"
markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "Ctrl27_MG",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)


top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_unk, file = "Ctrl27_MG_markers.csv")


# Subset MG and do a feature plot fig of the following genes:
# Page 1: CD14, CD163, MRC1, AIF1, P2RY12, TMEM119, SALL1, CX3CR1, HEXB, ARG1, RETNLB, YM1, 
# Page 2: CD68, CD11b, IL1B, ISG15, IFIT2, CXCL10, NOS2, TREM2, APOE, CLEC7A, SPP1, HIF1A

clusters_of_interest <- c("1", "4", "21", "Ctrl27_MG", "MG0", "BAM")

mg_subset <- subset(fb_seurat, idents = clusters_of_interest)

run_lasso(mg_subset)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
mg_subset$cell_type_v5 <- as.character(mg_subset$cell_type_v5)
mg_subset$cell_type_v5[gated_cells] <- "MG_sub"

# Convert back to factor if needed
mg_subset$cell_type_v5 <- as.factor(mg_subset$cell_type_v5)

mg_subset$cell_type_v6 <- Idents(mg_subset)

Idents(mg_subset) <- "cell_type_v5"
mg_subset <- subset(mg_subset, idents = "MG_sub")


markers1 <- c('CD14', 'CD163', 'MRC1', 'AIF1', 'P2RY12', 'TMEM119', 
              'SALL1', 'CX3CR1', 'HEXB', 'ARG1', 'MERTK', "C1QB")

markers2 <- c('CD68', "CD80", 'ITGAM', 'IL1B', 'ISG15', 'IFIT2', 'CXCL10', 
              'NOS2', 'TREM2', 'APOE', 'CLEC7A', 'SPP1')

markers3 <- c('MKI67', 'DAB2', 'MCM5', 'AXL', 'CRYBB1', 
              'CSF1', 'CXCR2', 'MAFB', 'MEF2A', "TLR2", "IGF1", "SERPINE1")

markers4 <- c('CD163', 'MRC1', 'CSF1R', 'LYVE1', 'FOLR2', 'P2RY12', 'SALL1', "CX3CR1", "CCR2")

FeaturePlot_scCustom(mg_subset, features = markers1, reduction = "umap", num_columns = 3)
ggsave("mg_fig_1.PDF", width = 10, height = 12, units = "in", dpi = 300)
FeaturePlot_scCustom(mg_subset, features = markers2, reduction = "umap", num_columns = 3)
ggsave("mg_fig_2.PDF", width = 10, height = 12, units = "in", dpi = 300)
FeaturePlot_scCustom(mg_subset, features = markers3, reduction = "umap", num_columns = 3)
ggsave("mg_fig_3.PDF", width = 10, height = 12, units = "in", dpi = 300)
FeaturePlot_scCustom(mg_subset, features = markers4, reduction = "umap", num_columns = 3)
ggsave("mg_fig_4.PDF", width = 10, height = 10, units = "in", dpi = 300)

fb_seurat <- RenameIdents(fb_seurat, "1" = "MG1", "4" = "MG1", "21" = "MG1", "Ctrl27_MG" = "MG1")

fb_seurat$cell_type_v6 <- Idents(fb_seurat)

fb_seurat <- FindSubCluster(
  fb_seurat,
  cluster = "CGE0",
  graph.name = "SCT_snn",
  resolution = 0.3
)


mge_markers <- list(MGE = c("LHX6", "NKX2-1", "SOX6", "MAF", "MAFB"))
cge_markers <- list(CGE = c("NR2F2", "PROX1", "SP8", "NR2F1", "HTR3A"))

fb_seurat <- AddModuleScore(fb_seurat, features = mge_markers, name = "MGE_score")
fb_seurat <- AddModuleScore(fb_seurat, features = cge_markers, name = "CGE_score")


saveRDS(fb_seurat, file = "fb_seurat_v3.RDS")

cge1_cells <- WhichCells(fb_seurat, idents = "CGE1")


fib_dot <- DotPlot_scCustom(fb_seurat, 
                            features = c("COL1A2", "VIM", "DCN", "FMOD", "PDGFRB", 
                                         "MYOF", "TAGLN", "ABCC9", "PECAM1", "VWF", 
                                         "CD34", "FLT1", "SOX2", "CD133", "CCDC153", 
                                         "TTR", "CLDN3", "CLDN4", "FOXJ1", "RFX2", 
                                         "RFX3", "DNALI1", "TUBA4B", "CETN2", "AQP4"), 
                            group.by = "cell_type_v2", remove_axis_titles = FALSE,  
                            flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

fib_dot_filtered <- fib_dot +
  scale_y_discrete(limits = c("FIB", "PC", "EC", "19", "Ep"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
fib_dot_filtered

ggsave("Misc.pdf", width = 10, height = 8)


run_lasso(fb_seurat)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
fb_seurat$cell_type_v6 <- as.character(fb_seurat$cell_type_v6)
fb_seurat$cell_type_v6[gated_cells] <- "ChP"

# Convert back to factor if needed
fb_seurat$cell_type_v6 <- as.factor(fb_seurat$cell_type_v6)

# Code to relabel to nearest neighbors
nn_graph <- fb_seurat@graphs$SCT_snn

# Get current cell type labels — adjust column name to match yours
cell_labels <- fb_seurat@meta.data$cell_type_v6
names(cell_labels) <- colnames(fb_seurat)

cells_19 <- which(cell_labels == "19")
new_labels <- cell_labels

for (cell_idx in cells_19) {
  
  # Get neighbors
  neighbors <- which(nn_graph[cell_idx, ] > 0)
  neighbor_labels <- cell_labels[neighbors]
  
  # Exclude other 19 cells from the vote
  neighbor_labels <- neighbor_labels[neighbor_labels != "19"]
  
  # Majority vote
  if (length(neighbor_labels) > 0) {
    new_labels[cell_idx] <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]
  }
}

# Store new labels in fb_seurat
fb_seurat@meta.data$cell_type_relabeled <- new_labels

# Drop unused factor levels
fb_seurat$cell_type_relabeled <- droplevels(fb_seurat$cell_type_relabeled)

# Second KNN pass — reassign any remaining "19" cells
nn_graph <- fb_seurat@graphs$SCT_snn

cell_labels <- fb_seurat@meta.data$cell_type_relabeled
names(cell_labels) <- colnames(fb_seurat)

cells_19 <- which(cell_labels == "19")
new_labels <- cell_labels

for (cell_idx in cells_19) {
  
  # Get neighbors
  neighbors <- which(nn_graph[cell_idx, ] > 0)
  neighbor_labels <- cell_labels[neighbors]
  
  # Exclude other 19 cells from the vote
  neighbor_labels <- neighbor_labels[neighbor_labels != "19"]
  
  # Majority vote
  if (length(neighbor_labels) > 0) {
    new_labels[cell_idx] <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]
  }
}

# Store new labels in fb_seurat
fb_seurat@meta.data$cell_type_relabeled_2 <- new_labels

# Drop unused factor levels
fb_seurat$cell_type_relabeled_2 <- droplevels(fb_seurat$cell_type_relabeled_2)

saveRDS(fb_seurat, file = "fb_seurat_v4.RDS")

run_lasso(fb_seurat)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
fb_seurat$cell_type_relabeled_2 <- as.character(fb_seurat$cell_type_relabeled_2)
fb_seurat$cell_type_relabeled_2[gated_cells] <- "EC"

# Convert back to factor if needed
fb_seurat$cell_type_relabeled_2 <- as.factor(fb_seurat$cell_type_relabeled_2)

saveRDS(fb_seurat, file = "fb_seurat_v5.RDS")

Idents(fb_seurat) <- "cell_type_relabeled_2"


# =============================================================================
# 10. LINEAGE SUBSET ANALYSES
# =============================================================================

# Immune subset
clusters_of_interest <- c("IMM", "MG0", "MG1", "BAM")

mg_subset <- subset(fb_seurat, idents = clusters_of_interest)

p1 <- DimPlot_scCustom(mg_subset, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, group.by = "cell_type_relabeled_2") +
  ggtitle("Immune subcluster")
save_umap(p1, "mg_subset_umap_orig.pdf")

saveRDS(mg_subset, file = "mg_subset_original.RDS")

mg_subset <- SCTransform(mg_subset, vars.to.regress = "percent.mt")
mg_subset <- RunPCA(mg_subset, npcs = 30)


mg_subset <- RunHarmony(
  mg_subset,
  group.by.vars = "sample",
  reduction.use = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)

mg_subset <- FindNeighbors(mg_subset, reduction = "harmony", dims = 1:20)
mg_subset <- FindClusters(mg_subset, resolution = 0.3)
mg_subset <- RunUMAP(mg_subset, reduction = "harmony", dims = 1:20)

p1 <- DimPlot_scCustom(mg_subset, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, group.by = "cell_type_relabeled_2") +
  ggtitle("Immune subcluster")
save_umap(p1, "mg_subset_umap.pdf")

saveRDS(mg_subset, file = "mg_subset_recluster.RDS")

#astro subset
clusters_of_interest <- c("UP/I", "vRG", "oRG", "oRG-G", "AS0", "AS1", "AS2")

astro_subset <- subset(fb_seurat, idents = clusters_of_interest)

saveRDS(astro_subset, file = "astro_subset_original.RDS")

astro_subset <- SCTransform(astro_subset, vars.to.regress = "percent.mt")
astro_subset <- RunPCA(astro_subset, npcs = 30)


astro_subset <- RunHarmony(
  astro_subset,
  group.by.vars = "sample",
  reduction.use = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)

astro_subset <- FindNeighbors(astro_subset, reduction = "harmony", dims = 1:20)
astro_subset <- FindClusters(astro_subset, resolution = 0.3)
astro_subset <- RunUMAP(astro_subset, reduction = "harmony", dims = 1:20)

p1 <- DimPlot_scCustom(astro_subset, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, group.by = "cell_type_relabeled_2") +
  ggtitle("Astro subcluster")
save_umap(p1, "astro_subset_umap.pdf")

saveRDS(astro_subset, file = "astro_subset_recluster.RDS")


Idents(fb_seurat) <- "cell_type_relabeled_2"
markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "oRG",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)

top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_unk, file = "oRG_markers.csv")

Idents(fb_seurat) <- "cell_type_relabeled_2"
markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "ChP",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)

top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_unk, file = "ChP_markers.csv")

Idents(fb_seurat) <- "cell_type_relabeled_2"
fb_seurat <- RenameIdents(fb_seurat, "oRG" = "APC")
fb_seurat$cell_type_v7 <- Idents(fb_seurat)

Idents(fb_seurat) <- "cell_type_v7"
markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "oRG-G",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)

top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_unk, file = "oRG-G_markers.csv")

Idents(fb_seurat) <- "cell_type_v7"
fb_seurat <- RenameIdents(fb_seurat, "AS2" = "AS3")
fb_seurat <- RenameIdents(fb_seurat, "AS1" = "AS2")
fb_seurat <- RenameIdents(fb_seurat, "AS0" = "AS1")
fb_seurat <- RenameIdents(fb_seurat, "oRG-G" = "AS0")
fb_seurat$cell_type_v7 <- Idents(fb_seurat)

run_lasso(fb_seurat)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
fb_seurat$cell_type_v7 <- as.character(fb_seurat$cell_type_v7)
fb_seurat$cell_type_v7[gated_cells] <- "PC"

# Convert back to factor if needed
fb_seurat$cell_type_v7 <- as.factor(fb_seurat$cell_type_v7)

Idents(fb_seurat) <- "cell_type_v7"


fb_seurat <- RenameIdents(fb_seurat, "AS1" = "AST")
fb_seurat <- RenameIdents(fb_seurat, "AS2" = "AS1")
fb_seurat <- RenameIdents(fb_seurat, "AST" = "AS2")
fb_seurat$cell_type_v7 <- Idents(fb_seurat)

Idents(fb_seurat) <- "cell_type_v7"
fb_seurat <- RenameIdents(fb_seurat, "AS2" = "oRG-G", "AS3" = "AS2")
fb_seurat$cell_type_v7 <- as.factor(as.character(Idents(fb_seurat)))

saveRDS(fb_seurat, file = "fb_seurat_v6.RDS")


# astro subset
astro_sub <- subset(fb_seurat,
                    subset = cell_type_v7 %in% c("APC", "AS0", "AS1", "AS2", "oRG-G"))

astro_sub <- astro_sub %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sample", reduction = "pca",
               reduction.save = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
astro_sub <- FindNeighbors(astro_sub, reduction = "harmony", dims = 1:30)
astro_sub <- FindClusters(astro_sub, resolution = 0.3)


saveRDS(astro_sub, file = "astro_sub.RDS")

vrg_dot <- DotPlot_scCustom(fb_seurat, 
                            features = c("MKI67", 
                                         "TOP2A", "PCNA", "ASPM", "CENPF", 
                                         "CDK1", "ANLN", "VIM", "S100A6", "HIF1A","MT2A", "MT3A",
                                         "HOPX", "FAM107A", "TNC",
                                         "S100B", "APOE", "TIMP3", 
                                         "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP",
                                         "AQP4", "GJA1", "AGT" ,"MLC1", "GLAST", "ATP1A2" ,"ATP1B2", "MFGE8", "GLUL" ,"ALDOC", "CKB",
                                         "NDRG2", "HEPACAM", "HES1" ,"HEY1", "SOX9", "ID3" ,"ALDH1A1",
                                          "FABP7", "PAX6", "LHX2", "MT3", "SFRP1",
                                         "CXCL14", "SPARCL", "CLU",  "PTN", "CST3", "CPE", "EDNRB", "BBOX1", "GDPD2", "ACSBG1", "HSD17B6"
                                         ), 
                            group.by = "cell_type_v7", remove_axis_titles = FALSE,  
                            flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

vrg_dot_filtered <- vrg_dot +
  scale_y_discrete(limits = c("APC", "AS0", "AS1", "AS2", "AS3"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
vrg_dot_filtered

ggsave("astro_dot.pdf", vrg_dot_filtered, width = 8, height = 12)

# =============================================================================
# 11. scVELO METADATA EXPORT
# =============================================================================

markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "AS2",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
)

top_markers_as2 <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_as2, file = "AS2_markers.csv")

save_subset_for_scvelo(fb_seurat, "full")

# oligo subset
oligo_sub <- subset(fb_seurat,
                    subset = cell_type_v7 %in% c("OPC", "OL1", "OL2", "COP1", "COP2", "vRG"))

oligo_sub <- oligo_sub %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sample", reduction = "pca",
               reduction.save = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
oligo_sub <- FindNeighbors(oligo_sub, reduction = "harmony", dims = 1:30)
oligo_sub <- FindClusters(oligo_sub, resolution = 0.3)


saveRDS(oligo_sub, file = "oligo_sub.RDS")

# MG subset
mg_sub <- subset(fb_seurat,
                    subset = cell_type_v7 %in% c("MG0", "MG1"))

mg_sub <- mg_sub %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sample", reduction = "pca",
               reduction.save = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
mg_sub <- FindNeighbors(mg_sub, reduction = "harmony", dims = 1:30)
mg_sub <- FindClusters(mg_sub, resolution = 0.3)

plot1 <- FeaturePlot_scCustom(mg_sub, features = c('IL1B', 'P2RY12', 'CX3CR1', 'SALL1', 'APOE', 'ISG15', 'LGALS3', 'AIF1', 'SPP1'), reduction = "umap")
ggsave("mg_feat.pdf", plot = plot1, width = 15, height = 10)

Idents(mg_sub) <- "cell_type_v7"
markers <- FindMarkers(
  object = mg_sub,
  ident.1 = "MG1",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
) %>%
  Add_Pct_Diff()

top_markers_MG1 <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  percent

plot4 <- FeaturePlot_scCustom(mg_sub, features = c('CD86'), reduction = "umap")
ggsave("mg_feat4.pdf", plot = plot4, width = 15, height = 10)
write_csv(top_markers_MG1, file = "MG1_markers.csv")

saveRDS(mg_sub, file = "mg_sub.RDS")

# BAM MG subset
bam_mg_sub <- subset(fb_seurat,
                    subset = cell_type_v7 %in% c("BAM", "MG0", "MG1"))

# Use harmony embeddings that already exist on these cells
bam_mg_sub <- RunUMAP(bam_mg_sub, reduction = "harmony", dims = 1:30,
                     reduction.name = "umap")
bam_mg_sub <- FindNeighbors(bam_mg_sub, reduction = "harmony", dims = 1:30)
bam_mg_sub <- FindClusters(bam_mg_sub, resolution = 0.3)


bam_mg_sub <- bam_mg_sub %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sample", reduction = "pca",
               reduction.save = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
bam_mg_sub <- FindNeighbors(bam_mg_sub, reduction = "harmony", dims = 1:30)
bam_mg_sub <- FindClusters(bam_mg_sub, resolution = 0.3)


saveRDS(bam_mg_sub, file = "bam_mg_sub.RDS")

# Neuro subset
neuro_sub <- subset(fb_seurat,
                    subset = cell_type_v7 %in% c("ImmN", "UP/I", "IPC", "ExN", "CGE0", "CGE1"))

neuro_sub <- neuro_sub %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sample", reduction = "pca",
               reduction.save = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
neuro_sub <- FindNeighbors(neuro_sub, reduction = "harmony", dims = 1:30)
neuro_sub <- FindClusters(neuro_sub, resolution = 0.3)


saveRDS(neuro_sub, file = "neuro_sub.RDS")


# ExN subset
exn_sub <- subset(fb_seurat,
                    subset = cell_type_v7 %in% c("ImmN", "IPC", "ExN"))

exn_sub <- exn_sub %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sample", reduction = "pca",
               reduction.save = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
exn_sub <- FindNeighbors(exn_sub, reduction = "harmony", dims = 1:30)
exn_sub <- FindClusters(exn_sub, resolution = 0.3)


saveRDS(exn_sub, file = "exn_sub.RDS")

save_subset_for_scvelo(astro_sub, "astro")
save_subset_for_scvelo(neuro_sub, "neuro")
save_subset_for_scvelo(exn_sub, "exn")
save_subset_for_scvelo(mg_sub, "mg")
save_subset_for_scvelo(oligo_sub, "oligo")

exn_dot <- DotPlot_scCustom(fb_seurat, 
                            features = c("DCX", "NEUROD1", "EOMES", "TBR1", "SOX5", 
                                         "POU3F2", "NR2F1", "SLC17A6", "NEUROD2", "SATB2", 
                                         "UNC5D", "GRIA2", "CUX1", "CUX2", "GRIN2B", 
                                         "SLC17A7", "CAMK2A", "NRGN", "NKX2-1", "LHX6", 
                                         "SOX6", "DLX1", "DLX5", "DLX6", "GAD1", 
                                         "GAD2", "PROX1", "CALB2", "HTR3A", "NRSN1", 
                                         "POU2F2", "RELN" ,"VIP", "CCK", "SST", 
                                         "NPY", "PVALB"), 
                             group.by = "cell_type_v7", remove_axis_titles = FALSE,  
                            flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

exn_dot_filtered <- exn_dot +
  scale_y_discrete(limits = c("IPC", "ExN", "ImmN", "CGE0", "CGE1"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
exn_dot_filtered

run_lasso(neuro_sub)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
neuro_sub$cell_type_v7 <- as.character(neuro_sub$cell_type_v7)
neuro_sub$cell_type_v7[gated_cells] <- "UNK"

# Convert back to factor if needed
neuro_sub$cell_type_v7 <- as.factor(neuro_sub$cell_type_v7)

Idents(neuro_sub) <- "cell_type_v7"

cux2_plot <- FeaturePlot_scCustom(neuro_sub, features = "CUX2", reduction = "umap")
ggsave("cux2_plot.pdf", plot = cux2_plot)
neuro_plot <- DimPlot(neuro_sub, group.by = "cell_type_v7", reduction = "umap")
ggsave("neuro_plot.pdf", plot = neuro_plot)

run_lasso(neuro_sub)
gated_cells <- readRDS("gated_cells.rds")

# Add gated cells as a new cell type in your existing cluster column
neuro_sub$cell_type_v7 <- as.character(neuro_sub$cell_type_v7)
neuro_sub$cell_type_v7[gated_cells] <- "CGE"

# Convert back to factor if needed
neuro_sub$cell_type_v7 <- as.factor(neuro_sub$cell_type_v7)

Idents(neuro_sub) <- "cell_type_v7"

neuro_sub <- FindSubCluster(
  neuro_sub,
  cluster = "CGE",
  graph.name = "SCT_snn",
  resolution = 0.05
)

Idents(neuro_sub) <- neuro_sub$cell_type_relabeled
# Transfer new labels back to the full object

# Convert to character first
fb_seurat$cell_type_v8 <- as.character(fb_seurat$cell_type_v7)

# Now the assignment will work
fb_seurat$cell_type_v8[colnames(neuro_sub)] <- as.character(Idents(neuro_sub))

fb_seurat$cell_type_v8 <- as.factor(fb_seurat$cell_type_v8)


# Visualize on the main UMAP to confirm
Idents(fb_seurat) <- "cell_type_v8"

Idents(neuro_sub) <- "cell_type_v7"

# KNN relabeling to reassign CGE1 cells by majority vote of neighbors
nn_graph <- neuro_sub@graphs$SCT_snn
cell_labels <- neuro_sub@meta.data$cell_type_relabeled
names(cell_labels) <- colnames(neuro_sub)

cells_CGE0 <- which(cell_labels == "CGE1")
new_labels <- cell_labels

for (cell_idx in cells_CGE0) {
  
  # Get neighbors
  neighbors <- which(nn_graph[cell_idx, ] > 0)
  neighbor_labels <- cell_labels[neighbors]
  
  # Exclude other cluster cells from the vote
  neighbor_labels <- neighbor_labels[neighbor_labels != "CGE1"]
  
  # Majority vote
  if (length(neighbor_labels) > 0) {
    new_labels[cell_idx] <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]
  }
}

# Store new labels in fb_seurat
neuro_sub@meta.data$cell_type_relabeled <- new_labels

# Drop unused factor levels
neuro_sub$cell_type_relabeled <- droplevels(neuro_sub$cell_type_relabeled)

Idents(neuro_sub) <- "cell_type_relabeled"
neuro_sub <- RenameIdents(neuro_sub, "CGE_0" = "CGE", "CGE_1" = "MGE")
neuro_sub$cell_type_v7 <- Idents(neuro_sub)

Idents(neuro_sub) <- neuro_sub$cell_type_v7
# Transfer new labels back to the full object

# Convert to character first
fb_seurat$cell_type_v8 <- as.character(fb_seurat$cell_type_v7)

# Now the assignment will work
fb_seurat$cell_type_v8[colnames(neuro_sub)] <- as.character(Idents(neuro_sub))

fb_seurat$cell_type_v8 <- as.factor(fb_seurat$cell_type_v8)


# Visualize on the main UMAP to confirm
Idents(fb_seurat) <- "cell_type_v8"


saveRDS(neuro_sub, file = "neuro_sub_v2.RDS")
saveRDS(fb_seurat, file = "fb_seurat_FINAL_v2.RDS")


markers <- FindMarkers(
  object = fb_seurat,
  ident.1 = "UNK",   # Your cluster of interest
  min.pct = 0.25,             # Min fraction of cells expressing the gene
  logfc.threshold = 0.25      # Min log fold-change threshold
) %>%
Add_Pct_Diff()


top_markers_unk <- markers %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 


write_csv(top_markers_unk, file = "UNK_fb_markers.csv")

Idents(fb_seurat) <- "cell_type_v8"
fb_seurat <- RenameIdents(fb_seurat, "UNK" = "MGE-SST")

fb_seurat$cell_type_v8 <- Idents(fb_seurat)

# convert to char and back to factor to get in alphabetical order
fb_seurat$cell_type_v8 <- as.character(fb_seurat$cell_type_v8)
fb_seurat$cell_type_v8 <- as.factor(fb_seurat$cell_type_v8)


# =============================================================================
# 12. PAPER FIGURES
# =============================================================================


# Generate the base plot without labels
p <- DimPlot_scCustom(
  fb_seurat,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  repel = TRUE,
  colors_use = DiscretePalette_scCustomize(
    num_colors = 36,
    palette = "varibow",
    shuffle_pal = FALSE
  ),
  color_seed = 42
)

# Check column names if needed: colnames(p$data)

# Compute label positions
label_data <- p$data %>%
  group_by(ident) %>%
  summarise(
    x = median(.data[[grep("_1$", colnames(p$data), value = TRUE)[1]]]),
    y = median(.data[[grep("_2$", colnames(p$data), value = TRUE)[1]]])
  )

# Plot with repelled shadow labels
p2 <- p + geom_text_repel(
  data = label_data,
  aes(x = x, y = y, label = ident),
  colour = "black",
  bg.color = "white",
  bg.r = 0.15,
  size = 4,
  fontface = "bold",
  max.overlaps = Inf,
  box.padding = 0.5,
  point.padding = 0.3,
  segment.color = "grey50",
  segment.size = 0.3,
  min.segment.length = 0.2
)+
  xlab("UMAP 1") +
  ylab("UMAP 2")
ggsave(
  "varibow_false.pdf",
  plot = p2,
  width = 15,
  height = 10,
  device = "pdf"
)

plot_umap_genes(neuro_sub, genes = c('EOMES', 'GAD1', 'SATB2'), outfile = "neuro_genes.pdf")

plot_umap_genes(neuro_sub, genes = c('LHX6', 'PROX1', 'TBR1'), outfile = "neuro_genes_2.pdf")

plot_umap_genes(oligo_sub, genes = c('PDGFRA', 'BCAS1', 'MAL'), outfile = "oligo_genes.pdf")

plot_umap_genes(astro_sub, genes = c('EGFR', 'GFAP', 'HOPX'), outfile = "astro_genes.pdf")

plot_umap_genes(mg_sub, genes = c('IL1B', 'SALL1'), outfile = "mg_genes.pdf")

colors <- get_atlas_colors(fb_seurat)
plot_umap_clusters(fb_seurat, custom_colors = colors, outfile = "umap.pdf")


Idents(neuro_sub) <- "cell_type_v8"
colors <- get_subset_colors(neuro_sub)
plot_umap_clusters(neuro_sub, custom_colors = colors, outfile = "neuro_umap.pdf")

neuro_sub$cell_type_v8 <- droplevels(neuro_sub$cell_type_v8)

colors <- get_subset_colors(neuro_sub, group_by = "cell_type_v8")

# Swap IPC and ImmN
colors[c("IPC", "ImmN")] <- colors[c("ImmN", "IPC")]

plot_umap_clusters(neuro_sub,
                   group_by      = "cell_type_v8",
                   custom_colors = colors,
                   outfile       = "neuro_umap.pdf")

plot_umap_clusters(neuro_sub,
                   custom_colors = get_subset_colors(neuro_sub),
                   outfile = "neuro_umap.pdf")

Idents(oligo_sub) <- "cell_type_v8"
colors <- get_subset_colors(oligo_sub)
plot_umap_clusters(oligo_sub, custom_colors = colors, outfile = "oligo_umap.pdf")

Idents(astro_sub) <- "cell_type_v8"
colors <- get_subset_colors(astro_sub)
plot_umap_clusters(astro_sub, custom_colors = colors, outfile = "astro_umap.pdf")

Idents(bam_mg_sub) <- "cell_type_v8"
colors <- get_subset_colors(bam_mg_sub)
plot_umap_clusters(bam_mg_sub, custom_colors = colors, outfile = "bam_mg_umap.pdf")

fb_seurat <- RenameIdents(fb_seurat, "MGE-SST" = "MGE-IN")

fb_seurat$cell_type_v8 <- Idents(fb_seurat)
# convert to char and back to factor to get in alphabetical order
fb_seurat$cell_type_v8 <- as.character(fb_seurat$cell_type_v8)
fb_seurat$cell_type_v8 <- as.factor(fb_seurat$cell_type_v8)
Idents(fb_seurat) <- "cell_type_v8"

saveRDS(fb_seurat, "fb_seurat_FINAL_v2.RDS")

Idents(neuro_sub) <- "cell_type_v8"
neuro_sub <- RenameIdents(neuro_sub, "MGE-SST" = "MGE-IN")
neuro_sub$cell_type_v8 <- Idents(neuro_sub)
# convert to char and back to factor to get in alphabetical order
neuro_sub$cell_type_v8 <- as.character(neuro_sub$cell_type_v8)
neuro_sub$cell_type_v8 <- as.factor(neuro_sub$cell_type_v8)
Idents(neuro_sub) <- "cell_type_v8"

saveRDS(neuro_sub, "neuro_sub_FINAL.RDS")

plot_umap_clusters(fb_seurat,
                     custom_colors = get_atlas_colors(fb_seurat),
                      outfile = "umap.pdf")

plot_harmony_integration(
  fb_seurat,
  outfile_combined = "harmony_combined.pdf",
  outfile_sample   = "harmony_sample.pdf",
  outfile_gest_age = "harmony_gest_age.pdf"
)

Idents(bam_mg_sub) <- "sample"
colors <- get_subset_colors(bam_mg_sub)
plot_umap_clusters(bam_mg_sub, custom_colors = colors, outfile = "bam_mg_sample_umap.pdf", group_by = "sample")

unique(fb_seurat$sample)

Idents(fb_seurat) <- "sample"
fb_seurat <- RenameIdents(fb_seurat, "Ctrl13" = "Media", "Ctrl24" = "Media", "Ctrl25" = "Media", "Ctrl26" = "Media", "Ctrl27" = "Media", "Sal7" = "Catheterized", "Sal9" = "Catheterized", "Sal10" = "Catheterized", "Sal13" = "Catheterized")

fb_seurat$treatment <- Idents(fb_seurat)

fb_seurat <- RenameIdents(fb_seurat, "Ctrl13" = "Ctrl5", "Ctrl24" = "Ctrl6", "Ctrl25" = "Ctrl7", "Ctrl26" = "Ctrl8", "Ctrl27" = "Ctrl9", "Sal7" = "Ctrl1", "Sal9" = "Ctrl2", "Sal10" = "Ctrl3", "Sal13" = "Ctrl4")

fb_seurat$id_project <- Idents(fb_seurat)

saveRDS(fb_seurat, 'fb_seurat_FINAL_v2.RDS')


# Step 3: Run scComp
sccomp_result <- fb_seurat@meta.data %>%
  rownames_to_column("cell_id") %>%
  sccomp_estimate(
    formula_composition = ~ treatment,
    .sample = sample,
    .cell_group = cell_type_v8,
    bimodal_mean_variability_association = TRUE,
    cores = 4
  )

# Step 4: Test Media vs Saline
sccomp_result <- sccomp_result %>%
  sccomp_test(contrasts = c("treatmentCatheterized"))

# Compute proportions
prop_data <- fb_seurat@meta.data %>%
  group_by(sample, cell_type_v8, treatment) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()


# Add significance marker to facet labels
sig_labels <- sccomp_result %>%
  select(cell_type_v8, c_FDR) %>%
  mutate(cell_type_label = ifelse(c_FDR < 0.05,
                                  paste0(cell_type_v8, " *"),
                                  cell_type_v8))

prop_data <- prop_data %>%
  left_join(sig_labels, by = "cell_type_v8")

plot_box <- ggplot(prop_data, aes(x = treatment, y = proportion, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ cell_type_label, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    strip.text  = element_text(size = 14, color = "black"),
    axis.title  = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  labs(x = "Treatment", y = "Proportion")

pdf("treat_sccomp_boxplot.pdf", width = 12, height = 10)
print(plot_box)
dev.off()

Idents(fb_seurat) <- "sample"
fb_seurat <- RenameIdents(fb_seurat, "Ctrl13" = "F", "Ctrl24" = "M", "Ctrl25" = "M", "Ctrl26" = "F", "Ctrl27" = "F", "Sal7" = "F", "Sal9" = "M", "Sal10" = "F", "Sal13" = "F")

fb_seurat$fet_sex <- Idents(fb_seurat)

saveRDS(fb_seurat, 'fb_seurat_FINAL_v2.RDS')

# Step 3: Run scComp
sccomp_result <- fb_seurat@meta.data %>%
  rownames_to_column("cell_id") %>%
  sccomp_estimate(
    formula_composition = ~ fet_sex,
    .sample = sample,
    .cell_group = cell_type_v8,
    bimodal_mean_variability_association = TRUE,
    cores = 4
  )

# Step 4: Test Media vs Saline
sccomp_result <- sccomp_result %>%
  sccomp_test(contrasts = c("fet_sexM"))




# Compute proportions
prop_data <- fb_seurat@meta.data %>%
  group_by(sample, cell_type_v8, fet_sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Add significance labels from fet_sex sccomp_result
sig_labels <- sccomp_result %>%
  select(cell_type_v8, c_FDR) %>%
  mutate(cell_type_label = ifelse(c_FDR < 0.05,
                                  paste0(cell_type_v8, " *"),
                                  cell_type_v8))

prop_data <- prop_data %>%
  left_join(sig_labels, by = "cell_type_v8")

plot_box <- ggplot(prop_data, aes(x = fet_sex, y = proportion, fill = fet_sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ cell_type_label, scales = "free_y", ncol = 5) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    strip.text  = element_text(size = 14, color = "black"),
    axis.title  = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  labs(x = "Fetal Sex", y = "Proportion")

pdf("fet_sex_sccomp_boxplot.pdf", width = 12, height = 10)
print(plot_box)
dev.off()

bam_mg_sub <- readRDS("bam_mg_sub.RDS")
Idents(bam_mg_sub) <- "sample"
bam_mg_sub <- RenameIdents(bam_mg_sub, "Ctrl13" = "CTRL5", "Ctrl24" = "CTRL6", "Ctrl25" = "CTRL7", "Ctrl26" = "CTRL8", "Ctrl27" = "CTRL9", "Sal7" = "CTRL1", "Sal9" = "CTRL2", "Sal10" = "CTRL3", "Sal13" = "CTRL4")

bam_mg_sub$id_project <- Idents(bam_mg_sub)
# convert to char and back to factor to get in alphabetical order
bam_mg_sub$id_project <- as.character(bam_mg_sub$id_project)
bam_mg_sub$id_project <- as.factor(bam_mg_sub$id_project)
Idents(bam_mg_sub) <- "id_project"
colors <- get_subset_colors(bam_mg_sub)
plot_umap_clusters(bam_mg_sub, custom_colors = colors, outfile = "bam_mg_sample_umap.pdf", group_by = "id_project")

mg_sub <- readRDS("mg_sub.RDS")
plot_umap_genes(mg_sub, genes = c('FGF14', 'TIMP1', 'CALD1'), outfile = "mg0_genes.pdf")
plot_umap_genes(mg_sub, genes = c('KCNQ3', 'SMAD3', 'PALD1'), outfile = "mg1_genes.pdf")

#Cell cycle scoring

# Load Seurat's built-in cell cycle gene lists (human genes, uppercase)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Run cell cycle scoring
fb_seurat <- CellCycleScoring(
  fb_seurat,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE  # TRUE would overwrite active identity with Phase
)

# Check results — adds S.Score, G2M.Score, and Phase columns to metadata
head(fb_seurat@meta.data[, c("S.Score", "G2M.Score", "Phase")])
table(fb_seurat$cell_type_v8, fb_seurat$Phase)

# 1. UMAP colored by cell cycle phase
p_phase <- DimPlot_scCustom(
  fb_seurat, reduction = "umap",
  group.by = "Phase",
  aspect_ratio = 0.75
) +
  ggtitle("Cell Cycle Phase") +
  theme(
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold")
  )

save_umap(p_phase, "cell_cycle_phase.pdf")

# 2. S and G2M scores as continuous feature plots (side by side)
p_s   <- FeaturePlot_scCustom(fb_seurat, features = "S.Score",   reduction = "umap")
p_g2m <- FeaturePlot_scCustom(fb_seurat, features = "G2M.Score", reduction = "umap")

save_umap(p_s | p_g2m, "cell_cycle_scores.pdf", width = 20)

# 3. Phase overlaid on your cell type UMAP (split view)
p_celltype <- DimPlot_scCustom(
  fb_seurat, reduction = "umap",
  group.by = "cell_type_v8",
  label = TRUE, repel = TRUE, aspect_ratio = 0.75
)
p_phase2 <- DimPlot_scCustom(
  fb_seurat, reduction = "umap",
  group.by = "Phase", aspect_ratio = 0.75
)

save_umap(p_celltype | p_phase2, "celltype_vs_phase.pdf", width = 24)

# 4. Phase composition per cell type (bar chart)
phase_df <- fb_seurat@meta.data %>%
  group_by(cell_type_v8, Phase) %>%
  tally() %>%
  group_by(cell_type_v8) %>%
  mutate(pct = n / sum(n))

ggplot(phase_df, aes(x = cell_type_v8, y = pct, fill = Phase)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion", title = "Cell Cycle Phase by Cell Type") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold")
  )

ggsave("cell_cycle_by_celltype.pdf", width = 10, height = 8)

# 5. Violin plot of scores per cell type
VlnPlot(fb_seurat, features = c("S.Score", "G2M.Score"),
        group.by = "cell_type_v8", pt.size = 0, ncol = 2)
ggsave("cell_cycle_vln.pdf", width = 20, height = 6)


fib_dot <- DotPlot_scCustom(fb_seurat, 
                            features = c('RPL37A', 'RPL32', 'RPS24', 'RPS2', 'RPL23A', "GREM1", "ENSMMUG00000012525", "HMGA2", "FN1", "COL1A2",  "VIM", "DCN", "FMOD", "PDGFRB", 
                                         "MYOF", "TAGLN", "ABCC9", "PECAM1", "VWF", 
                                         "CD34", "FLT1", "SOX2", "CD133", "CCDC153", 
                                         "TTR", "CLDN3", "CLDN4", 
                                         "FOLR1", "KL", "CLIC6", "AQP1", 
                                         "FOXJ1", "RFX2", 
                                         "RFX3", "DNALI1", "TUBA4B", "CETN2", "AQP4"
                                        
                                         ), 
                            group.by = "cell_type_v8", remove_axis_titles = FALSE,  
                            flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

fib_dot_filtered <- fib_dot +
  scale_y_discrete(limits = c("UP/I", "MPC", "FIB", "PC", "EC", "ChP", "Ep"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
fib_dot_filtered

ggsave("misc_dot.pdf", plot = fib_dot_filtered, width = 8, height = 10)

sessionInfo()

p1 <- FeaturePlot_scCustom(neuro_sub, features = "NR2F1", reduction = "umap")
ggsave("neuro_gene.pdf", plot = p1)

exn_dot <- DotPlot_scCustom(fb_seurat,
                            features = c("SLC17A7", "CAMK2A", "NRGN", "DCX", "NEUROD1", "EOMES", "TBR1", "SOX5", 
                                         "POU3F2", "NR2F1", "SLC17A6", "NEUROD2", "SATB2", 
                                         "UNC5D", "GRIA2", "CUX1", "CUX2", "GRIN2B","GAD1", 
                                         "GAD2", "DLX1", "DLX5", "DLX6", "MBIP", "ZSWIM5", "PROX1",
                                         "NKX2-1", "LHX6",
                                         "SOX6",  "CALB2", "ST18", "HTR3A", "NRSN1", 
                                         "POU2F2", "RELN" ,"VIP", "CCK", "SST", 
                                         "NPY", "PVALB"), 
                             group.by = "cell_type_v8", remove_axis_titles = FALSE,  
                            flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

exn_dot_filtered <- exn_dot +
  scale_y_discrete(limits = c( "ImmN","IPC", "ExN", "CGE", "MGE", "MGE-IN"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
exn_dot_filtered

ggsave("neuro_dot.pdf", plot = exn_dot_filtered, width = 8, height = 10)

p1 <- FeaturePlot_scCustom(fb_seurat, features = c('NES', 'GFAP', 'VIM', 'EGFR'), reduction = "umap")
ggsave("gene.pdf", plot = p1)

p1 <- FeaturePlot_scCustom(fb_seurat, features = "ENSMMUG00000056728", reduction = "umap")
ggsave("gene1.pdf", plot = p1)

Idents(fb_seurat) <- "cell_type_v8"
fb_seurat <- RenameIdents(fb_seurat, "vRG" = "pre-OPC")
fb_seurat$cell_type_v8 <- Idents(fb_seurat)
# convert to char and back to factor to get in alphabetical order
fb_seurat$cell_type_v8 <- as.character(fb_seurat$cell_type_v8)
fb_seurat$cell_type_v8 <- as.factor(fb_seurat$cell_type_v8)
Idents(fb_seurat) <- "cell_type_v8"

saveRDS(fb_seurat, file = "fb_seurat_FINAL_v2.RDS")

fb_seurat <- RenameIdents(fb_seurat, "APC" = "AST+", "AS0" = "AST0", "AS1" = "AST1", "AS2" = "AST2")
fb_seurat$cell_type_v8 <- Idents(fb_seurat)
fb_seurat$cell_type_v8 <- as.character(fb_seurat$cell_type_v8)
fb_seurat$cell_type_v8 <- as.factor(fb_seurat$cell_type_v8)
Idents(fb_seurat) <- "cell_type_v8"

saveRDS(fb_seurat, file = "fb_seurat_FINAL_v2.RDS")