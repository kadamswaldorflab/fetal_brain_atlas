# =============================================================================
# plot_harmony_integration.R
#
# Publication-ready figures demonstrating batch effect removal by Harmony.
# Produces three outputs:
#
#   1. Combined figure (sample row + gest_age row + LISI) -- for review
#   2. Sample-only figure + LISI -- main figure
#   3. Gest_age-only figure + LISI -- supplement figure
#
# Layout per single-variable figure:
#   Row 1:  PCA pre | UMAP pre | PCA post | UMAP post
#   Row 2:  iLISI violin -- pre vs post, PCA vs UMAP
#
# LISI computed in pure R via FNN (bundled with Seurat) -- no install needed.
#
# USAGE:
#   source("plot_harmony_integration.R")
#   plot_harmony_integration(fb_seurat,
#                            outfile_combined = "harmony_combined.pdf",
#                            outfile_sample   = "harmony_sample.pdf",
#                            outfile_gest_age = "harmony_gest_age.pdf")
#
# DEPENDENCIES: Seurat, ggplot2, patchwork, dplyr, FNN (bundled with Seurat)
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


# =============================================================================
# HELPERS
# =============================================================================

.get_embedding <- function(seurat_obj, reduction) {
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste0("Reduction '", reduction, "' not found. Available: ",
                paste(names(seurat_obj@reductions), collapse = ", ")))
  }
  as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
}


.make_scatter_panel <- function(embed_df, color_var, color_palette,
                                x_label, y_label, title,
                                pt_size = 0.3, alpha = 0.6,
                                legend = TRUE, aspect_ratio = 1) {
  p <- ggplot(embed_df, aes(x = x, y = y, color = .data[[color_var]])) +
    geom_point(size = pt_size, alpha = alpha, stroke = 0) +
    scale_color_manual(values = color_palette, name = color_var) +
    labs(title = title, x = x_label, y = y_label) +
    # Enforce square (or near-square) aspect ratio so panels aren't tall/skinny.
    # ratio = y_range / x_range * aspect_ratio gives visually correct shape.
    coord_fixed(ratio = aspect_ratio) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text       = element_text(size = 9),
      axis.title      = element_text(size = 10),
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      plot.margin     = margin(4, 4, 4, 4)
    )
  if (!legend) p <- p + theme(legend.position = "none")
  p
}


.build_embed_df <- function(seurat_obj, reduction, dims = 1:2, meta_cols) {
  emb    <- .get_embedding(seurat_obj, reduction)
  coords <- emb[, dims, drop = FALSE]
  colnames(coords) <- c("x", "y")
  meta   <- seurat_obj@meta.data[, meta_cols, drop = FALSE]
  cbind(coords, meta)
}


.make_palette <- function(vals, type = "sample") {
  n <- length(vals)
  if (type == "gest_age") {
    cols <- colorRampPalette(
      c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D")
    )(n)
  } else {
    base <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
      "#984EA3", "#A65628", "#F781BF", "#1ABC9C",
      "#D4AC0D", "#2C3E50", "#E67E22", "#16A085"
    )
    cols <- base[seq_len(min(n, length(base)))]
  }
  setNames(cols, vals)
}


.build_row <- function(panels, color_var, color_palette, pt_size, alpha,
                       aspect_ratio = 1) {
  lapply(seq_along(panels), function(i) {
    p <- panels[[i]]
    .make_scatter_panel(
      embed_df      = p$df,
      color_var     = color_var,
      color_palette = color_palette,
      x_label       = p$x,
      y_label       = p$y,
      title         = p$title,
      pt_size       = pt_size,
      alpha         = alpha,
      legend        = (i == length(panels)),
      aspect_ratio  = aspect_ratio
    )
  })
}


.assemble_figure <- function(scatter_row, lisi_panel, main_title,
                              lisi_height = 0.55) {
  row_patch <- wrap_plots(scatter_row, nrow = 1)
  if (!is.null(lisi_panel)) {
    fig <- row_patch / lisi_panel + plot_layout(heights = c(1, lisi_height))
  } else {
    fig <- row_patch
  }
  fig + plot_annotation(
    title = main_title,
    theme = theme(plot.title = element_text(size = 14, face = "bold",
                                            hjust = 0.5))
  )
}


.save_fig <- function(fig, outfile, width, height) {
  if (!is.null(outfile)) {
    ggsave(outfile, plot = fig, width = width, height = height,
           device = "pdf", limitsize = FALSE)
    message("Saved: ", outfile)
  }
}


# =============================================================================
# LISI -- pure R implementation (no external package required)
#
# Local Inverse Simpson's Index (Korsunsky et al. 2019, Nature Methods).
# For each cell, finds k nearest neighbors in the embedding, then computes
# the inverse Simpson's index over the batch labels of those neighbors.
#
# Score interpretation:
#   ~1          = neighborhood dominated by one sample (poor mixing)
#   ~n_samples  = neighborhood perfectly mixed across all samples
#
# Uses FNN::get.knn which ships with Seurat -- no additional install needed.
# =============================================================================

.ilisi_pure_r <- function(embeddings, batch_labels, k = 90) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package FNN not found. Install with: install.packages('FNN')")
  }
  embeddings <- as.matrix(embeddings)
  n_cells    <- nrow(embeddings)
  labels     <- as.character(batch_labels)
  k          <- min(k, n_cells - 1L)
  message(paste0("    kNN (k=", k, ") for ", n_cells, " cells..."))
  nn_idx <- FNN::get.knn(embeddings, k = k)$nn.index
  vapply(seq_len(n_cells), function(i) {
    props <- table(labels[nn_idx[i, ]]) / k
    1 / sum(props^2)
  }, numeric(1))
}


.compute_lisi_scores <- function(seurat_obj, batch_var, reductions_list) {
  meta      <- as.character(seurat_obj@meta.data[[batch_var]])
  n_batches <- length(unique(meta))
  k         <- min(90L, ncol(seurat_obj) - 1L)
  all_scores <- lapply(reductions_list, function(r) {
    message(paste0("  iLISI: ", r$label, " (", r$stage, ")"))
    emb    <- as.matrix(.get_embedding(seurat_obj, r$name)[, r$dims, drop = FALSE])
    scores <- .ilisi_pure_r(emb, meta, k = k)
    data.frame(lisi_score = scores, embedding = r$label, stage = r$stage)
  })
  result <- do.call(rbind, all_scores)
  message(paste0("  iLISI range: [", round(min(result$lisi_score), 2),
                 ", ", round(max(result$lisi_score), 2),
                 "]  |  perfect mixing = ", n_batches))
  result
}


.make_lisi_panel <- function(lisi_df, batch_var) {
  lisi_df$stage <- factor(lisi_df$stage, levels = c("Pre-Harmony", "Post-Harmony"))
  med_df <- lisi_df %>%
    group_by(embedding, stage) %>%
    summarise(med = median(lisi_score), .groups = "drop")
  ggplot(lisi_df, aes(x = stage, y = lisi_score, fill = stage)) +
    geom_violin(alpha = 0.7, color = NA, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 color = "black", fill = "white", alpha = 0.8) +
    geom_text(data = med_df,
              aes(x = stage, y = med, label = round(med, 2)),
              vjust = -0.6, size = 3.2, fontface = "bold", color = "black") +
    facet_wrap(~ embedding, nrow = 1) +
    scale_fill_manual(values = c("Pre-Harmony"  = "#D6604D",
                                 "Post-Harmony" = "#4393C3")) +
    labs(title = paste0("iLISI by ", batch_var,
                        "  (higher = better sample mixing)"),
         x = NULL, y = "iLISI score") +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text        = element_text(size = 10),
      axis.title.y     = element_text(size = 10),
      strip.text       = element_text(size = 10, face = "bold"),
      strip.background = element_blank(),
      legend.position  = "none",
      panel.spacing    = unit(0.8, "cm")
    )
}


# =============================================================================
# MAIN FUNCTION
# =============================================================================

plot_harmony_integration <- function(
    seurat_obj,
    sample_col          = "sample",
    gest_age_col        = "gest_age",
    reduction_pre_pca   = "pca",
    reduction_post_pca  = "harmony",
    reduction_pre_umap  = "umap_preharmony",
    reduction_post_umap = "umap",
    pca_dims            = 1:2,
    n_lisi_cells        = 5000,
    lisi_pca_dims       = 1:30,
    pt_size             = 0.3,
    alpha               = 0.6,
    aspect_ratio        = 1,
    outfile_combined    = "harmony_combined.pdf",
    outfile_sample      = "harmony_sample.pdf",
    outfile_gest_age    = "harmony_gest_age.pdf",
    width_combined      = 18,
    height_combined     = 16,
    width_single        = 18,
    height_single       = 10,
    save                = TRUE
) {
  for (col in c(sample_col, gest_age_col)) {
    if (!col %in% colnames(seurat_obj@meta.data))
      stop(paste("Metadata column not found:", col))
  }
  for (r in c(reduction_pre_pca, reduction_post_pca,
              reduction_pre_umap, reduction_post_umap)) {
    if (!r %in% names(seurat_obj@reductions))
      stop(paste0("Reduction not found: ", r))
  }

  seurat_obj@meta.data[[gest_age_col]] <-
    as.character(seurat_obj@meta.data[[gest_age_col]])

  samples   <- sort(unique(seurat_obj@meta.data[[sample_col]]))
  gest_ages <- sort(unique(seurat_obj@meta.data[[gest_age_col]]))
  sample_pal   <- .make_palette(samples,   "sample")
  gest_age_pal <- .make_palette(gest_ages, "gest_age")

  meta_cols <- c(sample_col, gest_age_col)
  pc_x <- paste0("PC", pca_dims[1])
  pc_y <- paste0("PC", pca_dims[2])

  panels <- list(
    list(df = .build_embed_df(seurat_obj, reduction_pre_pca,   pca_dims, meta_cols),
         x = pc_x, y = pc_y,         title = "PCA -- Pre-Harmony"),
    list(df = .build_embed_df(seurat_obj, reduction_post_pca,  pca_dims, meta_cols),
         x = pc_x, y = pc_y,         title = "PCA -- Post-Harmony"),
    list(df = .build_embed_df(seurat_obj, reduction_pre_umap,  1:2,      meta_cols),
         x = "UMAP 1", y = "UMAP 2", title = "UMAP -- Pre-Harmony"),
    list(df = .build_embed_df(seurat_obj, reduction_post_umap, 1:2,      meta_cols),
         x = "UMAP 1", y = "UMAP 2", title = "UMAP -- Post-Harmony")
  )

  sample_row   <- .build_row(panels, sample_col,   sample_pal,   pt_size, alpha, aspect_ratio)
  gest_age_row <- .build_row(panels, gest_age_col, gest_age_pal, pt_size, alpha, aspect_ratio)

  message("Computing iLISI scores (pure R)...")
  cells_use <- colnames(seurat_obj)
  if (!is.null(n_lisi_cells) && n_lisi_cells < length(cells_use)) {
    set.seed(42)
    cells_use <- sample(cells_use, n_lisi_cells)
    message("  Subsampled to ", n_lisi_cells, " cells.")
  }
  sub_obj <- seurat_obj[, cells_use]

  lisi_reductions <- list(
    list(name = reduction_pre_pca,   dims = lisi_pca_dims,
         label = "PCA",  stage = "Pre-Harmony"),
    list(name = reduction_post_pca,  dims = lisi_pca_dims,
         label = "PCA",  stage = "Post-Harmony"),
    list(name = reduction_pre_umap,  dims = 1:2,
         label = "UMAP", stage = "Pre-Harmony"),
    list(name = reduction_post_umap, dims = 1:2,
         label = "UMAP", stage = "Post-Harmony")
  )

  # iLISI computed by sample only.
  # gest_age is 1:1 confounded with sample (each animal = unique timepoint)
  # so iLISI by gest_age is mathematically identical to iLISI by sample.
  # Biological preservation of gestational age structure is shown visually
  # in the gest_age UMAP panels rather than as a separate LISI metric.
  lisi_sample <- tryCatch(
    .compute_lisi_scores(sub_obj, sample_col, lisi_reductions),
    error = function(e) { warning("iLISI failed: ", e$message); NULL }
  )
  lisi_panel_sample <- if (!is.null(lisi_sample))
    .make_lisi_panel(lisi_sample, sample_col) else NULL

  combined_fig <- {
    top <- wrap_plots(sample_row, nrow = 1) +
      plot_annotation(title = paste("Colored by", sample_col),
                      theme = theme(plot.title = element_text(size = 12,
                                    face = "bold", hjust = 0)))
    mid <- wrap_plots(gest_age_row, nrow = 1) +
      plot_annotation(title = paste("Colored by", gest_age_col),
                      theme = theme(plot.title = element_text(size = 12,
                                    face = "bold", hjust = 0)))
    if (!is.null(lisi_panel_sample)) {
      top / mid / lisi_panel_sample +
        plot_layout(heights = c(1, 1, 0.6)) +
        plot_annotation(title = "Harmony Batch Effect Correction",
                        theme = theme(plot.title = element_text(size = 15,
                                      face = "bold", hjust = 0.5)))
    } else {
      top / mid +
        plot_annotation(title = "Harmony Batch Effect Correction",
                        theme = theme(plot.title = element_text(size = 15,
                                      face = "bold", hjust = 0.5)))
    }
  }

  # Sample figure: scatter + iLISI (main figure)
  sample_fig <- .assemble_figure(
    sample_row, lisi_panel_sample,
    paste("Harmony Integration --", sample_col)
  )

  # Gest_age figure: scatter only, no LISI (supplement)
  gest_age_fig <- .assemble_figure(
    gest_age_row, NULL,
    paste("Harmony Integration --", gest_age_col)
  )

  if (save) {
    .save_fig(combined_fig, outfile_combined, width_combined, height_combined)
    .save_fig(sample_fig,   outfile_sample,   width_single,   height_single)
    .save_fig(gest_age_fig, outfile_gest_age, width_single,   height_single)
  }

  invisible(list(combined    = combined_fig,
                 sample      = sample_fig,
                 gest_age    = gest_age_fig,
                 lisi_sample = lisi_sample))
}
