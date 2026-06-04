# =============================================================================
# plot_umap_clusters.R
#
# Generate a labeled UMAP cluster plot from a Seurat object, with repelled
# shadow labels at cluster centroids. Wraps DimPlot_scCustom + ggrepel.
#
# USAGE:
#   source("plot_umap_clusters.R")
#
#   # Basic — varibow palette, auto color count
#   plot_umap_clusters(fb_seurat, outfile = "umap_clusters.pdf")
#
#   # Custom palette and number of colors
#   plot_umap_clusters(fb_seurat,
#                      palette    = "polychrome",
#                      num_colors = 26,
#                      outfile    = "umap_clusters.pdf")
#
#   # Pass your own hex colors directly
#   my_colors <- c("#E64B35", "#4DBBD5", "#00A087", ...)
#   plot_umap_clusters(fb_seurat,
#                      custom_colors = my_colors,
#                      outfile = "umap_clusters.pdf")
#
#   # Return ggplot object without saving
#   p <- plot_umap_clusters(fb_seurat, save = FALSE)
#
# DEPENDENCIES: Seurat, scCustomize, ggplot2, dplyr, ggrepel
# =============================================================================

library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(ggrepel)


plot_umap_clusters <- function(
    seurat_obj,
    reduction     = "umap",
    group_by      = NULL,
    palette       = "varibow",
    num_colors    = NULL,
    shuffle_pal   = FALSE,
    color_seed    = 42,
    custom_colors = NULL,
    pt_size       = 0.5,
    label         = TRUE,
    label_size    = 4,
    label_face    = "bold",
    label_color   = "black",
    label_bg      = "white",
    label_bg_r    = 0.15,
    box_padding   = 0.5,
    point_padding = 0.3,
    max_overlaps  = Inf,
    seg_color     = "grey50",
    seg_size      = 0.3,
    min_seg_len   = 0.2,
    title         = NULL,
    width         = 15,
    height        = 10,
    outfile       = NULL,
    save          = !is.null(outfile)
) {
  # ---------------------------------------------------------------------------
  # Arguments:
  #   seurat_obj     Seurat object
  #   reduction      Name of the reduction to plot. Default: "umap"
  #   group_by       Metadata column to color by. Default: NULL (uses active
  #                  identity set by Idents())
  #   palette        scCustomize palette name passed to
  #                  DiscretePalette_scCustomize(). Ignored if custom_colors
  #                  is supplied. Common options: "varibow", "polychrome",
  #                  "stepped", "ditto_seq", "turbo", "greenblue".
  #                  Default: "varibow"
  #   num_colors     Number of colors to generate. If NULL, auto-detected from
  #                  the number of unique identities in the plot.
  #   shuffle_pal    Whether to shuffle the palette order. Default: FALSE
  #   color_seed     Random seed for color shuffling. Default: 42
  #   custom_colors  Optional character vector of hex colors. If supplied,
  #                  overrides palette/num_colors entirely. Length must be >=
  #                  number of unique identities.
  #   pt_size        Scatter point size. Default: 0.5
  #   label          Whether to add repelled cluster labels. Default: TRUE
  #   label_size     Font size for cluster labels (ggplot units). Default: 4
  #   label_face     Font face for labels: "bold", "plain", "italic". Default: "bold"
  #   label_color    Label text color. Default: "black"
  #   label_bg       Label shadow/halo background color. Default: "white"
  #   label_bg_r     Radius of label shadow. Default: 0.15
  #   box_padding    ggrepel box padding. Default: 0.5
  #   point_padding  ggrepel point padding. Default: 0.3
  #   max_overlaps   ggrepel max overlaps. Default: Inf (show all labels)
  #   seg_color      Color of repel leader lines. Default: "grey50"
  #   seg_size       Width of repel leader lines. Default: 0.3
  #   min_seg_len    Minimum leader line length before it's drawn. Default: 0.2
  #   title          Plot title. Default: NULL (no title)
  #   width          Output width in inches. Default: 15
  #   height         Output height in inches. Default: 10
  #   outfile        Output filename (e.g. "umap.pdf"). If NULL and save=TRUE,
  #                  auto-generates from palette name.
  #   save           Whether to save. Defaults TRUE if outfile given, else FALSE.
  # ---------------------------------------------------------------------------

  # ---- Input validation ----
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste0("Reduction '", reduction, "' not found. Available: ",
                paste(names(seurat_obj@reductions), collapse = ", ")))
  }
  if (!is.null(group_by) && !(group_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste0("group_by column '", group_by, "' not found in metadata."))
  }

  # ---- Determine number of identities ----
  if (!is.null(group_by)) {
    n_idents <- length(unique(seurat_obj@meta.data[[group_by]]))
  } else {
    n_idents <- length(levels(Idents(seurat_obj)))
    if (n_idents == 0) n_idents <- length(unique(Idents(seurat_obj)))
  }

  if (is.null(num_colors)) num_colors <- n_idents

  if (!is.null(custom_colors) && length(custom_colors) < n_idents) {
    stop(paste0("custom_colors has ", length(custom_colors), " colors but ",
                n_idents, " identities found. Provide at least ", n_idents,
                " colors."))
  }

  # ---- Build color vector ----
  if (!is.null(custom_colors)) {
    colors_use <- custom_colors[seq_len(n_idents)]
  } else {
    colors_use <- DiscretePalette_scCustomize(
      num_colors  = num_colors,
      palette     = palette,
      shuffle_pal = shuffle_pal
    )
  }

  # ---- Base DimPlot ----
  dimplot_args <- list(
    seurat_obj,
    reduction  = reduction,
    label      = FALSE,        # we draw our own labels via ggrepel
    pt.size    = pt_size,
    repel      = TRUE,
    colors_use = colors_use,
    color_seed = color_seed
  )
  if (!is.null(group_by)) dimplot_args$group.by <- group_by

  p <- do.call(DimPlot_scCustom, dimplot_args)

  # ---- Compute label positions (median of each cluster) ----
  if (label) {
    # Detect which column in p$data holds the grouping variable.
    # When group.by is NULL, DimPlot uses "ident". When group.by is set,
    # DimPlot names the column after the metadata column (e.g. "cell_type_v8").
    plot_cols <- colnames(p$data)
    x_col <- grep("_1$", plot_cols, value = TRUE)[1]
    y_col <- grep("_2$", plot_cols, value = TRUE)[1]

    if (!is.null(group_by) && group_by %in% plot_cols) {
      group_col <- group_by
    } else if ("ident" %in% plot_cols) {
      group_col <- "ident"
    } else {
      warning("Could not detect grouping column in plot data — labels skipped. ",
              "Columns found: ", paste(plot_cols, collapse = ", "))
      label <- FALSE
      group_col <- NULL
    }

    if (is.na(x_col) || is.na(y_col)) {
      warning("Could not detect UMAP coordinate columns — labels skipped. ",
              "Columns found: ", paste(plot_cols, collapse = ", "))
      label <- FALSE
    }
  }

  if (label) {
    label_data <- p$data %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        x = median(.data[[x_col]]),
        y = median(.data[[y_col]]),
        .groups = "drop"
      )

    p <- p + geom_text_repel(
      data              = label_data,
      aes(x = x, y = y, label = .data[[group_col]]),
      colour            = label_color,
      bg.color          = label_bg,
      bg.r              = label_bg_r,
      size              = label_size,
      fontface          = label_face,
      max.overlaps      = max_overlaps,
      box.padding       = box_padding,
      point.padding     = point_padding,
      segment.color     = seg_color,
      segment.size      = seg_size,
      min.segment.length = min_seg_len
    )
  }

  # ---- Axis labels and title ----
  p <- p +
    xlab("UMAP 1") +
    ylab("UMAP 2")

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  # ---- Save ----
  if (save) {
    if (is.null(outfile)) {
      outfile <- paste0("umap_clusters_", palette, ".pdf")
    }
    ggsave(outfile, plot = p, width = width, height = height, device = "pdf")
    message("Saved: ", outfile)
  }

  invisible(p)
}
