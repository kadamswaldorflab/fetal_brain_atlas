# =============================================================================
# plot_umap_genes.R
#
# Plot multi-gene expression on a UMAP, with each gene colored distinctly.
# Cells expressing a gene (> 0) are plotted on top of a grey background of
# all cells. Intended for Seurat objects from the fetal brain atlas.
#
# USAGE:
#   source("plot_umap_genes.R")
#
#   # Basic — auto-assigns colors
#   plot_umap_genes(neuron_subset, genes = c("EOMES", "GAD1", "SATB2"),
#                   outfile = "fig_3b.pdf")
#
#   # Manual colors
#   plot_umap_genes(neuron_subset,
#                   genes  = c("EOMES", "GAD1", "SATB2"),
#                   colors = c("EOMES" = "#2ecc71",
#                              "GAD1"  = "#3498db",
#                              "SATB2" = "#e67e22"),
#                   outfile = "fig_3b.pdf")
#
#   # Return the ggplot object without saving (e.g. for further modification)
#   p <- plot_umap_genes(astro_subset, genes = c("HOPX", "AQP4", "GFAP"),
#                        save = FALSE)
#   p + ggtitle("My custom title")
#
# DEPENDENCIES: Seurat, ggplot2, dplyr, tidyr
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)


# Default color palette — maximally distinct hues for up to 10 genes.
# Order is designed so the first 3 (the most commonly used) are
# bright yellow, blue, and orange — high contrast on both white and dark
# UMAP backgrounds, and clearly distinct from each other.
.DEFAULT_GENE_COLORS <- c(
  "#F1C40F",  # bright yellow   (1st — high visibility on dark backgrounds)
  "#2980B9",  # strong blue     (2nd — maximum contrast with yellow)
  "#E67E22",  # orange          (3rd — distinct warm hue from both above)
  "#27AE60",  # green           (4th)
  "#8E44AD",  # purple          (5th)
  "#E74C3C",  # red             (6th)
  "#16A085",  # teal            (7th)
  "#D35400",  # deep orange     (8th)
  "#2C3E50",  # dark slate      (9th)
  "#F781BF"   # pink            (10th)
)


plot_umap_genes <- function(
    seurat_obj,
    genes,
    colors         = NULL,
    title          = NULL,
    umap_reduction = "umap",
    min_expr       = 0,
    bg_color       = "lightgrey",
    bg_size        = 0.5,
    bg_alpha       = 0.4,
    pt_size        = 0.8,
    pt_alpha       = 0.8,
    legend_pt_size = 4,
    outfile        = NULL,
    width          = 10,
    height         = 8,
    save           = !is.null(outfile)
) {
  # ---------------------------------------------------------------------------
  # Arguments:
  #   seurat_obj     Seurat object (any subset — uses its own UMAP embedding)
  #   genes          Character vector of gene names to plot
  #   colors         Named character vector mapping gene -> hex color.
  #                  If NULL, colors are auto-assigned from .DEFAULT_GENE_COLORS.
  #                  Partial specification is OK — missing genes get auto-colors.
  #   title          Plot title. Defaults to comma-joined gene names.
  #   umap_reduction Name of the UMAP reduction in seurat_obj@reductions.
  #                  Default: "umap"
  #   min_expr       Minimum expression threshold for a cell to be plotted as
  #                  expressing. Default 0 (any detectable expression).
  #   bg_color       Color for background (non-expressing) cells. Default "lightgrey"
  #   bg_size        Point size for background cells. Default 0.5
  #   bg_alpha       Transparency for background cells. Default 0.4
  #   pt_size        Point size for expressing cells. Default 0.8
  #   pt_alpha       Transparency for expressing cells. Default 0.8
  #   legend_pt_size Point size used in the legend. Default 4
  #   outfile        Output filename (e.g. "fig_3b.pdf"). If NULL and save=TRUE,
  #                  defaults to a name derived from the gene list.
  #   width          Figure width in inches. Default 10
  #   height         Figure height in inches. Default 8
  #   save           Whether to save the plot. Defaults to TRUE if outfile is
  #                  provided, FALSE otherwise.
  # ---------------------------------------------------------------------------

  # ---- Input validation ----
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (length(genes) == 0) {
    stop("`genes` must be a non-empty character vector.")
  }

  missing_genes <- genes[!genes %in% rownames(seurat_obj)]
  if (length(missing_genes) > 0) {
    stop(paste("Gene(s) not found in seurat_obj:",
               paste(missing_genes, collapse = ", ")))
  }

  if (!(umap_reduction %in% names(seurat_obj@reductions))) {
    stop(paste0("Reduction '", umap_reduction, "' not found. ",
                "Available: ", paste(names(seurat_obj@reductions), collapse = ", ")))
  }

  # ---- UMAP coordinates ----
  umap_coords <- as.data.frame(
    seurat_obj@reductions[[umap_reduction]]@cell.embeddings
  )

  # Normalize column names to lowercase umap_1 / umap_2 regardless of Seurat
  # version (which may produce UMAP_1/UMAP_2 or umap_1/umap_2)
  colnames(umap_coords) <- tolower(colnames(umap_coords))
  if (!all(c("umap_1", "umap_2") %in% colnames(umap_coords))) {
    stop(paste("Could not find UMAP_1/UMAP_2 columns. Found:",
               paste(colnames(umap_coords), collapse = ", ")))
  }
  umap_coords <- umap_coords[, c("umap_1", "umap_2")]

  # ---- Expression data ----
  expr_data <- FetchData(object = seurat_obj, vars = genes)
  plot_data <- cbind(umap_coords, expr_data)

  # Pivot to long format and filter to expressing cells
  long_data <- plot_data %>%
    pivot_longer(
      cols      = all_of(genes),
      names_to  = "gene",
      values_to = "expression"
    ) %>%
    filter(expression > min_expr) %>%
    mutate(gene = factor(gene, levels = genes))  # preserve gene order

  # Shuffle so no single gene systematically occludes others
  set.seed(42)
  shuffled_data <- long_data %>% sample_frac(1)

  # ---- Color assignment ----
  # Start with defaults, then overlay any user-supplied colors
  auto_colors <- setNames(
    .DEFAULT_GENE_COLORS[seq_along(genes)],
    genes
  )
  if (!is.null(colors)) {
    # Validate that user-supplied names are in genes
    bad_names <- names(colors)[!names(colors) %in% genes]
    if (length(bad_names) > 0) {
      warning(paste("Color names not in `genes` (ignored):",
                    paste(bad_names, collapse = ", ")))
    }
    # Override auto-colors with user-supplied ones
    auto_colors[names(colors)[names(colors) %in% genes]] <-
      colors[names(colors) %in% genes]
  }
  final_colors <- auto_colors

  # ---- Title ----
  if (is.null(title)) {
    title <- paste(genes, collapse = ", ")
  }

  # ---- Build plot ----
  p <- ggplot() +

    # Background: all cells in grey
    geom_point(
      data  = umap_coords,
      aes(x = umap_1, y = umap_2),
      color = bg_color,
      size  = bg_size,
      alpha = bg_alpha
    ) +

    # Foreground: expressing cells, colored by gene identity
    geom_point(
      data = shuffled_data,
      aes(x = umap_1, y = umap_2, color = gene),
      size  = pt_size,
      alpha = pt_alpha
    ) +

    scale_color_manual(values = final_colors) +

    guides(
      color = guide_legend(
        title      = "Gene",
        override.aes = list(size = legend_pt_size, alpha = pt_alpha)
      )
    ) +

    labs(
      title = title,
      x     = "UMAP 1",
      y     = "UMAP 2"
    ) +

    theme_minimal() +
    theme(
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "black"),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 10),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 11),
      legend.key       = element_blank()
    )

  # ---- Save ----
  if (save) {
    if (is.null(outfile)) {
      gene_str <- paste(genes[seq_len(min(3, length(genes)))], collapse = "_")
      outfile  <- paste0("umap_genes_", gene_str, ".pdf")
    }
    ggsave(outfile, plot = p, width = width, height = height, units = "in")
    message("Saved: ", outfile)
  }

  invisible(p)
}
