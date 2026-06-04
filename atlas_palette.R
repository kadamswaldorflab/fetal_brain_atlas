# =============================================================================
# atlas_palette.R
#
# Cell type color palette for the fetal brain SVZ atlas (v8 annotations).
# Identical to CELL_TYPE_PALETTE in the scVelo scripts — same cell type
# always gets the same color across Python and R figures.
#
# Colors are grouped by biological lineage:
#   Radial glia    → reds
#   Progenitors    → oranges
#   Astrocytes     → blues (light → dark)
#   OL lineage     → greens (light → dark)
#   Neurons        → purples
#   CGE            → salmon-pink
#   MGE            → magentas
#   Microglia      → blue-grays
#   Stromal        → violets / grays
#   Non-neural     → yellows / orange / teal
#
# USAGE:
#   source("atlas_palette.R")
#
#   # Use with plot_umap_clusters():
#   plot_umap_clusters(fb_seurat,
#                      custom_colors = get_atlas_colors(fb_seurat),
#                      outfile = "umap.pdf")
#
#   # Use with plot_umap_genes() — pass individual gene colors:
#   plot_umap_genes(astro_subset,
#                   genes  = c("HOPX", "AQP4", "GFAP"),
#                   colors = c("HOPX" = ATLAS_PALETTE[["APC"]],
#                              "AQP4" = ATLAS_PALETTE[["AS1"]],
#                              "GFAP" = ATLAS_PALETTE[["AS2"]]),
#                   outfile = "astro_markers.pdf")
# =============================================================================

ATLAS_PALETTE <- c(
  # --- Radial glia ---
  "vRG"     = "#145A32",   # dark forest green (goes with OL lineage)
  "oRG-G"   = "#1F618D",   # steel blue

  # --- Progenitor / intermediate states ---
  "APC"     = "#D6EAF8",   # very pale blue (matches astrocyte lineage)
  "UP/I"    = "#F0B27A",   # peach (undifferentiated/intermediate)
  "IPC"     = "#784212",   # dark brown-orange

  # --- Astrocyte lineage (light → dark blue, APC is lightest) ---
  "AS0"     = "#AED6F1",   # pale sky blue
  "AS1"     = "#5DADE2",   # cornflower blue
  "AS2"     = "#2874A6",   # medium-dark blue

  # --- Oligodendrocyte lineage (light → dark green, vRG is darkest) ---
  "OPC"     = "#A9DFBF",   # pale green
  "COP1"    = "#52BE80",   # medium-light green
  "COP2"    = "#27AE60",   # medium green
  "OL1"     = "#1E8449",   # dark green
  "OL2"     = "#0B5345",   # very dark green

  # --- Excitatory neurons (purples) ---
  "ExN"     = "#7D3C98",   # medium purple
  "ImmN"    = "#4A235A",   # dark purple

  # --- CGE interneurons ---
  "CGE"     = "#F1948A",   # salmon-pink

  # --- MGE interneurons ---
  "MGE"     = "#D81B60",   # magenta
  "MGE-IN" = "#880E4F",   # dark magenta (SST-positive subset)

  # --- Microglia ---
  "MG0"     = "#B2BABB",   # light blue-gray (homeostatic)
  "MG1"     = "#5D6D7E",   # dark blue-gray (reactive)

  # --- Other immune / myeloid ---
  "BAM"     = "#1ABC9C",   # teal

  # --- Non-neural populations ---
  "ChP"     = "#F4D03F",   # yellow (choroid plexus)
  "Ep"      = "#F0A500",   # gold (ependymal)
  "EC"      = "#FF7043",   # deep orange (endothelial)

  # --- Stromal populations (warm sandy family) ---
  "FIB"     = "#D5C4A1",   # pale sandy beige
  "PC"      = "#9E8B6A",   # medium sand
  "MPC"     = "#6D5D44"    # dark khaki/umber (mesenchymal progenitor)
)


# =============================================================================
# HIGH-CONTRAST SUBSET PALETTE
#
# For subset UMAPs (astrocyte lineage, OL lineage, etc.) where the lineage-
# based atlas colors are too similar at smaller plot scale. Assigns maximally
# distinct hues regardless of biological family — purely for visual separability.
#
# Usage:
#   plot_umap_clusters(astro_subset,
#                      custom_colors = get_subset_colors(astro_subset),
#                      outfile = "astro_umap.pdf")
# =============================================================================
SUBSET_PALETTE <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#FF7F00",  # orange
  "#984EA3",  # purple
  "#A65628",  # brown
  "#F781BF",  # pink
  "#1ABC9C",  # teal
  "#D4AC0D",  # gold (replaces yellow — visible on white)
  "#2C3E50"   # dark slate
)


#' Helper: extract identity order that matches what DimPlot will use.
#'
#' DimPlot_scCustom assigns colors positionally by factor level order.
#' This must match exactly — even one cell type out of position swaps colors.
#' Priority: (1) factor levels of the column, (2) Idents() levels,
#' (3) alphabetical sort as last resort.
.get_ident_order <- function(seurat_obj, group_by = NULL) {
  if (!is.null(group_by)) {
    col <- seurat_obj@meta.data[[group_by]]
    if (is.factor(col)) {
      # Use actual factor levels — this is what DimPlot reads
      return(levels(col))
    } else {
      # Column is character/numeric — sort alphabetically (Seurat default)
      return(sort(unique(as.character(col))))
    }
  } else {
    idents <- levels(Idents(seurat_obj))
    if (length(idents) == 0) {
      idents <- sort(unique(as.character(Idents(seurat_obj))))
    }
    return(idents)
  }
}


#' Get high-contrast colors ordered to match a Seurat object's identity levels.
#'
#' Uses SUBSET_PALETTE instead of ATLAS_PALETTE — best for subset UMAPs with
#' few cell types where lineage-based blues/greens are too similar.
#'
#' Color order matches DimPlot_scCustom's positional assignment exactly.
#' If group_by column is a factor in your Seurat object, set its levels before
#' calling this function and they will be respected.
#'
#' @param seurat_obj  A Seurat object.
#' @param group_by    Metadata column to use. If NULL, uses active Idents().
#'
#' @return Named character vector of hex colors, one per identity level.
get_subset_colors <- function(seurat_obj, group_by = NULL) {
  idents <- .get_ident_order(seurat_obj, group_by)

  if (length(idents) > length(SUBSET_PALETTE)) {
    warning(paste("More identities than palette colors — colors will recycle.",
                  "Consider using get_atlas_colors() instead."))
  }

  colors <- SUBSET_PALETTE[((seq_along(idents) - 1) %% length(SUBSET_PALETTE)) + 1]
  names(colors) <- idents

  message(paste("Subset color order:",
                paste(names(colors), collapse = " > ")))
  colors
}


#' Get atlas colors ordered to match a Seurat object's identity levels.
#'
#' DimPlot and DimPlot_scCustom assign colors positionally — the first color
#' goes to the first factor level, the second color to the second level, etc.
#' This function returns ATLAS_PALETTE values in the correct order for a given
#' Seurat object, with a fallback gray for any cell type not in the palette.
#'
#' Color order matches DimPlot_scCustom's positional assignment exactly.
#' If group_by column is a factor in your Seurat object, set its levels before
#' calling this function and they will be respected.
#'
#' @param seurat_obj  A Seurat object.
#' @param group_by    Metadata column to use. If NULL, uses active Idents().
#' @param fallback    Hex color to assign to unrecognized cell types. Default "#999999".
#'
#' @return Named character vector of hex colors, one per identity level.
#'
#' @examples
#'   colors <- get_atlas_colors(fb_seurat)
#'   plot_umap_clusters(fb_seurat, custom_colors = colors, outfile = "umap.pdf")
get_atlas_colors <- function(seurat_obj, group_by = NULL, fallback = "#999999") {
  idents <- .get_ident_order(seurat_obj, group_by)

  missing <- idents[!idents %in% names(ATLAS_PALETTE)]
  if (length(missing) > 0) {
    warning(paste("Cell types not in ATLAS_PALETTE (will use fallback color):",
                  paste(missing, collapse = ", ")))
  }

  colors <- ifelse(idents %in% names(ATLAS_PALETTE),
                   ATLAS_PALETTE[idents],
                   fallback)
  names(colors) <- idents
  colors
}
