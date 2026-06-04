"""
scVelo trajectory analysis on the FULL fetal brain atlas (SVZ, late gestation macaque).

Loads all loom files, attaches Seurat metadata + UMAP, computes stochastic-mode
RNA velocity, and writes plots + an h5ad checkpoint for downstream CellRank.

Run AFTER exporting Seurat metadata to CSV with columns:
    barcode (index), sample, <cell_type_col>, gest_age, umap_1, umap_2

Loom filename convention:  {SampleID}.loom   e.g. CTRL13.loom, SAL7.loom
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv
import scipy.sparse as sp

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# =============================================================================
# CELL TYPE COLOR PALETTE
# Lineage-grouped, consistent across full atlas and subset scripts.
# Colors are assigned by biological lineage family so related cell types share
# a hue, making the figure easier to interpret at a glance.
# Any cell type not listed here falls back to '#999999' (mid-gray).
# =============================================================================
CELL_TYPE_PALETTE = {
    # --- Radial glia ---
    # vRG → green family (goes with oligodendrocytes)
    'vRG':     '#145A32',   # dark forest green
    # oRG-G → blue family
    'oRG-G':   '#1F618D',   # steel blue (distinct from AS sky blues)

    # --- Progenitor / intermediate states ---
    # APC → blue family (matches astrocyte lineage)
    'APC':     '#D6EAF8',   # very pale blue (lightest in astrocyte family)
    'UP/I':    '#F0B27A',   # peach (undifferentiated/intermediate)
    'IPC':     '#784212',   # dark brown-orange

    # --- Astrocyte lineage (light → dark blue, APC is lightest) ---
    'AS0':     '#AED6F1',   # pale sky blue
    'AS1':     '#5DADE2',   # cornflower blue
    'AS2':     '#2874A6',   # medium-dark blue

    # --- Oligodendrocyte lineage (light → dark green, vRG is darkest) ---
    'OPC':     '#A9DFBF',   # pale green
    'COP1':    '#52BE80',   # medium-light green
    'COP2':    '#27AE60',   # medium green
    'OL1':     '#1E8449',   # dark green
    'OL2':     '#0B5345',   # very dark green

    # --- Excitatory neurons (purples) ---
    'ExN':     '#7D3C98',   # medium purple
    'ImmN':    '#4A235A',   # dark purple

    # --- CGE interneurons ---
    'CGE':     '#F1948A',   # salmon-pink

    # --- MGE interneurons ---
    'MGE':     '#D81B60',   # magenta
    'MGE-SST': '#880E4F',   # dark magenta (SST-positive subset)

    # --- Microglia ---
    'MG0':     '#B2BABB',   # light blue-gray (homeostatic)
    'MG1':     '#5D6D7E',   # dark blue-gray (reactive)

    # --- Other immune / myeloid ---
    'BAM':     '#1ABC9C',   # teal

    # --- Non-neural populations ---
    'ChP':     '#F4D03F',   # yellow (choroid plexus)
    'Ep':      '#F0A500',   # gold (ependymal)
    'EC':      '#FF7043',   # deep orange (endothelial)

    # --- Stromal populations (same warm sandy family) ---
    'FIB':     '#D5C4A1',   # pale sandy beige
    'PC':      '#9E8B6A',   # medium sand
    'MPC':     '#6D5D44',   # dark khaki/umber (mesenchymal progenitor)
}


def get_palette_for_adata(adata, cell_type_col):
    """
    Return ordered list of hex colors matching adata.obs[cell_type_col].cat.categories.
    Falls back to '#999999' for any cell type not in CELL_TYPE_PALETTE.
    """
    if not hasattr(adata.obs[cell_type_col], 'cat'):
        adata.obs[cell_type_col] = adata.obs[cell_type_col].astype('category')
    cats = list(adata.obs[cell_type_col].cat.categories)
    return [CELL_TYPE_PALETTE.get(c, '#999999') for c in cats]


# =============================================================================
# CONFIG
# =============================================================================
# Directory containing per-sample .loom files produced by velocyto run10x.
# Each file should be named {SAMPLE_ID}.loom (e.g. CTRL13.loom, SAL7.loom).
LOOM_DIR = '/path/to/loom_files'

OUTPUT_DIR = 'output_full'

# CSV exported from the Seurat object (see fb_v4.R: save_subset_for_scvelo).
# Required columns: barcode (index), sample, <CELL_TYPE_COL>, umap_1, umap_2.
METADATA_FILE = 'full_metadata.csv'
CELL_TYPE_COL = 'cell_type_v8'

os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def extract_sample_id(filename):
    """
    Extract sample ID from loom filename.

    Filename format: {CTRL|SAL}{number}.loom
    Examples: CTRL13.loom -> CTRL13, SAL7.loom -> SAL7
    """
    basename = os.path.basename(filename).replace('.loom', '')
    match = re.match(r'^(CTRL|SAL)(\d+)$', basename)
    if match:
        return f"{match.group(1)}{match.group(2)}"
    return None


def load_all_loom_files(loom_dir):
    """Load all loom files and return concatenated AnnData with cleaned barcodes."""
    print("=" * 60)
    print("LOADING ALL LOOM FILES")
    print("=" * 60)

    loom_files = sorted(glob.glob(f'{loom_dir}/*.loom'))
    print(f"\nFound {len(loom_files)} loom files")

    adatas = []
    for loom_file in loom_files:
        sample_name = os.path.basename(loom_file).replace('.loom', '')
        sample_id = extract_sample_id(loom_file)

        if sample_id is None:
            print(f"WARNING: could not parse sample ID from {sample_name}, skipping")
            continue

        print(f"Loading {sample_name} -> {sample_id}...")
        adata_loom = sc.read_loom(loom_file)
        adata_loom.obs['loom_sample'] = sample_name
        adata_loom.obs['sample_id'] = sample_id
        adata_loom.var_names_make_unique()
        adatas.append(adata_loom)

    print("\nConcatenating all loom files...")
    if len(adatas) == 1:
        adata_all = adatas[0]
    else:
        adata_all = sc.concat(adatas, axis=0, join='outer')

    print(f"Combined loom data: {adata_all.shape[0]} cells, {adata_all.shape[1]} genes")

    # velocyto barcode format: "sample_name:BARCODEx" -> "BARCODE-1"
    def extract_barcode(full_name):
        parts = full_name.split(':')
        if len(parts) == 2:
            return parts[1].rstrip('x') + '-1'
        return full_name

    adata_all.obs['barcode_clean'] = [extract_barcode(x) for x in adata_all.obs_names]

    print(f"\nBarcode conversion example:")
    print(f"  Original: {adata_all.obs_names[0]}")
    print(f"  Cleaned:  {adata_all.obs['barcode_clean'].iloc[0]}")

    return adata_all


def attach_seurat_metadata(adata_all, metadata_file):
    """
    Match loom cells against Seurat metadata via {sample_id}_{barcode} keys.
    Returns subsetted AnnData with metadata + UMAP attached.
    """
    metadata = pd.read_csv(metadata_file, index_col='barcode')
    print(f"Loaded metadata: {metadata.shape[0]} cells")

    # Seurat barcode format expected: {sample}_{BARCODE-1}
    metadata['barcode_clean'] = metadata.index.str.split('_').str[-1]
    metadata['sample_id'] = metadata['sample']

    print(f"Example metadata barcode: {metadata.index[0]}")
    print(f"  Extracted barcode: {metadata['barcode_clean'].iloc[0]}")
    print(f"  Sample ID (raw):   {metadata['sample_id'].iloc[0]}")

    # Normalize case to lowercase on both sides before matching.
    # Loom filenames are all-caps (CTRL13, SAL7) but Seurat stores
    # sample IDs in title case (Ctrl13, Sal7) — lowercase both to reconcile.
    metadata['match_key'] = (metadata['sample_id'].str.lower()
                             + '_' + metadata['barcode_clean'])
    adata_all.obs['match_key'] = (adata_all.obs['sample_id'].str.lower()
                                  + '_' + adata_all.obs['barcode_clean'])

    print(f"  Match key (normalized): {metadata['match_key'].iloc[0]}")

    common_keys = set(metadata['match_key']).intersection(set(adata_all.obs['match_key']))
    print(f"\nFound {len(common_keys)} matching cells")

    metadata_matched = metadata[metadata['match_key'].isin(common_keys)].copy()
    match_key_to_obs = dict(zip(adata_all.obs['match_key'], adata_all.obs_names))
    adata_obs_names = [match_key_to_obs[key] for key in metadata_matched['match_key']]

    adata = adata_all[adata_obs_names, :].copy()
    assert len(adata) == len(metadata_matched)

    pct = 100 * len(metadata_matched) / len(metadata)
    print(f"Matched {len(metadata_matched)} / {len(metadata)} cells ({pct:.1f}%)")

    if len(metadata_matched) == 0:
        raise RuntimeError("No cells matched between loom and metadata.")

    # Use Seurat barcodes as obs_names
    adata.obs_names = metadata_matched.index.values

    # Transfer all metadata columns
    for col in metadata_matched.columns:
        if col not in ('barcode_clean', 'sample_id', 'match_key'):
            adata.obs[col] = metadata_matched[col].values

    # Attach UMAP
    if 'umap_1' in metadata_matched.columns and 'umap_2' in metadata_matched.columns:
        adata.obsm['X_umap'] = metadata_matched[['umap_1', 'umap_2']].values
        print("Added UMAP coordinates from Seurat")
    else:
        raise RuntimeError("Metadata missing umap_1/umap_2 columns.")

    return adata


def diagnostic_check_counts(adata):
    """Sanity-check that spliced/unspliced layers are raw integer counts."""
    spliced_raw = adata.layers['spliced']
    unspliced_raw = adata.layers['unspliced']
    spliced_vals = spliced_raw.data if sp.issparse(spliced_raw) else spliced_raw.flatten()
    unspliced_vals = unspliced_raw.data if sp.issparse(unspliced_raw) else unspliced_raw.flatten()
    total = adata.shape[0] * adata.shape[1]

    print("\n--- PRE-NORMALIZATION DIAGNOSTIC ---")
    print(f"Cells: {adata.shape[0]}, Genes: {adata.shape[1]}")
    print(f"Spliced   nonzero mean: {spliced_vals.mean():.2f}  "
          f"max: {spliced_vals.max():.0f}  "
          f"% nonzero: {100*len(spliced_vals)/total:.2f}%  "
          f"int: {np.all(spliced_vals == spliced_vals.astype(int))}")
    print(f"Unspliced nonzero mean: {unspliced_vals.mean():.2f}  "
          f"max: {unspliced_vals.max():.0f}  "
          f"% nonzero: {100*len(unspliced_vals)/total:.2f}%  "
          f"int: {np.all(unspliced_vals == unspliced_vals.astype(int))}")
    print("--- END DIAGNOSTIC ---\n")


def filter_invalid_velocity_cells(adata):
    """
    Remove cells with no finite velocity values across velocity genes.

    KEY INSIGHT (from placenta v2 script): scVelo only computes velocity for a
    subset of genes (velocity genes). All other columns in the velocity layer
    are NaN BY DESIGN. Using .all(axis=1) would filter out every cell. Instead,
    identify velocity genes (cols with any finite value) and filter cells based
    only on those.
    """
    velocity_layer = adata.layers['velocity']
    vel_dense = velocity_layer.toarray() if sp.issparse(velocity_layer) else np.array(velocity_layer)

    velocity_gene_mask = np.isfinite(vel_dense).any(axis=0)
    n_velocity_genes = velocity_gene_mask.sum()
    print(f"Velocity genes identified: {n_velocity_genes} / {vel_dense.shape[1]}")

    if n_velocity_genes == 0:
        raise RuntimeError("No velocity genes — splicing signal too weak.")

    vel_subset = vel_dense[:, velocity_gene_mask]
    frac_finite = np.isfinite(vel_subset).mean(axis=1)
    print(f"Per-cell finite velocity fraction — "
          f"mean: {frac_finite.mean():.3f}  "
          f"median: {np.median(frac_finite):.3f}  "
          f"min: {frac_finite.min():.3f}")

    valid_mask = np.isfinite(vel_subset).any(axis=1)
    n_invalid = (~valid_mask).sum()

    if n_invalid > 0:
        pct = 100 * n_invalid / len(adata)
        print(f"Removing {n_invalid} cells with no finite velocity ({pct:.1f}%)")
        adata = adata[valid_mask, :].copy()
        print(f"Remaining cells: {len(adata)}")
    else:
        print("All cells have valid velocity values.")

    return adata


# =============================================================================
# PLOTTING
# =============================================================================
def _save(fig_path):
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_path}")


def make_publication_celltype_stream(adata, cell_type_col, output_dir,
                                     prefix, velocity_mode='stochastic'):
    """
    Publication-quality velocity stream plot colored by cell type.

    Outputs:
      {prefix}_stream_celltype_pub.pdf  — primary figure (vector, 300 dpi)
      {prefix}_stream_celltype_pub.svg  — editable vector copy

    Design decisions:
    - No internal title (use figure caption)
    - Lineage-grouped color palette from CELL_TYPE_PALETTE
    - Larger points, thicker stream lines for print legibility
    - Legend in right margin, 10pt font
    - Axis labels styled as 'UMAP 1' / 'UMAP 2' (no underscore)
    - Spines cleaned to left + bottom only
    """
    print("\n--- Publication stream plot (cell type) ---")

    umap = adata.obsm['X_umap']
    xlim = (umap[:, 0].min() - 1, umap[:, 0].max() + 1)
    ylim = (umap[:, 1].min() - 1, umap[:, 1].max() + 1)

    # Register colors in adata so scVelo picks them up
    colors = get_palette_for_adata(adata, cell_type_col)
    adata.uns[f'{cell_type_col}_colors'] = colors

    n_cats = len(adata.obs[cell_type_col].cat.categories)

    # Scale figure width with number of cell types so legend fits
    fig_w = 10 + max(0, (n_cats - 12) * 0.15)
    fig_h = 8.5

    for fmt in ('pdf', 'svg'):
        try:
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))

            scv.pl.velocity_embedding_stream(
                adata,
                basis='umap',
                color=cell_type_col,
                legend_loc='right margin',
                legend_fontsize=10,
                legend_fontweight='normal',
                size=40,               # larger scatter points for print
                alpha=0.85,
                linewidth=1.0,         # slightly thicker for visibility
                density=2.0,           # increased: more arrows on full atlas
                arrow_size=1.5,        # slightly larger arrowheads
                title='',              # no internal title — use caption
                ax=ax,
                show=False,
            )

            # Axis styling
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel('UMAP 1', fontsize=12, labelpad=6)
            ax.set_ylabel('UMAP 2', fontsize=12, labelpad=6)
            ax.tick_params(labelsize=10)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            out_path = f'{output_dir}/{prefix}_stream_celltype_pub.{fmt}'
            fig.savefig(out_path, dpi=300, bbox_inches='tight',
                        format=fmt)
            plt.close(fig)
            print(f"Saved: {out_path}")

        except Exception as e:
            print(f"  Publication stream ({fmt}) failed: {e}")
            plt.close()


def make_plots(adata, cell_type_col, output_dir, prefix='full_atlas',
               velocity_mode='stochastic'):
    """Generate stream, arrow, confidence, and pseudotime plots."""
    print("\nGenerating plots...")

    umap = adata.obsm['X_umap']
    xlim = (umap[:, 0].min() - 1, umap[:, 0].max() + 1)
    ylim = (umap[:, 1].min() - 1, umap[:, 1].max() + 1)
    print(f"Plot limits: x={xlim}, y={ylim}")

    print("Recomputing velocity embedding for visualization...")
    scv.tl.velocity_embedding(adata, basis='umap')

    # ---- PUBLICATION STREAM PLOT (cell type, PDF + SVG) ----
    make_publication_celltype_stream(
        adata, cell_type_col, output_dir, prefix, velocity_mode
    )

    # ---- DIAGNOSTIC STREAM PLOTS (SVG) ----
    print("\n--- Diagnostic stream plots (SVG) ---")
    plot_specs = [
        (cell_type_col, 'celltype', 'Cell Type'),
        ('sample', 'sample', 'Sample'),
    ]
    if 'gest_age' in adata.obs.columns:
        plot_specs.append(('gest_age', 'gest_age', 'Gestational Age'))

    for color, suffix, title in plot_specs:
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            scv.pl.velocity_embedding_stream(
                adata, basis='umap', color=color,
                legend_loc='right margin',
                title=f'{prefix} ({velocity_mode}) — Velocity Stream by {title}',
                legend_fontsize=6 if color == 'sample' else 8,
                ax=ax, show=False,
            )
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            _save(f'{output_dir}/{prefix}_stream_{suffix}.svg')
        except Exception as e:
            print(f"Error stream by {color}: {e}")
            plt.close()

    # ---- ARROW PLOTS (PDF) ----
    print("\n--- Arrow plots (PDF) ---")
    for color, suffix, title in plot_specs:
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            scv.pl.velocity_embedding(
                adata, basis='umap', color=color,
                arrow_length=3, arrow_size=2,
                legend_loc='right margin',
                title=f'{prefix} ({velocity_mode}) — Velocity Arrows by {title}',
                legend_fontsize=6 if color == 'sample' else 8,
                ax=ax, show=False,
            )
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            _save(f'{output_dir}/{prefix}_arrows_{suffix}.pdf')
        except Exception as e:
            print(f"Error arrows by {color}: {e}")
            plt.close()

    # ---- CONFIDENCE ----
    print("\n--- Velocity confidence ---")
    scv.tl.velocity_confidence(adata)
    try:
        fig, ax = plt.subplots(figsize=(10, 8))
        conf = adata.obs['velocity_confidence'].values
        sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=conf, cmap='coolwarm',
                             s=5, alpha=0.8,
                             vmin=np.percentile(conf, 5),
                             vmax=np.percentile(conf, 95))
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_xlabel('UMAP_1'); ax.set_ylabel('UMAP_2')
        ax.set_title(f'{prefix} ({velocity_mode}) — Velocity Confidence')
        plt.colorbar(sc_plot, ax=ax, label='velocity_confidence')
        _save(f'{output_dir}/{prefix}_confidence.pdf')
    except Exception as e:
        print(f"Error confidence: {e}")
        plt.close()

    # ---- PSEUDOTIME ----
    print("\n--- Velocity pseudotime ---")
    scv.tl.velocity_pseudotime(adata)
    try:
        fig, ax = plt.subplots(figsize=(10, 8))
        pt = adata.obs['velocity_pseudotime'].values
        sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=pt, cmap='gnuplot', s=5, alpha=0.8)
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_xlabel('UMAP_1'); ax.set_ylabel('UMAP_2')
        ax.set_title(f'{prefix} ({velocity_mode}) — Velocity Pseudotime')
        plt.colorbar(sc_plot, ax=ax, label='velocity_pseudotime')
        _save(f'{output_dir}/{prefix}_pseudotime.pdf')
    except Exception as e:
        print(f"Error pseudotime: {e}")
        plt.close()

    # ---- TOP VELOCITY GENES ----
    print("\n--- Top velocity genes ---")
    if len(adata.obs[cell_type_col].unique()) > 1:
        try:
            scv.tl.rank_velocity_genes(adata, groupby=cell_type_col, min_corr=0.3)
            df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
            print(f"\nTop velocity genes by {cell_type_col} (head):")
            print(df.head(10))
            df.to_csv(f'{output_dir}/{prefix}_top_velocity_genes.csv')
            print(f"Saved: {output_dir}/{prefix}_top_velocity_genes.csv")
        except Exception as e:
            print(f"Error ranking velocity genes: {e}")


# =============================================================================
# MAIN PIPELINE
# =============================================================================
def run_velocity_full_atlas(adata_all, metadata_file, cell_type_col, output_dir):
    """Run scVelo on the full atlas (stochastic mode)."""
    print(f"\n{'='*60}")
    print("FULL ATLAS — scVelo (stochastic)")
    print(f"{'='*60}\n")

    adata = attach_seurat_metadata(adata_all, metadata_file)
    diagnostic_check_counts(adata)

    print("Preprocessing (filter_and_normalize, moments)...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print("Computing RNA velocity (stochastic)...")
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)

    adata = filter_invalid_velocity_cells(adata)
    if len(adata) == 0:
        raise RuntimeError("No valid cells remain after velocity filtering.")

    make_plots(adata, cell_type_col, output_dir,
               prefix='full_atlas', velocity_mode='stochastic')

    out_h5ad = f'{output_dir}/full_atlas_velocity.h5ad'
    adata.write(out_h5ad)
    print(f"\nSaved: {out_h5ad}")
    return adata


# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == '__main__':
    print("scVelo trajectory analysis — FULL FETAL BRAIN ATLAS")
    print(f"Loom dir: {LOOM_DIR}")
    print(f"Output:   {OUTPUT_DIR}\n")

    adata_all = load_all_loom_files(LOOM_DIR)

    print("\n" + "=" * 60)
    print("LOOMS LOADED — RUNNING FULL ATLAS VELOCITY")
    print("=" * 60)

    run_velocity_full_atlas(
        adata_all=adata_all,
        metadata_file=METADATA_FILE,
        cell_type_col=CELL_TYPE_COL,
        output_dir=OUTPUT_DIR,
    )

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
