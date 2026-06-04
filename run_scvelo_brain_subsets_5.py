"""
scVelo trajectory analysis on LINEAGE SUBSETS of the fetal brain atlas.

Each subset (astrocyte lineage, OL lineage, etc.) gets its own scVelo run
on a re-embedded UMAP, with stochastic-mode velocity. Output h5ad checkpoints
feed into the CellRank script.

Run AFTER exporting per-subset Seurat metadata to CSVs with columns:
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
# Identical to the full atlas script — same cell type always gets same color.
# =============================================================================
CELL_TYPE_PALETTE = {
    # --- Radial glia ---
    'vRG':     '#145A32',   # dark forest green (goes with OL lineage)
    'oRG-G':   '#1F618D',   # steel blue
    # --- Progenitor / intermediate states ---
    'APC':     '#D6EAF8',   # very pale blue (matches astrocyte lineage)
    'UP/I':    '#F0B27A',
    'IPC':     '#784212',
    # --- Astrocyte lineage (light → dark blue) ---
    'AS0':     '#AED6F1',
    'AS1':     '#5DADE2',
    'AS2':     '#2874A6',
    # --- Oligodendrocyte lineage (light → dark green) ---
    'OPC':     '#A9DFBF',
    'COP1':    '#52BE80',
    'COP2':    '#27AE60',
    'OL1':     '#1E8449',
    'OL2':     '#0B5345',
    # --- Excitatory neurons (purples) ---
    'ExN':     '#7D3C98',
    'ImmN':    '#4A235A',
    # --- CGE interneurons ---
    'CGE':     '#F1948A',
    # --- MGE interneurons ---
    'MGE':     '#D81B60',
    'MGE-SST': '#880E4F',
    # --- Microglia ---
    'MG0':     '#B2BABB',
    'MG1':     '#5D6D7E',
    # --- Other immune / myeloid ---
    'BAM':     '#1ABC9C',
    # --- Non-neural populations ---
    'ChP':     '#F4D03F',
    'Ep':      '#F0A500',
    'EC':      '#FF7043',
    # --- Stromal populations (same warm sandy family) ---
    'FIB':     '#D5C4A1',   # pale sandy beige
    'PC':      '#9E8B6A',   # medium sand
    'MPC':     '#6D5D44',   # dark khaki/umber
}


def get_palette_for_adata(adata, cell_type_col):
    """Return ordered hex colors for adata.obs[cell_type_col].cat.categories.
    Uses lineage-based CELL_TYPE_PALETTE — for full atlas plots only."""
    if not hasattr(adata.obs[cell_type_col], 'cat'):
        adata.obs[cell_type_col] = adata.obs[cell_type_col].astype('category')
    cats = list(adata.obs[cell_type_col].cat.categories)
    return [CELL_TYPE_PALETTE.get(c, '#999999') for c in cats]


# High-contrast palette for subset trajectory plots.
# Identical to SUBSET_PALETTE in atlas_palette.R — assigns maximally distinct
# hues positionally so nearby clusters are always clearly different regardless
# of biological lineage. Used in make_publication_celltype_stream for subsets.
SUBSET_PALETTE = [
    '#E41A1C',  # red
    '#377EB8',  # blue
    '#4DAF4A',  # green
    '#FF7F00',  # orange
    '#984EA3',  # purple
    '#A65628',  # brown
    '#F781BF',  # pink
    '#1ABC9C',  # teal
    '#D4AC0D',  # gold
    '#2C3E50',  # dark slate
]


def get_subset_palette_for_adata(adata, cell_type_col):
    """
    Return high-contrast colors for each cell type category.

    If a name-based color map was built from category_order (stored in
    adata.uns['{col}_color_map']), uses that — colors are looked up by
    cell type NAME so they survive category reordering during preprocessing.

    Falls back to positional assignment if no map is stored.
    """
    if not hasattr(adata.obs[cell_type_col], 'cat'):
        adata.obs[cell_type_col] = adata.obs[cell_type_col].astype('category')
    cats = list(adata.obs[cell_type_col].cat.categories)

    color_map_key = f'{cell_type_col}_color_map'
    if color_map_key in adata.uns:
        color_map = adata.uns[color_map_key]
        colors = [color_map.get(c, '#999999') for c in cats]
        print(f"  Color assignment (by name):")
        for cat, col in zip(cats, colors):
            print(f"    {cat} → {col}")
        return colors
    else:
        # Fallback: positional assignment
        n = len(SUBSET_PALETTE)
        return [SUBSET_PALETTE[i % n] for i in range(len(cats))]


# =============================================================================
# CONFIG
# =============================================================================
# Directory containing per-sample .loom files produced by velocyto run10x.
# Each file should be named {SAMPLE_ID}.loom (e.g. CTRL13.loom, SAL7.loom).
LOOM_DIR = '/path/to/loom_files'

OUTPUT_DIR = 'output_subsets'
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# HELPERS  (mirror the full-atlas script for consistency)
# =============================================================================
def extract_sample_id(filename):
    """Extract sample ID from CTRL/SAL loom filenames."""
    basename = os.path.basename(filename).replace('.loom', '')
    match = re.match(r'^(CTRL|SAL)(\d+)$', basename)
    if match:
        return f"{match.group(1)}{match.group(2)}"
    return None


def load_all_loom_files(loom_dir):
    """Load and concatenate all loom files."""
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
            print(f"WARNING: could not parse {sample_name}, skipping")
            continue

        print(f"Loading {sample_name} -> {sample_id}...")
        adata_loom = sc.read_loom(loom_file)
        adata_loom.obs['loom_sample'] = sample_name
        adata_loom.obs['sample_id'] = sample_id
        adata_loom.var_names_make_unique()
        adatas.append(adata_loom)

    print("\nConcatenating...")
    adata_all = adatas[0] if len(adatas) == 1 else sc.concat(adatas, axis=0, join='outer')
    print(f"Combined: {adata_all.shape[0]} cells, {adata_all.shape[1]} genes")

    def extract_barcode(full_name):
        parts = full_name.split(':')
        if len(parts) == 2:
            return parts[1].rstrip('x') + '-1'
        return full_name

    adata_all.obs['barcode_clean'] = [extract_barcode(x) for x in adata_all.obs_names]
    print(f"\nBarcode example: {adata_all.obs_names[0]} -> {adata_all.obs['barcode_clean'].iloc[0]}")
    return adata_all


def attach_seurat_metadata(adata_all, metadata_file):
    """Match loom cells against Seurat metadata via {sample_id}_{barcode} keys."""
    metadata = pd.read_csv(metadata_file, index_col='barcode')
    print(f"Loaded metadata: {metadata.shape[0]} cells")

    metadata['barcode_clean'] = metadata.index.str.split('_').str[-1]
    metadata['sample_id'] = metadata['sample']

    # Normalize case to lowercase on both sides — loom filenames are all-caps
    # (CTRL13, SAL7) but Seurat stores sample IDs in title case (Ctrl13, Sal7).
    metadata['match_key'] = (metadata['sample_id'].str.lower()
                             + '_' + metadata['barcode_clean'])
    adata_all.obs['match_key'] = (adata_all.obs['sample_id'].str.lower()
                                  + '_' + adata_all.obs['barcode_clean'])

    common_keys = set(metadata['match_key']).intersection(set(adata_all.obs['match_key']))
    print(f"Matching cells: {len(common_keys)}")

    metadata_matched = metadata[metadata['match_key'].isin(common_keys)].copy()
    match_key_to_obs = dict(zip(adata_all.obs['match_key'], adata_all.obs_names))
    adata_obs_names = [match_key_to_obs[key] for key in metadata_matched['match_key']]

    adata = adata_all[adata_obs_names, :].copy()
    pct = 100 * len(metadata_matched) / len(metadata)
    print(f"Matched {len(metadata_matched)} / {len(metadata)} cells ({pct:.1f}%)")

    if len(metadata_matched) == 0:
        raise RuntimeError("No cells matched.")

    adata.obs_names = metadata_matched.index.values

    for col in metadata_matched.columns:
        if col not in ('barcode_clean', 'sample_id', 'match_key'):
            adata.obs[col] = metadata_matched[col].values

    if 'umap_1' in metadata_matched.columns and 'umap_2' in metadata_matched.columns:
        adata.obsm['X_umap'] = metadata_matched[['umap_1', 'umap_2']].values
        print("Added UMAP coordinates from Seurat")
    else:
        raise RuntimeError("Metadata missing umap_1/umap_2 columns.")

    return adata


def diagnostic_check_counts(adata):
    """Sanity-check that spliced/unspliced layers are raw integer counts."""
    spliced = adata.layers['spliced']
    unspliced = adata.layers['unspliced']
    s_vals = spliced.data if sp.issparse(spliced) else spliced.flatten()
    u_vals = unspliced.data if sp.issparse(unspliced) else unspliced.flatten()
    total = adata.shape[0] * adata.shape[1]

    print("\n--- PRE-NORMALIZATION DIAGNOSTIC ---")
    print(f"Cells: {adata.shape[0]}, Genes: {adata.shape[1]}")
    print(f"Spliced   nonzero mean: {s_vals.mean():.2f}  max: {s_vals.max():.0f}  "
          f"% nonzero: {100*len(s_vals)/total:.2f}%  int: {np.all(s_vals == s_vals.astype(int))}")
    print(f"Unspliced nonzero mean: {u_vals.mean():.2f}  max: {u_vals.max():.0f}  "
          f"% nonzero: {100*len(u_vals)/total:.2f}%  int: {np.all(u_vals == u_vals.astype(int))}")
    print("--- END DIAGNOSTIC ---\n")


def filter_invalid_velocity_cells(adata):
    """Remove cells with no finite velocity values (see full-atlas script for rationale)."""
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
          f"mean: {frac_finite.mean():.3f}  median: {np.median(frac_finite):.3f}  "
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


def _save(fig_path):
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_path}")


def make_publication_celltype_stream(adata, cell_type_col, output_dir,
                                     subset_name, velocity_mode='stochastic'):
    """
    Publication-quality velocity stream plot colored by cell type.
    See full atlas script for full documentation.
    """
    print("\n--- Publication stream plot (cell type) ---")

    umap = adata.obsm['X_umap']
    xlim = (umap[:, 0].min() - 1, umap[:, 0].max() + 1)
    ylim = (umap[:, 1].min() - 1, umap[:, 1].max() + 1)

    # Use high-contrast positional palette for subsets — lineage-based colors
    # are too similar on smaller plots (e.g. all blues in astrocyte lineage).
    colors = get_subset_palette_for_adata(adata, cell_type_col)
    adata.uns[f'{cell_type_col}_colors'] = colors

    n_cats = len(adata.obs[cell_type_col].cat.categories)
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
                size=40,
                alpha=0.85,
                linewidth=0.8,
                density=1.0,
                arrow_size=1.2,
                title='',
                ax=ax,
                show=False,
            )

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel('UMAP 1', fontsize=12, labelpad=6)
            ax.set_ylabel('UMAP 2', fontsize=12, labelpad=6)
            ax.tick_params(labelsize=10)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            out_path = f'{output_dir}/{subset_name}_stream_celltype_pub.{fmt}'
            fig.savefig(out_path, dpi=300, bbox_inches='tight', format=fmt)
            plt.close(fig)
            print(f"Saved: {out_path}")

        except Exception as e:
            print(f"  Publication stream ({fmt}) failed: {e}")
            plt.close()


def make_plots(adata, cell_type_col, output_dir, subset_name,
               velocity_mode='stochastic'):
    """Generate stream, arrow, confidence, and pseudotime plots for a subset."""
    print("\nGenerating plots...")

    umap = adata.obsm['X_umap']
    xlim = (umap[:, 0].min() - 1, umap[:, 0].max() + 1)
    ylim = (umap[:, 1].min() - 1, umap[:, 1].max() + 1)

    print("Recomputing velocity embedding for visualization...")
    scv.tl.velocity_embedding(adata, basis='umap')

    # ---- PUBLICATION STREAM PLOT (cell type, PDF + SVG) ----
    make_publication_celltype_stream(
        adata, cell_type_col, output_dir, subset_name, velocity_mode
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
                title=f'{subset_name} ({velocity_mode}) — Stream by {title}',
                legend_fontsize=6 if color == 'sample' else 8,
                ax=ax, show=False,
            )
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            _save(f'{output_dir}/{subset_name}_stream_{suffix}.svg')
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
                title=f'{subset_name} ({velocity_mode}) — Arrows by {title}',
                legend_fontsize=6 if color == 'sample' else 8,
                ax=ax, show=False,
            )
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            _save(f'{output_dir}/{subset_name}_arrows_{suffix}.pdf')
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
        ax.set_title(f'{subset_name} ({velocity_mode}) — Velocity Confidence')
        plt.colorbar(sc_plot, ax=ax, label='velocity_confidence')
        _save(f'{output_dir}/{subset_name}_confidence.pdf')
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
        ax.set_title(f'{subset_name} ({velocity_mode}) — Velocity Pseudotime')
        plt.colorbar(sc_plot, ax=ax, label='velocity_pseudotime')
        _save(f'{output_dir}/{subset_name}_pseudotime.pdf')
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
            df.to_csv(f'{output_dir}/{subset_name}_top_velocity_genes.csv')
            print(f"Saved: {output_dir}/{subset_name}_top_velocity_genes.csv")
        except Exception as e:
            print(f"Error ranking velocity genes: {e}")


# =============================================================================
# PIPELINE
# =============================================================================
def run_velocity_subset(adata_all, metadata_file, subset_name, cell_type_col,
                        output_dir, category_order=None):
    """
    Run scVelo (stochastic mode) on a subset.

    Parameters
    ----------
    adata_all : AnnData
        Pre-loaded concatenated loom data.
    metadata_file : str
        Path to subset's CSV (barcode-indexed, with sample/cell_type/UMAP).
    subset_name : str
        Used for output filename prefixes (e.g. 'astro_lineage', 'OL_lineage').
    cell_type_col : str
        Annotation column to color plots and rank velocity genes by.
    output_dir : str
        Where to write outputs.
    category_order : list, optional
        Explicit category order for cell_type_col. Must match the factor level
        order in your Seurat object so the subset palette assigns the same
        color to each cell type in both the R UMAP and the Python stream plot.
        Get the correct order with: levels(your_subset$cell_type_col) in R.
        If None, uses the default alphabetical order from pandas Categorical.
    """
    print(f"\n{'='*60}")
    print(f"SUBSET: {subset_name} — scVelo (stochastic)")
    print(f"{'='*60}\n")

    adata = attach_seurat_metadata(adata_all, metadata_file)

    # Apply category order so palette colors match the R UMAP exactly.
    # The subset palette assigns colors positionally — category order must
    # be identical between R and Python for colors to be consistent.
    import pandas as pd
    if not hasattr(adata.obs[cell_type_col], 'cat') or \
       not hasattr(adata.obs[cell_type_col].dtype, 'categories'):
        adata.obs[cell_type_col] = adata.obs[cell_type_col].astype('category')

    if category_order is not None:
        current = list(adata.obs[cell_type_col].cat.categories)
        missing = [c for c in category_order if c not in current]
        extra   = [c for c in current if c not in category_order]
        if missing:
            print(f"  WARNING: category_order contains unknown categories: {missing}")
        full_order = [c for c in category_order if c in current] + extra
        adata.obs[cell_type_col] = adata.obs[cell_type_col].cat.reorder_categories(full_order)
        print(f"  Category order set to: {full_order}")
    else:
        full_order = sorted(adata.obs[cell_type_col].cat.categories.tolist())
        print(f"  Category order (alphabetical default): {full_order}")

    # Build name→color map NOW before preprocessing can reset category order.
    # get_subset_palette_for_adata will look up colors by name from this map,
    # so the correct color follows each cell type regardless of later reordering.
    n = len(SUBSET_PALETTE)
    color_map = {cat: SUBSET_PALETTE[i % n] for i, cat in enumerate(full_order)}
    adata.uns[f'{cell_type_col}_color_map'] = color_map
    print(f"  Color map stored: { {k: v for k, v in color_map.items()} }")

    diagnostic_check_counts(adata)

    print("Preprocessing (filter_and_normalize, moments)...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print("Computing RNA velocity (stochastic)...")
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)

    # Snapshot the color map before filtering (copy may drop uns entries)
    color_map_key = f'{cell_type_col}_color_map'
    saved_color_map = adata.uns.get(color_map_key, None)

    adata = filter_invalid_velocity_cells(adata)
    if len(adata) == 0:
        raise RuntimeError(f"No valid cells remain for {subset_name}.")

    # Restore color map if filtering dropped it
    if saved_color_map is not None and color_map_key not in adata.uns:
        adata.uns[color_map_key] = saved_color_map
        print(f"  Color map restored after cell filtering.")

    make_plots(adata, cell_type_col, output_dir, subset_name=subset_name,
               velocity_mode='stochastic')

    out_h5ad = f'{output_dir}/{subset_name}_velocity.h5ad'
    adata.write(out_h5ad)
    print(f"\nSaved: {out_h5ad}")
    return adata


# =============================================================================
# MAIN — define your subsets here
# =============================================================================
if __name__ == '__main__':
    print("scVelo trajectory analysis — FETAL BRAIN SUBSETS")
    print(f"Loom dir: {LOOM_DIR}")
    print(f"Output:   {OUTPUT_DIR}\n")

    adata_all = load_all_loom_files(LOOM_DIR)

    print("\n" + "=" * 60)
    print("LOOMS LOADED — PROCESSING SUBSETS")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # Define subsets here. Add/remove/comment as needed.
    # Each subset = one Seurat metadata CSV exported separately.
    # -------------------------------------------------------------------------
    subsets = [
        # Astrocyte lineage: RG -> APC -> immature astrocyte -> astrocyte
        dict(
            metadata_file='astro_metadata.csv',
            subset_name='astro_lineage',
            cell_type_col='cell_type_v8',
        ),
        # Oligodendrocyte lineage: RG -> pre-OPC -> OPC -> pre-OL -> OL
        dict(
            metadata_file='oligo_metadata.csv',
            subset_name='ol_lineage',
            cell_type_col='cell_type_v8',
        ),
        # Neuronal lineage: RG -> IPC -> neuron
        # category_order = present cell types in the order they appear in the
        # full atlas factor levels: levels(neuro_sub$cell_type_v8)
        # CGE(6), ExN(12), ImmN(14), IPC(15), MGE(18), MGE-SST(19)
        dict(
            metadata_file='exn_metadata.csv',
            subset_name='exn_lineage',
            cell_type_col='cell_type_v8',
            category_order=['CGE', 'ExN', 'ImmN', 'IPC', 'MGE', 'MGE-SST'],
        ),

        dict(
            metadata_file='mg_metadata.csv',
            subset_name='mg_lineage',
            cell_type_col='cell_type_v8',
        ),

        dict(
            metadata_file='neuro_metadata.csv',
            subset_name='neuro_lineage',
            cell_type_col='cell_type_v8',
        ),
    ]

    for subset in subsets:
        try:
            run_velocity_subset(
                adata_all=adata_all,
                output_dir=OUTPUT_DIR,
                **subset,
            )
        except Exception as e:
            print(f"\nFAILED for {subset['subset_name']}: {e}")
            print("Continuing to next subset...\n")

    print("\n" + "=" * 60)
    print("ALL SUBSETS DONE")
    print("=" * 60)
