"""
plot_mg_specific_genes.py

Standalone script to generate publication-quality pseudotime expression trend
plots for specific genes of interest (IL1B, CD86) for MG0 and MG1 clusters,
using the microglia CellRank h5ad output.

Reuses _plot_cluster_driver_trends_direct from run_cellrank_brain_5.py so the
figures match the style of the existing cluster driver plots exactly.

USAGE:
    python plot_mg_specific_genes.py

Edit H5AD_FILE, OUTPUT_DIR, CELL_TYPE_COL, and GENES_OF_INTEREST as needed.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr, rankdata, binned_statistic as bstat
from matplotlib.lines import Line2D
import warnings

# =============================================================================
# CONFIG
# =============================================================================
H5AD_FILE    = 'output_cellrank/mg_lineage/mg_lineage_cellrank.h5ad'  # EDIT
OUTPUT_DIR   = 'output_cellrank/mg_lineage'
CELL_TYPE_COL = 'cell_type_v8'

# Clusters to generate figures for
CLUSTERS_OF_INTEREST = ['MG0', 'MG1']

# Genes to plot for each cluster
GENES_OF_INTEREST = ['IL1B', 'CD86']

os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# PUBLICATION PLOT FUNCTION
# (copied from run_cellrank_brain_5.py to keep this script standalone)
# =============================================================================
def plot_specific_genes(adata, X_full, genes_full, genes_to_plot,
                        membership, cluster_name, output_dir,
                        subset_name, driver_stats=None, n_bins=30):
    """
    Publication-quality pseudotime trend plots for specific genes of interest.

    Matches the style of _plot_cluster_driver_trends_direct:
    - 300 dpi PDF output
    - In-cluster cells red, outside blue
    - Weighted trend line (black)
    - Spearman rho + q-value annotated per panel if driver_stats provided
    - Spines cleaned, fonts scaled for journal column width
    """
    # ---- Style constants (match run_cellrank_brain_5.py exactly) ----
    LABEL_FS    = 11
    TICK_FS     = 10
    TITLE_FS    = 12
    ANNOT_FS    = 9
    SUPTITLE_FS = 13
    POINT_SIZE  = 8
    POINT_ALPHA = 0.45
    TREND_LW    = 2.0
    DPI         = 300
    COLOR_IN    = '#E64B35'   # in-cluster (vermillion)
    COLOR_OUT   = '#4DBBD5'   # outside (steel blue)

    pt = adata.obs['velocity_pseudotime'].values.astype(float)
    valid = np.isfinite(pt)
    if valid.sum() < 10:
        print(f"  Too few cells with valid pseudotime for {cluster_name}, skipping.")
        return

    # Filter to genes present in var_names
    available_idx = [genes_full.index(g) for g in genes_to_plot
                     if g in genes_full]
    available = [genes_full[i] for i in available_idx]
    missing = [g for g in genes_to_plot if g not in genes_full]
    if missing:
        print(f"  WARNING: genes not in var_names (skipped): {missing}")
    if not available:
        print(f"  No requested genes available in {cluster_name}, skipping.")
        return

    pt_v    = pt[valid]
    X_v     = X_full[np.ix_(valid, available_idx)]
    mem_v   = membership[valid]

    bins        = np.linspace(pt_v.min(), pt_v.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    ncols = min(3, len(available))
    nrows = (len(available) + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.5 * ncols, 3.8 * nrows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    n_scatter   = min(500, len(pt_v))
    scatter_idx = np.random.choice(len(pt_v), n_scatter, replace=False)
    point_colors = np.where(mem_v[scatter_idx] > 0.5, COLOR_IN, COLOR_OUT)

    for gi, (gene, col_idx) in enumerate(zip(available, range(len(available)))):
        ax   = axes_flat[gi]
        expr = X_v[:, col_idx]

        # Weighted trend (weight by membership — emphasises in-cluster cells)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            wsum, _, _ = bstat(pt_v, expr * mem_v, statistic='sum', bins=bins)
            wbin, _, _ = bstat(pt_v, mem_v,        statistic='sum', bins=bins)
        with np.errstate(invalid='ignore', divide='ignore'):
            trend = np.where(wbin > 0, wsum / wbin, np.nan)

        ax.scatter(pt_v[scatter_idx], expr[scatter_idx],
                   c=point_colors, s=POINT_SIZE, alpha=POINT_ALPHA,
                   linewidths=0, rasterized=True)

        vb = np.isfinite(trend)
        if vb.sum() >= 3:
            ax.plot(bin_centers[vb], trend[vb],
                    color='black', lw=TREND_LW, zorder=5)

        ax.set_title(gene, fontsize=TITLE_FS, fontweight='bold', pad=4)

        row_i = gi // ncols
        col_i = gi % ncols
        ax.set_ylabel('Normalized expression', fontsize=LABEL_FS) \
            if col_i == 0 else ax.set_ylabel('')
        ax.set_xlabel('Pseudotime', fontsize=LABEL_FS) \
            if (row_i == nrows - 1 or gi >= len(available) - ncols) \
            else ax.set_xlabel('')

        ax.tick_params(labelsize=TICK_FS)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Annotate with Spearman rho + q-value if stats available
        if driver_stats is not None and gene in driver_stats.index:
            rho  = driver_stats.loc[gene, 'correlation']
            qval = driver_stats.loc[gene, 'qvalue']
            q_str = ('q < 0.001' if qval < 0.001
                     else f'q = {qval:.3f}' if qval < 0.01
                     else f'q = {qval:.2f}')
            ax.text(0.97, 0.97, f'ρ = {rho:.2f}, {q_str}',
                    transform=ax.transAxes,
                    ha='right', va='top', fontsize=ANNOT_FS, color='#333333',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              alpha=0.7, edgecolor='none'))

    for gi in range(len(available), len(axes_flat)):
        axes_flat[gi].set_visible(False)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_IN, markersize=7,
               label=cluster_name),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_OUT, markersize=7,
               label='Other clusters'),
        Line2D([0], [0], color='black', lw=2, label='Weighted trend'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=ANNOT_FS + 1, frameon=False,
               bbox_to_anchor=(0.5, -0.02))

    readable_subset = subset_name.replace('_', ' ').title()
    fig.suptitle(
        f'{readable_subset} — {", ".join(genes_to_plot)} in {cluster_name}',
        fontsize=SUPTITLE_FS, fontweight='bold', y=1.01,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])

    safe = cluster_name.replace('/', '_').replace(' ', '_')
    gene_str = '_'.join(g.replace('/', '_') for g in genes_to_plot)
    path = f'{output_dir}/{subset_name}_{safe}_{gene_str}.pdf'
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


# =============================================================================
# COMPUTE SPEARMAN STATS FOR ANNOTATION
# =============================================================================
def compute_gene_stats(X_full, genes_full, genes_to_plot, membership):
    """
    Compute Spearman correlation + BH-corrected q-value for each gene vs
    cluster membership, for use in panel annotations.
    """
    results = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for gene in genes_to_plot:
            if gene not in genes_full:
                continue
            idx  = genes_full.index(gene)
            expr = X_full[:, idx]
            if expr.std() < 1e-8:
                continue
            rho, pval = spearmanr(expr, membership)
            results.append({'gene': gene, 'correlation': rho, 'pvalue': float(pval)})

    if not results:
        return None

    df = pd.DataFrame(results).set_index('gene')
    n     = len(df)
    ranks = rankdata(df['pvalue'].values)
    df['qvalue'] = np.minimum(1.0, df['pvalue'].values * n / ranks)
    return df


# =============================================================================
# LOOM GENE INJECTION
# Loads expression for genes missing from the scVelo HVG-filtered h5ad
# directly from the original loom files and adds them to adata.
# =============================================================================
LOOM_DIR = '/gscratch/kawaldorflab/jcorn427/fet_brain/loom_files'  # EDIT

def inject_genes_from_loom(adata, genes_to_inject, loom_dir, cell_type_col):
    """
    Load spliced counts for specific genes from loom files and add them
    to adata so they can be plotted alongside HVG-filtered genes.

    Uses the barcode index already in adata.obs_names to match cells.
    Injected genes are added as new columns in adata.X (dense slice only —
    they are not added to adata.var permanently to avoid breaking sparse ops).
    Returns a supplemental dict: {gene: expr_array} for injected genes.
    """
    import glob
    import scipy.sparse as sp2

    loom_files = sorted(glob.glob(f'{loom_dir}/*.loom'))
    if not loom_files:
        raise FileNotFoundError(f"No loom files found in {loom_dir}")

    # Build barcode → loom cell index lookup across all loom files
    print(f"  Loading loom files to inject: {genes_to_inject}")
    gene_exprs = {g: np.zeros(adata.shape[0], dtype=np.float32)
                  for g in genes_to_inject}
    matched_total = 0

    for loom_path in loom_files:
        import scanpy as sc2
        loom_adata = sc2.read_loom(loom_path)
        loom_adata.var_names_make_unique()

        # Check which requested genes are present
        available = [g for g in genes_to_inject if g in loom_adata.var_names]
        if not available:
            continue

        # Clean loom barcodes to match adata.obs_names format
        # loom: "SAMPLE:BARCODEx" → "SAMPLE_BARCODE-1"
        sample_name = os.path.basename(loom_path).replace('.loom', '')
        def clean_bc(bc):
            parts = bc.split(':')
            if len(parts) == 2:
                return f"{sample_name}_{parts[1].rstrip('x')}-1"
            return bc

        loom_bcs = [clean_bc(bc) for bc in loom_adata.obs_names]

        # Also try lowercase sample matching (CTRL13 vs Ctrl13)
        loom_bc_lower = [bc.lower() for bc in loom_bcs]
        adata_bc_lower = [bc.lower() for bc in adata.obs_names]

        bc_to_adata_idx = {bc: i for i, bc in enumerate(adata_bc_lower)}

        # Get spliced expression for available genes
        spliced = loom_adata.layers['spliced']
        if sp2.issparse(spliced):
            spliced = spliced.toarray()

        n_matched = 0
        for loom_i, loom_bc in enumerate(loom_bc_lower):
            if loom_bc in bc_to_adata_idx:
                adata_i = bc_to_adata_idx[loom_bc]
                for gene in available:
                    gene_i = list(loom_adata.var_names).index(gene)
                    gene_exprs[gene][adata_i] = spliced[loom_i, gene_i]
                n_matched += 1

        matched_total += n_matched
        print(f"    {sample_name}: matched {n_matched} cells, "
              f"genes: {available}")

    print(f"  Total cells matched: {matched_total} / {adata.shape[0]}")

    # Log-normalize injected counts to match scVelo-processed expression
    for gene in gene_exprs:
        raw = gene_exprs[gene]
        # Simple library-size normalization + log1p to match scVelo output
        cell_totals = np.array(adata.X.sum(axis=1)).flatten() \
            if sp2.issparse(adata.X) else adata.X.sum(axis=1)
        median_total = np.median(cell_totals[cell_totals > 0])
        with np.errstate(divide='ignore', invalid='ignore'):
            gene_exprs[gene] = np.log1p(raw / (cell_totals + 1e-8) * median_total)
        print(f"  {gene}: max={gene_exprs[gene].max():.3f}, "
              f"nonzero={np.sum(gene_exprs[gene] > 0)} cells")

    return gene_exprs



if __name__ == '__main__':
    print(f"Loading: {H5AD_FILE}")
    if not os.path.exists(H5AD_FILE):
        raise FileNotFoundError(
            f"{H5AD_FILE} not found.\n"
            f"Make sure run_cellrank_brain_5.py has completed for mg_lineage."
        )

    adata = sc.read_h5ad(H5AD_FILE)
    print(f"Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
    print(f"Cell types present: {sorted(adata.obs[CELL_TYPE_COL].unique())}")

    if 'velocity_pseudotime' not in adata.obs.columns:
        raise RuntimeError(
            "velocity_pseudotime not found in adata.obs. "
            "Make sure the h5ad was saved after scVelo was run."
        )

    # Get dense expression matrix for HVG genes
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X     = np.array(X, dtype=float)
    genes = list(adata.var_names)

    # Check which requested genes are missing from the HVG set
    missing_genes = [g for g in GENES_OF_INTEREST if g not in genes]

    # Inject missing genes from loom files
    injected = {}
    if missing_genes:
        print(f"\nGenes missing from HVG set: {missing_genes}")
        print("Injecting from loom files...")
        try:
            injected = inject_genes_from_loom(
                adata, missing_genes, LOOM_DIR, CELL_TYPE_COL
            )
        except Exception as e:
            print(f"WARNING: loom injection failed: {e}")
            print("Missing genes will be skipped.")

    # Build augmented expression matrix including injected genes
    if injected:
        injected_array = np.stack([injected[g] for g in missing_genes], axis=1)
        X_aug     = np.hstack([X, injected_array])
        genes_aug = genes + missing_genes
        print(f"\nAugmented matrix: {X_aug.shape[1]} genes "
              f"({len(genes)} HVG + {len(missing_genes)} injected)")
    else:
        X_aug     = X
        genes_aug = genes

    for cluster in CLUSTERS_OF_INTEREST:
        if cluster not in adata.obs[CELL_TYPE_COL].values:
            print(f"WARNING: cluster '{cluster}' not found, skipping.")
            continue

        print(f"\n--- {cluster} ---")
        membership = (adata.obs[CELL_TYPE_COL] == cluster).astype(float).values
        n_in = int(membership.sum())
        print(f"  Cells in cluster: {n_in} / {len(membership)}")

        driver_stats = compute_gene_stats(
            X_aug, genes_aug, GENES_OF_INTEREST, membership
        )
        if driver_stats is not None:
            print(f"  Gene stats:")
            print(driver_stats[["correlation", "qvalue"]].to_string())

        plot_specific_genes(
            adata         = adata,
            X_full        = X_aug,
            genes_full    = genes_aug,
            genes_to_plot = GENES_OF_INTEREST,
            membership    = membership,
            cluster_name  = cluster,
            output_dir    = OUTPUT_DIR,
            subset_name   = "mg_lineage",
            driver_stats  = driver_stats,
        )

    print("\nDone.")
