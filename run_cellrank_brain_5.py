"""
CellRank 2 fate analysis on scVelo h5ad outputs.

Reads h5ad checkpoints written by run_scvelo_brain_full.py or
run_scvelo_brain_subsets.py and runs:

    1. Velocity-directed PAGA (cluster-level lineage graph)
    2. Terminal & initial state identification (GPCCA estimator)
    3. Fate probabilities (per-cell probability of reaching each terminal state)
    4. Lineage drivers (genes correlated with each fate)
    5. Gene trends along pseudotime (per lineage)

Designed for CellRank 2.x.

USAGE:
    Configure ANALYSES below (one entry per h5ad you want to process).
    Each analysis points to one h5ad, sets terminal states (auto or manual),
    and chooses which marker genes to track for gene trends.

    For the astrocyte lineage, terminal_states='auto' usually finds the right
    endpoint, but you can pin it to specific cell type names if needed.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import scvelo as scv

# NOTE on petsc4py + HPC clusters:
# petsc4py from conda-forge is compiled against its own bundled MPI which
# conflicts with the system MPI on Slurm nodes, causing MPI_Abort if init()
# is called. Do NOT call petsc4py.init() here.
# CellRank will detect petsc4py automatically if the environment is compatible.
# If not, the lgmres/bicgstab scipy cascade in run_fate_probabilities handles it.
import cellrank as cr

# CellRank 2 kernel + estimator imports
from cellrank.kernels import VelocityKernel, ConnectivityKernel
from cellrank.estimators import GPCCA

warnings.filterwarnings('ignore', category=FutureWarning)

scv.settings.verbosity = 3
sc.settings.verbosity = 3
sc.set_figure_params(dpi=100, frameon=False, figsize=(8, 6))


# =============================================================================
# CONFIG — one entry per subset you want to analyze
# =============================================================================
ANALYSES = [
    dict(
        h5ad_file='output_subsets/astro_lineage_velocity.h5ad',
        subset_name='astro_lineage',
        cell_type_col='cell_type_v8',
        terminal_states='auto',
        initial_states='auto',
        gene_trend_markers=[
            'FABP7', 'HOPX', 'VIM', 'NDRG2', 'MT3',
            'AQP4', 'GFAP', 'GJA1', 'SLC1A2', 'SLC1A3',
            'AGT', 'MLC1', 'GLUL', 'SPARCL1', 'MFGE8',
        ],
        n_drivers=50,
        n_macrostates="auto",
    ),
    dict(
        h5ad_file='output_subsets/ol_lineage_velocity.h5ad',
        subset_name='ol_lineage',
        cell_type_col='cell_type_v8',
        terminal_states='auto',
        initial_states='auto',
        gene_trend_markers=[
            'PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'SOX10',
            'BCAS1', 'GPR17', 'ENPP6',
            'MBP', 'PLP1', 'MOG', 'MAG', 'CLDN11', 'UGT8',
        ],
        n_drivers=50,
        n_macrostates='auto',
    ),
    dict(
        h5ad_file='output_subsets/exn_lineage_velocity.h5ad',
        subset_name='exn_lineage',
        cell_type_col='cell_type_v8',
        terminal_states='auto',
        initial_states='auto',
        gene_trend_markers=[
            'VIM', 'SOX2', 'HES1',          # progenitor
            'EOMES', 'PPP1R17',              # IPC
            'DCX', 'NEUROD2', 'NEUROD6',     # immature neuron
            'RBFOX3', 'SYP', 'GAD1', 'GAD2', # mature neuron
            'SLC17A7', 'TBR1',               # excitatory
        ],
        n_drivers=50,
        n_macrostates='auto',
    ),
    dict(
        h5ad_file='output_subsets/neuron_lineage_velocity.h5ad',
        subset_name='neuron_lineage',
        cell_type_col='cell_type_v8',
        terminal_states='auto',
        initial_states='auto',
        # Explicit order to match Seurat factor levels — fixes IPC/ImmN color swap.
        # Check your Seurat object with: levels(neuron_subset$cell_type_v8)
        # and paste that order here.
        category_order=['CGE', 'ExN', 'IPC', 'ImmN', 'MGE', 'MGE-IN'],
        gene_trend_markers=[
            'VIM', 'SOX2', 'HES1',           # radial glia / progenitor
            'EOMES', 'PPP1R17',               # IPC
            'DCX', 'NEUROD2', 'NEUROD6',      # immature neuron
            'RBFOX3', 'SYP',                  # mature neuron
            'SLC17A7', 'TBR1',               # excitatory
            'GAD1', 'GAD2', 'DLX2',          # inhibitory (CGE/MGE)
            'SST', 'LHX6',                    # MGE-specific
        ],
        n_drivers=50,
        n_macrostates='auto',
    ),
    dict(
        h5ad_file='output_subsets/mg_lineage_velocity.h5ad',
        subset_name='mg_lineage',
        cell_type_col='cell_type_v8',
        terminal_states=['MG0', 'MG1'],
        initial_states='auto',
        gene_trend_markers=[
            'CX3CR1', 'P2RY12', 'TMEM119',  # homeostatic MG
            'C1QB', 'TYROBP', 'TREM2',
            'AIF1', 'ITGAM', 'PTPRC',
            'MKI67', 'TOP2A',                # proliferating
        ],
        n_drivers=50,
        n_macrostates='auto',
    ),

    # -------- FULL ATLAS (optional; large data, slower) --------
    # dict(
    #     h5ad_file='output_full/full_atlas_velocity.h5ad',
    #     subset_name='full_atlas',
    #     cell_type_col='cell_type_v1',
    #     terminal_states='auto',
    #     initial_states='auto',
    #     gene_trend_markers=['FABP7', 'AQP4', 'GFAP', 'PDGFRA', 'MBP'],
    #     n_drivers=50,
    #     n_macrostates='auto',
    # ),
]

OUTPUT_DIR_BASE = 'output_cellrank'
os.makedirs(OUTPUT_DIR_BASE, exist_ok=True)


# =============================================================================
# CELLRANK PIPELINE
# =============================================================================
def build_kernel(adata, velocity_weight=0.8):
    """
    Build a combined velocity + connectivity kernel.

    Pure velocity kernels can be noisy in regions with weak splicing signal;
    mixing in 20% connectivity stabilizes the transition matrix without
    overriding velocity-derived directionality.
    """
    print(f"\nBuilding kernel ({velocity_weight:.0%} velocity + "
          f"{1-velocity_weight:.0%} connectivity)...")

    vk = VelocityKernel(adata)
    vk.compute_transition_matrix()

    ck = ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    combined = velocity_weight * vk + (1 - velocity_weight) * ck
    combined.compute_transition_matrix()

    return combined


def run_paga_velocity(adata, cell_type_col, output_dir, subset_name):
    """
    PAGA with velocity-directed edges.

    Uses scVelo's velocity-aware PAGA: edges are weighted by net velocity flow
    between clusters, giving directed lineage relationships at cluster level.
    """
    print("\n--- PAGA (velocity-directed) ---")
    try:
        # scvelo's paga uses the velocity_graph to direct edges.
        # NOTE: use_time_prior is omitted — it requires velocity_pseudotime to be
        # stored as a graph prior in a specific internal format that is only set up
        # in certain scVelo code paths. Omitting it is safe; PAGA still uses the
        # velocity_graph for edge directionality via vkey='velocity'.
        scv.tl.paga(
            adata,
            groups=cell_type_col,
            vkey='velocity',
        )

        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))
        scv.pl.paga(
            adata,
            basis='umap',
            size=50,
            alpha=0.5,
            min_edge_width=2,
            node_size_scale=1.5,
            title=f'{subset_name} — Velocity-directed PAGA',
            ax=ax,
            show=False,
        )
        plt.tight_layout()
        fig_path = f'{output_dir}/{subset_name}_paga.pdf'
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {fig_path}")

        # Save transition probabilities as CSV
        if 'paga' in adata.uns and 'transitions_confidence' in adata.uns['paga']:
            trans = adata.uns['paga']['transitions_confidence']
            cats = adata.obs[cell_type_col].cat.categories \
                if hasattr(adata.obs[cell_type_col], 'cat') \
                else sorted(adata.obs[cell_type_col].unique())
            trans_df = pd.DataFrame(trans.toarray() if hasattr(trans, 'toarray') else trans,
                                    index=cats, columns=cats)
            csv_path = f'{output_dir}/{subset_name}_paga_transitions.csv'
            trans_df.to_csv(csv_path)
            print(f"Saved: {csv_path}")

    except Exception as e:
        print(f"PAGA failed: {e}")


def run_terminal_states(g, adata, cell_type_col, output_dir, subset_name,
                        terminal_states, initial_states, n_macrostates):
    """
    Identify macrostates, then designate terminal & initial states.

    Returns the GPCCA estimator with terminal states set.
    """
    print("\n--- Macrostates & terminal states (GPCCA) ---")

    # Schur decomposition first — gives us the eigenvalue spectrum
    g.compute_schur(n_components=20)
    try:
        # CellRank 2 plot_spectrum() returns None in headless/non-interactive
        # environments. Call it then grab whatever matplotlib has open.
        g.plot_spectrum(real_only=True)
        fig = plt.gcf()
        fig_path = f'{output_dir}/{subset_name}_spectrum.pdf'
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close('all')
        print(f"Saved: {fig_path}")
    except Exception as e:
        plt.close('all')
        print(f"Spectrum plot skipped: {e}")

    # Number of macrostates
    if n_macrostates == 'auto':
        # Use the eigengap heuristic — look for a gap in the real Schur spectrum
        try:
            eigvals = np.real(g.eigendecomposition['D'])
            gaps = np.diff(eigvals[:15])
            # Largest gap among the first ~10 (skip first since it's always near 1)
            n_macrostates = int(np.argmax(np.abs(gaps[1:11])) + 2)
            print(f"Auto-selected n_macrostates = {n_macrostates} (largest eigengap)")
        except Exception:
            n_macrostates = 5
            print(f"Eigengap detection failed, defaulting to {n_macrostates}")

    print(f"Computing {n_macrostates} macrostates...")
    g.compute_macrostates(n_states=n_macrostates, cluster_key=cell_type_col)

    # Plot macrostates
    try:
        fig = g.plot_macrostates(which='all', basis='umap', show=False)
        if fig is not None:
            fig.tight_layout()
            fig.savefig(f'{output_dir}/{subset_name}_macrostates.pdf',
                        dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_dir}/{subset_name}_macrostates.pdf")
    except Exception as e:
        print(f"Macrostate plot skipped: {e}")

    # Terminal states
    print("\nDesignating terminal states...")

    # GPCCA appends suffixes to macrostate names (e.g. AS3 -> AS3_1, AS3_2)
    # and may split clusters across multiple macrostates. Always work from the
    # actual macrostate names rather than hardcoded cluster labels.
    available_macrostates = list(g.macrostates.cat.categories)
    print(f"Available macrostates: {available_macrostates}")

    if terminal_states == 'auto':
        try:
            g.predict_terminal_states()
            print(f"Auto-selected terminal states: {list(g.terminal_states.cat.categories)}")
        except ValueError as e:
            if "No macrostates have been selected" in str(e):
                # GPCCA couldn't auto-identify terminal states — common for
                # non-directed lineages (e.g. microglia activation states)
                # where all macrostates are equally metastable.
                # Fall back: use all macrostates as terminal.
                print(f"  WARNING: Auto terminal state selection failed: {e}")
                print(f"  Falling back: setting ALL macrostates as terminal.")
                g.set_terminal_states(
                    states=available_macrostates,
                    allow_overlap=True,
                )
                print(f"  Terminal states set to: {available_macrostates}")
            else:
                raise
    else:
        # terminal_states is a list of cluster name prefixes (e.g. ['AS1','AS2','AS3']).
        # Match against actual macrostate names which may have suffixes like '_1', '_2'.
        # This way ['AS1', 'AS2', 'AS3'] correctly matches ['AS1', 'AS3_1', 'AS3_2'].
        matched = []
        for requested in terminal_states:
            hits = [m for m in available_macrostates
                    if m == requested or m.startswith(requested + '_')]
            if hits:
                matched.extend(hits)
                print(f"  '{requested}' -> {hits}")
            else:
                print(f"  WARNING: '{requested}' matched no macrostates, skipping")

        if not matched:
            # Last resort: exclude anything containing the initial_states prefixes
            # and UP/I-style intermediate clusters
            exclude_prefixes = (['UP/I', 'UPI'] +
                                ([initial_states] if isinstance(initial_states, str)
                                 else (initial_states if initial_states != 'auto' else [])))
            matched = [m for m in available_macrostates
                       if not any(m.startswith(p) for p in exclude_prefixes)]
            print(f"  No requested states matched. Falling back to: {matched}")

        print(f"Setting terminal states: {matched}")
        g.set_terminal_states(states=matched, allow_overlap=True)

    try:
        fig = g.plot_macrostates(which='terminal', basis='umap', show=False)
        if fig is not None:
            fig.tight_layout()
            fig.savefig(f'{output_dir}/{subset_name}_terminal_states.pdf',
                        dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_dir}/{subset_name}_terminal_states.pdf")
    except Exception as e:
        print(f"Terminal-state plot skipped: {e}")

    # Initial states
    print("\nDesignating initial states...")
    try:
        if initial_states == 'auto':
            g.predict_initial_states(allow_overlap=True)
            print(f"Auto-selected initial states: {list(g.initial_states.cat.categories)}")
        else:
            # Same prefix-matching as terminal states above
            matched_init = []
            for requested in ([initial_states] if isinstance(initial_states, str)
                              else initial_states):
                hits = [m for m in available_macrostates
                        if m == requested or m.startswith(requested + '_')]
                if hits:
                    matched_init.extend(hits)
                    print(f"  '{requested}' -> {hits}")
                else:
                    print(f"  WARNING: '{requested}' matched no macrostates, skipping")

            if matched_init:
                print(f"Setting initial states: {matched_init}")
                g.set_initial_states(states=matched_init, allow_overlap=True)
                print(f"Set initial states: {matched_init}")
            else:
                # Fallback: auto-predict when no manual states matched
                print("  No manual initial states matched — falling back to auto-prediction")
                g.predict_initial_states(allow_overlap=True)
                print(f"Auto-selected initial states: {list(g.initial_states.cat.categories)}")

        fig = g.plot_macrostates(which='initial', basis='umap', show=False)
        if fig is not None:
            fig.tight_layout()
            fig.savefig(f'{output_dir}/{subset_name}_initial_states.pdf',
                        dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_dir}/{subset_name}_initial_states.pdf")
    except Exception as e:
        print(f"Initial-state assignment skipped: {e}")

    return g


def _inject_lineage(g, fate_matrix, term_names):
    """
    Inject a scipy-computed fate_matrix into the GPCCA estimator as a Lineage object.

    CellRank downstream methods (compute_lineage_drivers, gene_trends) check
    g.fate_probabilities, which is only populated by g.compute_fate_probabilities().
    Since we bypass that call (it ignores solver= without petsc4py), we inject
    directly into g._fate_probabilities and adata.obsm['lineages_fwd'].
    """
    try:
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as _plt
        from cellrank._utils._lineage import Lineage

        cmap = _plt.get_cmap('tab10')
        colors = [mcolors.to_hex(cmap(i % 10)) for i in range(len(term_names))]

        lin = Lineage(
            fate_matrix,
            names=list(term_names),
            colors=colors,
        )
        g.adata.obsm['lineages_fwd'] = lin
        g._fate_probabilities = lin
        print("  Injected fate probabilities into g._fate_probabilities "
              "and adata.obsm['lineages_fwd']")
    except Exception as e:
        print(f"  WARNING: Lineage injection failed ({e}). "
              f"Lineage drivers and gene trends will be skipped.")


def _plot_fate_probs(adata, fate_matrix, term_names, output_dir, subset_name):
    """Plot fate probability UMAPs — shared by both computation paths."""
    umap = adata.obsm['X_umap']
    xlim = (umap[:, 0].min() - 1, umap[:, 0].max() + 1)
    ylim = (umap[:, 1].min() - 1, umap[:, 1].max() + 1)
    n_terms = len(term_names)
    ncols = min(3, n_terms)
    nrows = (n_terms + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.array(axes).flatten() if n_terms > 1 else [axes]
    for k, tname in enumerate(term_names):
        ax = axes[k]
        vals = fate_matrix[:, k]
        sc = ax.scatter(umap[:, 0], umap[:, 1], c=vals,
                        cmap='viridis', s=3, alpha=0.7, vmin=0, vmax=1)
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_title(f'Fate: {tname}', fontsize=10)
        ax.set_xlabel('UMAP_1'); ax.set_ylabel('UMAP_2')
        plt.colorbar(sc, ax=ax)
    for k in range(n_terms, len(axes)):
        axes[k].set_visible(False)
    fig.suptitle(f'{subset_name} — Fate Probabilities', fontsize=12)
    fig.tight_layout()
    fig_path = f'{output_dir}/{subset_name}_fate_probabilities.pdf'
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fig_path}")


def run_fate_probabilities(g, output_dir, subset_name):
    """
    Compute per-cell fate probabilities.

    CellRank's compute_fate_probabilities() ignores the solver= parameter
    when petsc4py is unavailable — it always falls back to gmres internally
    regardless of what is passed. gmres without a preconditioner diverges on
    ill-conditioned transition matrices.

    Fix: extract the transition matrix and terminal state membership directly
    from the estimator and compute absorption probabilities ourselves using
    scipy.sparse.linalg.spsolve (direct LU), which has no convergence issues.

    Math: for a Markov chain with absorbing states, the absorption probability
    matrix Q satisfies (I - T_QQ) @ Q = T_QR, where T_QQ is the transient-to-
    transient block and T_QR is the transient-to-absorbing block.
    """
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    print("\n--- Fate probabilities (direct scipy solver) ---")

    # --- Extract transition matrix and terminal state membership ---
    T = g.transition_matrix          # sparse (n_cells x n_cells)
    memberships = g.macrostates_memberships  # (n_cells x n_macrostates) Lineage obj

    # Terminal state column indices in the membership matrix
    term_names = list(g.terminal_states.cat.categories)
    all_macro_names = list(memberships.names)
    print(f"Terminal states: {term_names}")

    term_idx = [all_macro_names.index(n) for n in term_names]
    M = memberships.X  # dense (n_cells x n_macrostates)

    # Classify absorbing vs transient using a SOFT threshold on terminal membership.
    # A cell is absorbing for terminal state k if its membership in k is >= threshold.
    # Using hard dominant-state assignment fails when all macrostates are terminal
    # (every cell gets classified absorbing, leaving no transient cells).
    ABSORBING_THRESHOLD = 0.5
    term_M = M[:, term_idx]  # (n_cells x n_terminal_states)
    is_absorbing = term_M.max(axis=1) >= ABSORBING_THRESHOLD
    is_transient = ~is_absorbing

    n_absorbing = is_absorbing.sum()
    n_transient = is_transient.sum()
    print(f"  Absorbing threshold: {ABSORBING_THRESHOLD}")
    print(f"  Absorbing cells: {n_absorbing}, Transient cells: {n_transient}")

    # If still all absorbing (very tight macrostates), lower the threshold
    if n_transient == 0:
        ABSORBING_THRESHOLD = 0.8
        is_absorbing = term_M.max(axis=1) >= ABSORBING_THRESHOLD
        is_transient = ~is_absorbing
        n_absorbing = is_absorbing.sum()
        n_transient = is_transient.sum()
        print(f"  Retrying with threshold {ABSORBING_THRESHOLD}: "
              f"absorbing={n_absorbing}, transient={n_transient}")

    # Final fallback: use macrostate memberships directly as fate probabilities.
    # This is valid when macrostates aren't well-separated enough for the absorbing
    # chain formulation — memberships are a reasonable soft approximation.
    if n_transient == 0 or n_absorbing == 0:
        print("  WARNING: Cannot separate absorbing/transient cells. "
              "Using macrostate memberships directly as fate probability approximation.")
        fate_matrix = term_M / term_M.sum(axis=1, keepdims=True).clip(min=1e-10)
        fate_df = pd.DataFrame(fate_matrix, index=g.adata.obs_names, columns=term_names)
        for col in term_names:
            g.adata.obs[f'fate_prob_{col}'] = fate_df[col].values
        csv_path = f'{output_dir}/{subset_name}_fate_probabilities.csv'
        fate_df.to_csv(csv_path)
        print(f"  Saved: {csv_path}")
        _inject_lineage(g, fate_matrix, term_names)
        _plot_fate_probs(g.adata, fate_matrix, term_names, output_dir, subset_name)
        return

    # Reorder: transient first, absorbing second (standard absorbing chain form)
    trans_idx = np.where(is_transient)[0]
    abs_idx = np.where(is_absorbing)[0]

    T_dense = np.array(T.toarray() if sp.issparse(T) else T)

    T_QQ = T_dense[np.ix_(trans_idx, trans_idx)]  # transient -> transient
    T_QR = T_dense[np.ix_(trans_idx, abs_idx)]    # transient -> absorbing

    print(f"  T_QQ: {T_QQ.shape}, T_QR: {T_QR.shape}")

    # Solve: (I - T_QQ) @ B = T_QR  for absorption probability matrix B
    # B[i, j] = probability that a random walk from transient cell i
    # is absorbed by absorbing cell j
    print("  Solving absorption system with scipy spsolve (direct LU)...")
    I_QQ = np.eye(T_QQ.shape[0])
    A = I_QQ - T_QQ   # (I - T_QQ)
    A_sparse = sp.csr_matrix(A)

    # Solve column by column (one per terminal state region)
    n_terms = T_QR.shape[1]
    B = np.zeros((n_transient, n_terms))
    for j in range(n_terms):
        rhs = T_QR[:, j]
        B[:, j] = spla.spsolve(A_sparse, rhs)

    # Normalize rows to sum to 1 (handle any small numerical errors)
    row_sums = B.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1, row_sums)
    B = B / row_sums

    print(f"  Solved. B range: [{B.min():.4f}, {B.max():.4f}]  "
          f"Row sum range: [{B.sum(axis=1).min():.4f}, {B.sum(axis=1).max():.4f}]")

    # --- Map fate probs back to per-terminal-state probabilities ---
    n_cells = T_dense.shape[0]
    fate_matrix = np.zeros((n_cells, len(term_names)))

    # Absorbing cells: assign fate prob proportional to their terminal membership
    for k in range(len(term_names)):
        fate_matrix[abs_idx, k] = term_M[abs_idx, k]

    # Normalize absorbing cells
    abs_row_sums = fate_matrix[abs_idx].sum(axis=1, keepdims=True).clip(min=1e-10)
    fate_matrix[abs_idx] = fate_matrix[abs_idx] / abs_row_sums

    # Transient cells: distribute B columns to terminal macrostates
    # Each column of B corresponds to an absorbing cell; aggregate by terminal state
    abs_term_membership = term_M[abs_idx]  # (n_absorbing x n_terminal_states)
    for k in range(len(term_names)):
        weights = abs_term_membership[:, k]
        fate_matrix[trans_idx, k] = B @ weights

    # Renormalize all rows
    row_sums = fate_matrix.sum(axis=1, keepdims=True).clip(min=1e-10)
    fate_matrix = fate_matrix / row_sums

    # Store in adata.obs
    fate_df = pd.DataFrame(
        fate_matrix,
        index=g.adata.obs_names,
        columns=term_names,
    )
    for col in term_names:
        g.adata.obs[f'fate_prob_{col}'] = fate_df[col].values

    # Inject into estimator so downstream CellRank methods work
    _inject_lineage(g, fate_matrix, term_names)

    # Save CSV
    csv_path = f'{output_dir}/{subset_name}_fate_probabilities.csv'
    fate_df.to_csv(csv_path)
    print(f"  Saved: {csv_path}")

    # Plot via shared helper
    _plot_fate_probs(g.adata, fate_matrix, term_names, output_dir, subset_name)

    try:
        fig = g.plot_fate_probabilities(basis='umap', same_plot=False, show=False)
        if fig is not None:
            fig.tight_layout()
            fig.savefig(f'{output_dir}/{subset_name}_fate_probabilities_cr.pdf',
                        dpi=300, bbox_inches='tight')
            plt.close(fig)
    except Exception as e:
        print(f"  CellRank fate prob plot skipped: {e}")

    # Save fate probabilities as CSV (per-cell, per-lineage)
    try:
        fate_df = g.fate_probabilities.X
        lineage_names = list(g.fate_probabilities.names)
        fate_pdf = pd.DataFrame(
            fate_df,
            index=g.adata.obs_names,
            columns=[f'fate_prob_{n}' for n in lineage_names],
        )
        csv_path = f'{output_dir}/{subset_name}_fate_probabilities.csv'
        fate_pdf.to_csv(csv_path)
        print(f"Saved: {csv_path}")
    except Exception as e:
        print(f"Fate probability CSV skipped: {e}")


def run_lineage_drivers(g, output_dir, subset_name, n_drivers=50):
    """Find genes correlated with commitment to each terminal fate."""
    print(f"\n--- Lineage drivers (top {n_drivers} per lineage) ---")
    try:
        drivers = g.compute_lineage_drivers()  # DataFrame: genes x lineage stats

        # Save full table
        full_path = f'{output_dir}/{subset_name}_lineage_drivers_full.csv'
        drivers.to_csv(full_path)
        print(f"Saved: {full_path}")

        # Save top N per lineage
        lineage_names = [c.split('_')[0] for c in drivers.columns
                         if c.endswith('_corr')]
        top_records = []
        for lin in lineage_names:
            corr_col = f'{lin}_corr'
            qval_col = f'{lin}_qval'
            if corr_col not in drivers.columns:
                continue
            top = drivers.nlargest(n_drivers, corr_col)[[corr_col, qval_col]].copy()
            top.columns = ['correlation', 'qvalue']
            top['lineage'] = lin
            top['gene'] = top.index
            top_records.append(top.reset_index(drop=True))

        if top_records:
            top_df = pd.concat(top_records, ignore_index=True)
            top_path = f'{output_dir}/{subset_name}_lineage_drivers_top.csv'
            top_df.to_csv(top_path, index=False)
            print(f"Saved: {top_path}")

            # Plot top drivers per lineage using custom binned approach
            # (cr.pl.gene_trends fails with injected Lineage objects)
            for lin in lineage_names:
                lin_safe = lin.replace('/', '_').replace('\\', '_').replace(' ', '_')
                corr_col = f'{lin}_corr'
                if corr_col not in drivers.columns:
                    continue
                top_genes = [g for g in drivers.nlargest(6, corr_col).index
                             if g in g.adata.var_names] if False else \
                            [gene for gene in drivers.nlargest(6, corr_col).index
                             if gene in g.adata.var_names]
                if not top_genes:
                    continue
                try:
                    # Re-use the custom trend plotting logic
                    import scipy.sparse as sp2
                    from scipy.stats import binned_statistic as bstat
                    import warnings as _w
                    fate_col = f'fate_prob_{lin}'
                    if fate_col not in g.adata.obs.columns:
                        continue
                    pt = g.adata.obs['velocity_pseudotime'].values.astype(float)
                    valid = np.isfinite(pt)
                    fate_w = g.adata.obs[fate_col].values.astype(float)
                    use = (fate_w > 0.05) & valid
                    if use.sum() < 10:
                        continue
                    pt_use = pt[use]
                    w_use = fate_w[use]
                    X_drv = g.adata[:, top_genes].X
                    if sp2.issparse(X_drv):
                        X_drv = X_drv.toarray()
                    X_drv = np.array(X_drv[use], dtype=float)
                    bins = np.linspace(pt_use.min(), pt_use.max(), 31)
                    bin_centers = (bins[:-1] + bins[1:]) / 2
                    ncols = min(3, len(top_genes))
                    nrows = (len(top_genes) + ncols - 1) // ncols
                    fig, axes = plt.subplots(nrows, ncols,
                                             figsize=(5*ncols, 3.5*nrows),
                                             squeeze=False)
                    for gi, gene in enumerate(top_genes):
                        ax = axes.flatten()[gi]
                        expr = X_drv[:, gi]
                        with _w.catch_warnings():
                            _w.simplefilter('ignore')
                            wsum, _, _ = bstat(pt_use, expr*w_use, 'sum', bins=bins)
                            wbin, _, _ = bstat(pt_use, w_use, 'sum', bins=bins)
                        with np.errstate(invalid='ignore', divide='ignore'):
                            trend = np.where(wbin > 0, wsum/wbin, np.nan)
                        idx = np.random.choice(len(pt_use), min(300, len(pt_use)), replace=False)
                        ax.scatter(pt_use[idx], expr[idx], c=w_use[idx],
                                   cmap='Blues', s=4, alpha=0.4, vmin=0, vmax=1)
                        vb = np.isfinite(trend)
                        if vb.sum() >= 3:
                            ax.plot(bin_centers[vb], trend[vb], color='crimson', lw=2)
                        ax.set_title(gene, fontsize=9, fontweight='bold')
                        ax.set_xlabel('Pseudotime', fontsize=7)
                        ax.set_ylabel('Expression', fontsize=7)
                    for gi in range(len(top_genes), len(axes.flatten())):
                        axes.flatten()[gi].set_visible(False)
                    fig.suptitle(f'{subset_name} — Top Drivers: {lin}', fontsize=10)
                    fig.tight_layout()
                    fig.savefig(f'{output_dir}/{subset_name}_drivers_{lin_safe}.pdf',
                                dpi=200, bbox_inches='tight')
                    plt.close(fig)
                    print(f"Saved: {output_dir}/{subset_name}_drivers_{lin_safe}.pdf")
                except Exception as e:
                    print(f"Driver gene plot for {lin} skipped: {e}")

    except Exception as e:
        print(f"Lineage drivers failed: {e}")


def run_gene_trends(g, marker_genes, output_dir, subset_name):
    """
    Plot expression trends along pseudotime using custom scipy binning.

    CellRank's cr.pl.gene_trends fails with injected Lineage objects because
    its internal fitting pipeline expects fate_probabilities set via
    compute_fate_probabilities(). This implementation bypasses that entirely:
    for each lineage, it weights cells by their fate probability, bins
    pseudotime, computes weighted mean expression per bin, and plots the trend.
    """
    from scipy.stats import binned_statistic
    import scipy.sparse as sp
    import warnings

    print(f"\n--- Gene trends along pseudotime (custom) ---")

    available = [g_ for g_ in marker_genes if g_ in g.adata.var_names]
    missing = set(marker_genes) - set(available)
    if missing:
        print(f"Markers not in var_names (skipped): {sorted(missing)}")
    if not available:
        print("No marker genes available, skipping gene trends.")
        return

    adata = g.adata

    if 'velocity_pseudotime' not in adata.obs.columns:
        print("  velocity_pseudotime not in adata.obs, skipping.")
        return

    pt = adata.obs['velocity_pseudotime'].values.astype(float)
    valid = np.isfinite(pt)
    print(f"  Cells with valid pseudotime: {valid.sum()} / {len(pt)}")
    if valid.sum() < 10:
        print("  Too few valid cells, skipping.")
        return

    fate_cols = [c for c in adata.obs.columns if c.startswith('fate_prob_')]
    if not fate_cols:
        print("  No fate_prob_ columns found, skipping.")
        return
    lineage_names = [c.replace('fate_prob_', '') for c in fate_cols]
    print(f"  Lineages: {lineage_names}")
    print(f"  Plotting trends for {len(available)} markers...")

    X = adata[:, available].X
    if sp.issparse(X):
        X = X.toarray()
    X = np.array(X, dtype=float)

    n_bins = 30

    def sanitize(name):
        """Replace filesystem-unsafe characters for use in filenames."""
        return name.replace('/', '_').replace('\\', '_').replace(' ', '_')

    # ---- PER-LINEAGE TREND PLOTS ----
    for lin_name, fate_col in zip(lineage_names, fate_cols):
        lin_safe = sanitize(lin_name)
        fate_w = adata.obs[fate_col].values.astype(float)
        use = (fate_w > 0.05) & valid
        if use.sum() < 10:
            print(f"  Skipping {lin_name}: too few cells with fate_prob > 0.05")
            continue

        pt_use = pt[use]
        X_use = X[use]
        w_use = fate_w[use]

        bins = np.linspace(pt_use.min(), pt_use.max(), n_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        ncols = min(4, len(available))
        nrows = (len(available) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(5 * ncols, 3.5 * nrows),
                                 squeeze=False)
        axes_flat = axes.flatten()

        for gi, gene in enumerate(available):
            ax = axes_flat[gi]
            expr = X_use[:, gi]

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                wsum, _, _ = binned_statistic(pt_use, expr * w_use,
                                              statistic='sum', bins=bins)
                wbin, _, _ = binned_statistic(pt_use, w_use,
                                              statistic='sum', bins=bins)

            with np.errstate(invalid='ignore', divide='ignore'):
                trend = np.where(wbin > 0, wsum / wbin, np.nan)

            # Scatter (subsample)
            idx = np.random.choice(len(pt_use),
                                   min(400, len(pt_use)), replace=False)
            sc = ax.scatter(pt_use[idx], expr[idx], c=w_use[idx],
                            cmap='Blues', s=4, alpha=0.4, vmin=0, vmax=1)

            valid_bins = np.isfinite(trend)
            if valid_bins.sum() >= 3:
                ax.plot(bin_centers[valid_bins], trend[valid_bins],
                        color='crimson', linewidth=2, zorder=5)

            ax.set_title(gene, fontsize=9, fontweight='bold')
            ax.set_xlabel('Pseudotime', fontsize=7)
            ax.set_ylabel('Expression', fontsize=7)
            ax.tick_params(labelsize=6)

        for gi in range(len(available), len(axes_flat)):
            axes_flat[gi].set_visible(False)

        fig.suptitle(f'{subset_name} — Gene Trends ({lin_name} fate)',
                     fontsize=11, fontweight='bold')
        fig.tight_layout()
        path = f'{output_dir}/{subset_name}_gene_trends_{lin_safe}.pdf'
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {path}")

    # ---- HEATMAP: genes × pseudotime bins, sorted by peak ----
    for lin_name, fate_col in zip(lineage_names, fate_cols):
        lin_safe = sanitize(lin_name)
        fate_w = adata.obs[fate_col].values.astype(float)
        use = (fate_w > 0.05) & valid
        if use.sum() < 10:
            continue

        pt_use = pt[use]
        X_use = X[use]
        w_use = fate_w[use]

        n_hbins = 50
        bins = np.linspace(pt_use.min(), pt_use.max(), n_hbins + 1)

        heatmap = np.zeros((len(available), n_hbins))
        for gi in range(len(available)):
            expr = X_use[:, gi]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                wsum, _, _ = binned_statistic(pt_use, expr * w_use,
                                              statistic='sum', bins=bins)
                wbin, _, _ = binned_statistic(pt_use, w_use,
                                              statistic='sum', bins=bins)
            with np.errstate(invalid='ignore', divide='ignore'):
                heatmap[gi] = np.where(wbin > 0, wsum / wbin, 0)

        # Z-score rows, sort by peak bin
        row_mean = heatmap.mean(axis=1, keepdims=True)
        row_std = heatmap.std(axis=1, keepdims=True).clip(min=1e-8)
        heatmap_z = (heatmap - row_mean) / row_std
        order = np.argsort(np.argmax(heatmap_z, axis=1))
        heatmap_z = heatmap_z[order]
        genes_sorted = [available[i] for i in order]

        fig, ax = plt.subplots(figsize=(12, max(3, len(available) * 0.4)))
        im = ax.imshow(heatmap_z, aspect='auto', cmap='RdBu_r',
                       vmin=-2, vmax=2, interpolation='nearest')
        ax.set_yticks(range(len(genes_sorted)))
        ax.set_yticklabels(genes_sorted, fontsize=8)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        ticks = np.linspace(0, n_hbins - 1, 5).astype(int)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{bin_centers[t]:.2f}' for t in ticks], fontsize=8)
        ax.set_xlabel('Pseudotime', fontsize=10)
        ax.set_title(f'{subset_name} — Marker Heatmap ({lin_name} fate, '
                     f'z-scored, sorted by peak)', fontsize=10)
        plt.colorbar(im, ax=ax, label='z-score')
        fig.tight_layout()
        path = f'{output_dir}/{subset_name}_heatmap_{lin_safe}.pdf'
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {path}")






# =============================================================================
# CLUSTER-SPECIFIC DRIVERS
# =============================================================================
def _plot_cluster_driver_trends(adata, top_genes, cluster_name, output_dir,
                                 subset_name, n_bins=30, n_scatter=300):
    """
    Plot expression trends along pseudotime for top drivers of a specific cluster.
    Cells are weighted by their binary cluster membership. Same custom binned
    approach used in run_gene_trends — no CellRank fitting machinery involved.
    """
    import scipy.sparse as sp2
    from scipy.stats import binned_statistic as bstat
    import warnings as _w

    if 'velocity_pseudotime' not in adata.obs.columns:
        print(f"    velocity_pseudotime not found, skipping trend plots.")
        return

    pt = adata.obs['velocity_pseudotime'].values.astype(float)
    valid = np.isfinite(pt)
    membership = (adata.obs.iloc[:, 0] == cluster_name).astype(float).values  # placeholder
    # Recompute properly from obs — passed in as weight array via caller
    # (see run_cluster_specific_drivers which calls this with correct membership)

    use = valid & (membership > 0)
    in_cluster = membership.astype(bool) & valid

    if in_cluster.sum() < 5:
        print(f"    Too few cells in cluster for trend plots.")
        return

    available = [g for g in top_genes if g in adata.var_names]
    if not available:
        return

    X = adata[:, available].X
    if sp2.issparse(X):
        X = X.toarray()
    X = np.array(X, dtype=float)

    pt_valid = pt[valid]
    X_valid = X[valid]
    membership_valid = membership[valid]

    bins = np.linspace(pt_valid.min(), pt_valid.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    ncols = min(3, len(available))
    nrows = (len(available) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 3.5 * nrows),
                             squeeze=False)
    axes_flat = axes.flatten()

    for gi, gene in enumerate(available):
        ax = axes_flat[gi]
        expr = X_valid[:, gi]

        with _w.catch_warnings():
            _w.simplefilter('ignore')
            wsum, _, _ = bstat(pt_valid, expr * membership_valid,
                               statistic='sum', bins=bins)
            wbin, _, _ = bstat(pt_valid, membership_valid,
                               statistic='sum', bins=bins)

        with np.errstate(invalid='ignore', divide='ignore'):
            trend = np.where(wbin > 0, wsum / wbin, np.nan)

        # Color points by cluster membership
        idx = np.random.choice(len(pt_valid),
                               min(n_scatter, len(pt_valid)), replace=False)
        colors = ['#d62728' if membership_valid[i] > 0.5 else '#aec7e8'
                  for i in idx]
        ax.scatter(pt_valid[idx], expr[idx], c=colors, s=4, alpha=0.5)

        valid_bins = np.isfinite(trend)
        if valid_bins.sum() >= 3:
            ax.plot(bin_centers[valid_bins], trend[valid_bins],
                    color='black', linewidth=2, zorder=5)

        ax.set_title(gene, fontsize=9, fontweight='bold')
        ax.set_xlabel('Pseudotime', fontsize=7)
        ax.set_ylabel('Expression', fontsize=7)
        ax.tick_params(labelsize=6)

    for gi in range(len(available), len(axes_flat)):
        axes_flat[gi].set_visible(False)

    safe = cluster_name.replace('/', '_').replace(' ', '_')
    fig.suptitle(f'{subset_name} — Top Drivers: {cluster_name} '
                 f'(red = in cluster)', fontsize=10, fontweight='bold')
    fig.tight_layout()
    path = f'{output_dir}/{subset_name}_cluster_drivers_{safe}_trends.pdf'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")


def run_cluster_specific_drivers(adata, cell_type_col, output_dir,
                                  subset_name, n_top=50):
    """
    Compute and plot lineage drivers for every cluster in cell_type_col.

    Uses Spearman correlation between gene expression and binary cluster
    membership (1 = in cluster, 0 = not), with Benjamini-Hochberg correction.
    Works for any cluster regardless of whether it is a GPCCA terminal state.

    Outputs per cluster:
      - CSV: ranked gene list with correlation and q-value
      - PDF: trend plots for top 6 driver genes along velocity_pseudotime,
             with in-cluster cells highlighted in red
    """
    import scipy.sparse as sp
    from scipy.stats import spearmanr, rankdata

    print(f"\n--- Cluster-specific drivers (all clusters) ---")

    if 'velocity_pseudotime' not in adata.obs.columns:
        print("  velocity_pseudotime not found, skipping.")
        return

    clusters = sorted(adata.obs[cell_type_col].unique())
    print(f"  Computing drivers for {len(clusters)} clusters: {clusters}")

    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.array(X, dtype=float)
    genes = list(adata.var_names)

    # Pre-filter to variable genes (std > 0) to skip constant genes
    gene_stds = X.std(axis=0)
    variable_mask = gene_stds > 1e-8
    X_var = X[:, variable_mask]
    genes_var = [g for g, m in zip(genes, variable_mask) if m]
    print(f"  Variable genes: {variable_mask.sum()} / {len(genes)}")

    for cluster in clusters:
        safe = cluster.replace('/', '_').replace(' ', '_')
        membership = (adata.obs[cell_type_col] == cluster).astype(float).values
        n_in = int(membership.sum())
        n_out = len(membership) - n_in

        print(f"\n  {cluster}: {n_in} cells in cluster, {n_out} outside")

        if n_in < 10 or n_out < 10:
            print(f"    Skipping — too few cells on one side.")
            continue

        # Spearman correlation for each variable gene vs membership
        results = []
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            for gi in range(X_var.shape[1]):
                rho, pval = spearmanr(X_var[:, gi], membership)
                results.append({'gene': genes_var[gi],
                                'correlation': rho,
                                'pvalue': float(pval)})

        df = pd.DataFrame(results).set_index('gene')

        # Benjamini-Hochberg FDR correction
        n = len(df)
        pvals = df['pvalue'].values
        ranks = rankdata(pvals)
        df['qvalue'] = np.minimum(1.0, pvals * n / ranks)
        df = df.sort_values('correlation', ascending=False)

        # Save full ranked list
        csv_path = f'{output_dir}/{subset_name}_cluster_drivers_{safe}.csv'
        df.to_csv(csv_path)
        print(f"    Saved: {csv_path}")

        # Top N summary
        n_sig = (df['qvalue'] < 0.05).sum()
        print(f"    Significant drivers (q<0.05): {n_sig}")
        print(f"    Top 10 positive drivers:")
        print(df.head(10)[['correlation', 'qvalue']].to_string())

        # Top drivers for trend plots (top 6 positively correlated)
        top_positive = df[df['correlation'] > 0].head(6).index.tolist()
        top_negative = df[df['correlation'] < 0].tail(6).index.tolist()

        # Positive drivers (genes enriched IN cluster)
        if top_positive:
            _plot_cluster_driver_trends_direct(
                adata=adata,
                X_full=X,
                genes_full=genes,
                top_genes=top_positive,
                membership=membership,
                cluster_name=cluster,
                direction='positive',
                output_dir=output_dir,
                subset_name=subset_name,
                driver_stats=df,
            )

        # Negative drivers (genes depleted in cluster = enriched outside)
        if top_negative:
            _plot_cluster_driver_trends_direct(
                adata=adata,
                X_full=X,
                genes_full=genes,
                top_genes=top_negative,
                membership=membership,
                cluster_name=cluster,
                direction='negative',
                output_dir=output_dir,
                subset_name=subset_name,
                driver_stats=df,
            )


def _plot_cluster_driver_trends_direct(adata, X_full, genes_full, top_genes,
                                        membership, cluster_name, direction,
                                        output_dir, subset_name, n_bins=30,
                                        driver_stats=None):
    """
    Publication-quality pseudotime expression trend plots for driver genes.

    Changes from exploratory version:
    - Font sizes scaled for journal column width (labels 11pt, titles 12pt)
    - 300 dpi output
    - Legend instead of suptitle annotation for cluster coloring
    - Spearman rho + q-value annotated per panel
    - Top/right spines removed
    - Larger point size and better alpha
    - Human-readable titles (underscores replaced)

    driver_stats : pd.DataFrame, optional
        DataFrame indexed by gene with 'correlation' and 'qvalue' columns.
        Used to annotate each panel with rho and q. Pass the full drivers
        DataFrame from run_cluster_specific_drivers.
    """
    from scipy.stats import binned_statistic as bstat
    from matplotlib.lines import Line2D
    import warnings as _w

    # ---- publication style constants ----
    LABEL_FS   = 11
    TICK_FS    = 10
    TITLE_FS   = 12
    ANNOT_FS   = 9
    SUPTITLE_FS = 13
    POINT_SIZE  = 8
    POINT_ALPHA = 0.45
    TREND_LW    = 2.0
    DPI         = 300
    # Color-blind friendly: vermillion / steel blue
    COLOR_IN    = '#E64B35'   # in-cluster
    COLOR_OUT   = '#4DBBD5'   # outside cluster

    pt = adata.obs['velocity_pseudotime'].values.astype(float)
    valid = np.isfinite(pt)
    if valid.sum() < 10:
        return

    available_idx = [genes_full.index(g) for g in top_genes if g in genes_full]
    available = [genes_full[i] for i in available_idx]
    if not available:
        return

    pt_v = pt[valid]
    X_v = X_full[np.ix_(valid, available_idx)]
    mem_v = membership[valid]

    bins = np.linspace(pt_v.min(), pt_v.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    ncols = min(3, len(available))
    nrows = (len(available) + ncols - 1) // ncols

    # Scale panel size for readability at journal width
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.5 * ncols, 3.8 * nrows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    n_scatter = min(500, len(pt_v))
    scatter_idx = np.random.choice(len(pt_v), n_scatter, replace=False)
    point_colors = np.where(mem_v[scatter_idx] > 0.5, COLOR_IN, COLOR_OUT)

    for gi, (gene, col_idx) in enumerate(zip(available, range(len(available)))):
        ax = axes_flat[gi]
        expr = X_v[:, col_idx]

        # Weighted trend
        w = mem_v if direction == 'positive' else (1 - mem_v)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            wsum, _, _ = bstat(pt_v, expr * w, statistic='sum', bins=bins)
            wbin, _, _ = bstat(pt_v, w, statistic='sum', bins=bins)
        with np.errstate(invalid='ignore', divide='ignore'):
            trend = np.where(wbin > 0, wsum / wbin, np.nan)

        # Scatter
        ax.scatter(pt_v[scatter_idx], expr[scatter_idx],
                   c=point_colors, s=POINT_SIZE, alpha=POINT_ALPHA,
                   linewidths=0, rasterized=True)

        # Trend line
        vb = np.isfinite(trend)
        if vb.sum() >= 3:
            ax.plot(bin_centers[vb], trend[vb],
                    color='black', lw=TREND_LW, zorder=5)

        # Gene title
        ax.set_title(gene, fontsize=TITLE_FS, fontweight='bold', pad=4)

        # Axis labels — only leftmost column gets y-label, only bottom row gets x-label
        row_i = gi // ncols
        col_i = gi % ncols
        if col_i == 0:
            ax.set_ylabel('Normalized expression', fontsize=LABEL_FS)
        else:
            ax.set_ylabel('')
        if row_i == nrows - 1 or gi >= len(available) - ncols:
            ax.set_xlabel('Pseudotime', fontsize=LABEL_FS)
        else:
            ax.set_xlabel('')

        ax.tick_params(labelsize=TICK_FS)

        # Remove top/right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Annotate with Spearman rho and q-value if stats provided
        if driver_stats is not None and gene in driver_stats.index:
            rho = driver_stats.loc[gene, 'correlation']
            qval = driver_stats.loc[gene, 'qvalue']
            # Format q-value
            if qval < 0.001:
                q_str = 'q < 0.001'
            elif qval < 0.01:
                q_str = f'q = {qval:.3f}'
            else:
                q_str = f'q = {qval:.2f}'
            annot = f'ρ = {rho:.2f}, {q_str}'
            ax.text(0.97, 0.97, annot,
                    transform=ax.transAxes,
                    ha='right', va='top',
                    fontsize=ANNOT_FS,
                    color='#333333',
                    bbox=dict(boxstyle='round,pad=0.2',
                              facecolor='white', alpha=0.7,
                              edgecolor='none'))

    # Hide unused panels
    for gi in range(len(available), len(axes_flat)):
        axes_flat[gi].set_visible(False)

    # Legend (replaces suptitle annotation)
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_IN, markersize=7,
               label=f'{cluster_name}'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_OUT, markersize=7,
               label='Other clusters'),
        Line2D([0], [0], color='black', lw=2,
               label='Weighted trend'),
    ]
    fig.legend(handles=legend_elements,
               loc='lower center',
               ncol=3,
               fontsize=ANNOT_FS + 1,
               frameon=False,
               bbox_to_anchor=(0.5, -0.02))

    # Clean suptitle — no internal notes, human-readable names
    readable_subset = subset_name.replace('_', ' ').title()
    dir_label_human = ('enriched in' if direction == 'positive'
                       else 'depleted in')
    fig.suptitle(
        f'{readable_subset} — Top genes {dir_label_human} {cluster_name}',
        fontsize=SUPTITLE_FS,
        fontweight='bold',
        y=1.01,
    )

    fig.tight_layout(rect=[0, 0.04, 1, 1])

    safe = cluster_name.replace('/', '_').replace(' ', '_')
    dir_label = 'enriched_in' if direction == 'positive' else 'depleted_in'
    path = (f'{output_dir}/{subset_name}_cluster_drivers_'
            f'{safe}_{dir_label}.pdf')
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")



def run_cellrank_for_subset(config):
    """Full CellRank 2 pipeline for one subset."""
    h5ad = config['h5ad_file']
    subset_name = config['subset_name']
    cell_type_col = config['cell_type_col']

    print(f"\n{'='*60}")
    print(f"CELLRANK 2 — {subset_name}")
    print(f"Reading: {h5ad}")
    print(f"{'='*60}")

    if not os.path.exists(h5ad):
        print(f"ERROR: {h5ad} does not exist. Run scVelo first.")
        return

    output_dir = f'{OUTPUT_DIR_BASE}/{subset_name}'
    os.makedirs(output_dir, exist_ok=True)

    adata = sc.read_h5ad(h5ad)
    print(f"Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")

    # Make sure cell_type_col is categorical (CellRank expects it)
    if not isinstance(adata.obs[cell_type_col].dtype, pd.CategoricalDtype):
        adata.obs[cell_type_col] = adata.obs[cell_type_col].astype('category')

    # Apply any custom category order specified in config.
    # The subset palette assigns colors POSITIONALLY, so category order must
    # match the order used in the R UMAP (i.e. the Seurat factor level order).
    # If not set, defaults to sorted alphabetical — which may mismatch Seurat.
    custom_order = config.get('category_order', None)
    if custom_order is not None:
        current = list(adata.obs[cell_type_col].cat.categories)
        missing = [c for c in custom_order if c not in current]
        extra   = [c for c in current if c not in custom_order]
        if missing:
            print(f"  WARNING: category_order contains unknown categories: {missing}")
        # Append any categories not listed in custom_order at the end
        full_order = custom_order + extra
        full_order = [c for c in full_order if c in current]
        adata.obs[cell_type_col] = adata.obs[cell_type_col].cat.reorder_categories(full_order)
        print(f"  Category order set to: {full_order}")
    else:
        print(f"  Category order (default alphabetical): "
              f"{sorted(adata.obs[cell_type_col].unique())}")

    # 1. Velocity-directed PAGA
    run_paga_velocity(adata, cell_type_col, output_dir, subset_name)

    # 2. Build kernel (velocity + connectivity)
    kernel = build_kernel(adata, velocity_weight=0.8)

    # 3. GPCCA estimator → macrostates → terminal/initial states
    g = GPCCA(kernel)
    g = run_terminal_states(
        g, adata, cell_type_col, output_dir, subset_name,
        terminal_states=config['terminal_states'],
        initial_states=config['initial_states'],
        n_macrostates=config['n_macrostates'],
    )

    # 4. Fate probabilities
    run_fate_probabilities(g, output_dir, subset_name)

    # 5. Lineage drivers (fate-probability based)
    run_lineage_drivers(g, output_dir, subset_name,
                        n_drivers=config['n_drivers'])

    # 6. Cluster-specific drivers (Spearman, all clusters)
    run_cluster_specific_drivers(
        adata=g.adata,
        cell_type_col=cell_type_col,
        output_dir=output_dir,
        subset_name=subset_name,
        n_top=config.get('n_drivers', 50),
    )

    # 7. Gene trends along pseudotime
    run_gene_trends(g, config['gene_trend_markers'], output_dir, subset_name)

    # Save final AnnData with all CellRank annotations
    out_h5ad = f'{output_dir}/{subset_name}_cellrank.h5ad'
    g.adata.write(out_h5ad)
    print(f"\nSaved: {out_h5ad}")
    print(f"\n{subset_name} — DONE\n")


# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == '__main__':
    print("CellRank 2 fate analysis on fetal brain scVelo outputs")
    print(f"Analyses queued: {len(ANALYSES)}\n")

    for config in ANALYSES:
        try:
            run_cellrank_for_subset(config)
        except Exception as e:
            print(f"\nFAILED for {config['subset_name']}: {e}")
            import traceback
            traceback.print_exc()
            print("Continuing to next analysis...\n")

    print("\n" + "=" * 60)
    print("ALL CELLRANK ANALYSES DONE")
    print("=" * 60)
