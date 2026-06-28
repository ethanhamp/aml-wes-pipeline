"""
normalization_diagnostics.py
============================
Characterize the normalization mismatch between BeatAML training data and
the Eisfeld expression matrix, restricted to the LASSO signature gene set.

Outputs (saved to Figures/):
  norm_diag_mean_scatter.png  — per-gene mean: BeatAML (x) vs Eisfeld (y)
                                colored by whether gene has nonzero LASSO coef
  norm_diag_distributions.png — density of all expression values across
                                signature genes, BeatAML vs Eisfeld overlaid

Printed stats:
  - Spearman r of per-gene means
  - Median absolute difference in means
  - % genes where Eisfeld mean > BeatAML mean
"""

# ---------------------------------------------------------------------------
# 0.  PATHS
# ---------------------------------------------------------------------------
import os
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # headless — no display needed
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import spearmanr, gaussian_kde

BASE_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
BEATAML_FILE = Path(os.environ.get("BEATAML_DIR", str(Path(os.environ.get("ARHG_BASE_DIR", "..")).parent / "BeatAML"))) / "beataml_waves1to4_norm_exp_dbgap.txt"
EISFELD_FILE  = BASE_DIR / "eisfeld_expression.tsv"
MODEL_PKL     = BASE_DIR / "ras_lasso_model.pkl"
FEAT_GENES    = BASE_DIR / "signature_feature_genes.txt"
FIGURES_DIR   = BASE_DIR / "Figures"
FIGURES_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# 1.  LOAD GENE LIST
# ---------------------------------------------------------------------------
print("=" * 65)
print("normalization_diagnostics.py")
print("=" * 65)

with open(FEAT_GENES) as fh:
    feature_genes = [line.strip() for line in fh if line.strip()]
feature_set = set(feature_genes)
print(f"\nSignature gene list : {len(feature_genes):,} genes loaded")

# ---------------------------------------------------------------------------
# 2.  EXTRACT LASSO NONZERO GENES
#     LogisticRegressionCV with L1/ElasticNet penalty stores coef_ as
#     shape (1, n_features) for binary classification.  We align by the
#     same feature_genes ordering used during training.
# ---------------------------------------------------------------------------
with open(MODEL_PKL, "rb") as fh:
    model = pickle.load(fh)

coef_vec = model.coef_.ravel()             # shape: (n_features,) = (15167,)
# feature_genes list IS the model's feature order (same file used in training)
lasso_nonzero = set(
    g for g, c in zip(feature_genes, coef_vec) if c != 0.0
)
print(f"LASSO nonzero genes : {len(lasso_nonzero):,} of {len(feature_genes):,}")

# ---------------------------------------------------------------------------
# 3.  LOAD BEATAML — genes as rows, filter to signature genes only
#     Columns: stable_id | display_label | description | biotype | samples...
#     We index on the first column (stable_id = Ensembl ID), but genes in
#     signature_feature_genes.txt are HGNC symbols.  The second column
#     (display_label) holds the gene symbol — use that to subset.
# ---------------------------------------------------------------------------
print("\nLoading BeatAML signature gene rows ...")

# Read just enough to understand column structure, then filter on load
# Strategy: read with index_col=0, then use display_label (col index 0 after
# setting stable_id as index) to find signature genes.
beataml_raw = pd.read_csv(
    BEATAML_FILE,
    sep="\t",
    index_col=0,          # stable_id becomes row index
    low_memory=False
)
print(f"  Full BeatAML shape : {beataml_raw.shape}  (genes x samples+annotation)")

# Columns 0,1,2 are annotation (display_label, description, biotype)
# Sample expression data starts at column 3
annot_cols  = list(beataml_raw.columns[:3])    # display_label, description, biotype
sample_cols = list(beataml_raw.columns[3:])

print(f"  Annotation columns : {annot_cols}")
print(f"  Sample columns     : {len(sample_cols):,} samples")

# Filter rows where display_label (gene symbol) is in signature
beataml_raw["gene_symbol"] = beataml_raw[annot_cols[0]]   # display_label
beataml_sig = beataml_raw[beataml_raw["gene_symbol"].isin(feature_set)].copy()
beataml_sig = beataml_sig.set_index("gene_symbol")[sample_cols]

# Handle duplicate gene symbols (keep first occurrence)
beataml_sig = beataml_sig[~beataml_sig.index.duplicated(keep="first")]
beataml_sig = beataml_sig.apply(pd.to_numeric, errors="coerce")

print(f"  Signature genes found in BeatAML : {beataml_sig.shape[0]:,}")
print(f"  BeatAML sample count             : {beataml_sig.shape[1]:,}")

# ---------------------------------------------------------------------------
# 4.  LOAD EISFELD — genes as rows (GeneID column)
#     First column is gene ID (gene symbols in this file).
# ---------------------------------------------------------------------------
print("\nLoading Eisfeld signature gene rows ...")

eisfeld_raw = pd.read_csv(
    EISFELD_FILE,
    sep="\t",
    index_col=0,           # GeneID becomes row index
    low_memory=False
)
print(f"  Full Eisfeld shape : {eisfeld_raw.shape}  (genes x samples)")

# Filter to signature genes
eisfeld_sig = eisfeld_raw[eisfeld_raw.index.isin(feature_set)].copy()
eisfeld_sig = eisfeld_sig[~eisfeld_sig.index.duplicated(keep="first")]
eisfeld_sig = eisfeld_sig.apply(pd.to_numeric, errors="coerce")

print(f"  Signature genes found in Eisfeld : {eisfeld_sig.shape[0]:,}")
print(f"  Eisfeld sample count             : {eisfeld_sig.shape[1]:,}")

# ---------------------------------------------------------------------------
# 5.  COMPUTE PER-GENE MEAN AND VARIANCE
# ---------------------------------------------------------------------------
print("\nComputing per-gene statistics ...")

# Align to common genes found in BOTH datasets
common_genes = sorted(set(beataml_sig.index) & set(eisfeld_sig.index))
print(f"  Genes in both datasets     : {len(common_genes):,}")
print(f"  Only in BeatAML            : {len(set(beataml_sig.index) - set(eisfeld_sig.index)):,}")
print(f"  Only in Eisfeld            : {len(set(eisfeld_sig.index) - set(beataml_sig.index)):,}")

ba_sub  = beataml_sig.loc[common_genes]
ei_sub  = eisfeld_sig.loc[common_genes]

stats = pd.DataFrame({
    "gene"         : common_genes,
    "ba_mean"      : ba_sub.mean(axis=1).values,
    "ei_mean"      : ei_sub.mean(axis=1).values,
    "ba_var"       : ba_sub.var(axis=1).values,
    "ei_var"       : ei_sub.var(axis=1).values,
    "lasso_nonzero": [g in lasso_nonzero for g in common_genes],
}).set_index("gene")

# ---------------------------------------------------------------------------
# 6.  PRINT SUMMARY STATISTICS
# ---------------------------------------------------------------------------
print("\n" + "=" * 65)
print("SUMMARY STATISTICS")
print("=" * 65)

rho, pval = spearmanr(stats["ba_mean"], stats["ei_mean"])
print(f"  Spearman r (per-gene means)    : {rho:.4f}  (p = {pval:.2e})")

mad = np.median(np.abs(stats["ei_mean"] - stats["ba_mean"]))
print(f"  Median absolute diff in means  : {mad:.4f}")

pct_higher = 100 * (stats["ei_mean"] > stats["ba_mean"]).mean()
print(f"  % genes where Eisfeld > BeatAML: {pct_higher:.1f}%")

# Additional context
print(f"\n  BeatAML mean (across sig genes): {stats['ba_mean'].mean():.4f}")
print(f"  Eisfeld mean (across sig genes): {stats['ei_mean'].mean():.4f}")
print(f"  Global mean shift (Ei - BA)    : {(stats['ei_mean'] - stats['ba_mean']).mean():.4f}")
print(f"  BeatAML median               : {stats['ba_mean'].median():.4f}")
print(f"  Eisfeld median               : {stats['ei_mean'].median():.4f}")

# Variance comparison
print(f"\n  BeatAML mean variance          : {stats['ba_var'].mean():.4f}")
print(f"  Eisfeld mean variance          : {stats['ei_var'].mean():.4f}")

# Check for large-positive-offset signature — the core question
print(f"\n  Among LASSO nonzero genes ({stats['lasso_nonzero'].sum()} common with both datasets):")
lasso_stats = stats[stats["lasso_nonzero"]]
if len(lasso_stats) > 0:
    rho_l, pval_l = spearmanr(lasso_stats["ba_mean"], lasso_stats["ei_mean"])
    mad_l = np.median(np.abs(lasso_stats["ei_mean"] - lasso_stats["ba_mean"]))
    pct_l = 100 * (lasso_stats["ei_mean"] > lasso_stats["ba_mean"]).mean()
    print(f"    Spearman r                   : {rho_l:.4f}  (p = {pval_l:.2e})")
    print(f"    Median absolute diff in means: {mad_l:.4f}")
    print(f"    % genes Eisfeld > BeatAML    : {pct_l:.1f}%")
print("=" * 65)

# ---------------------------------------------------------------------------
# 7.  PLOT A — MEAN SCATTER
#     BeatAML mean (x) vs Eisfeld mean (y) per gene.
#     Grey = LASSO zero, orange = LASSO nonzero.
#     Identity line (y = x) in black dashed.
# ---------------------------------------------------------------------------
print("\nGenerating plots ...")

fig, ax = plt.subplots(figsize=(7, 7))

# Split into zero-coef and nonzero-coef genes
zero_mask    = ~stats["lasso_nonzero"]
nonzero_mask =  stats["lasso_nonzero"]

ax.scatter(
    stats.loc[zero_mask,    "ba_mean"],
    stats.loc[zero_mask,    "ei_mean"],
    s=8, alpha=0.35, color="#999999", linewidths=0, label="LASSO coef = 0"
)
ax.scatter(
    stats.loc[nonzero_mask, "ba_mean"],
    stats.loc[nonzero_mask, "ei_mean"],
    s=22, alpha=0.80, color="#E07B39", linewidths=0.4,
    edgecolors="#7A3A0A", label="LASSO coef ≠ 0"
)

# Identity line over data range
xmin = min(stats["ba_mean"].min(), stats["ei_mean"].min())
xmax = max(stats["ba_mean"].max(), stats["ei_mean"].max())
ax.plot([xmin, xmax], [xmin, xmax], "k--", linewidth=1.0, zorder=3, label="Identity (y = x)")

ax.set_xlabel("BeatAML per-gene mean expression", fontsize=12)
ax.set_ylabel("Eisfeld per-gene mean expression", fontsize=12)
ax.set_title(
    f"Per-gene mean expression: BeatAML vs Eisfeld\n"
    f"Signature genes (n={len(common_genes):,} common)  |  Spearman r = {rho:.3f}",
    fontsize=11
)
ax.legend(fontsize=9, framealpha=0.7)

# Annotate global offset
ax.text(
    0.04, 0.97,
    f"Median |Ei - BA| = {mad:.2f}\n% Eisfeld > BeatAML = {pct_higher:.0f}%",
    transform=ax.transAxes, fontsize=8.5, va="top",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9)
)

plt.tight_layout()
out_scatter = FIGURES_DIR / "norm_diag_mean_scatter.png"
fig.savefig(out_scatter, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_scatter}")

# ---------------------------------------------------------------------------
# 8.  PLOT B — GLOBAL DISTRIBUTION OVERLAY
#     Density of ALL expression values across all signature genes and all
#     samples, BeatAML vs Eisfeld.
#     This makes the global shift (or scale difference) immediately visible.
# ---------------------------------------------------------------------------

# Flatten to 1D arrays (cap at 500k points per dataset to keep KDE fast)
def _sample_flat(df, cap=500_000, seed=42):
    vals = df.values.ravel().astype(float)
    vals = vals[np.isfinite(vals)]
    if len(vals) > cap:
        rng = np.random.default_rng(seed)
        vals = rng.choice(vals, size=cap, replace=False)
    return vals

ba_flat = _sample_flat(ba_sub)
ei_flat = _sample_flat(ei_sub)

print(f"  BeatAML values for density (sample): {len(ba_flat):,}  "
      f"range [{ba_flat.min():.2f}, {ba_flat.max():.2f}]")
print(f"  Eisfeld values for density (sample): {len(ei_flat):,}  "
      f"range [{ei_flat.min():.2f}, {ei_flat.max():.2f}]")

def _kde_curve(vals, n_points=500):
    """Return x, y arrays for a KDE curve over vals."""
    lo, hi = np.percentile(vals, 0.1), np.percentile(vals, 99.9)
    xs = np.linspace(lo, hi, n_points)
    kde = gaussian_kde(vals, bw_method="scott")
    return xs, kde(xs)

ba_xs, ba_ys = _kde_curve(ba_flat)
ei_xs, ei_ys = _kde_curve(ei_flat)

fig, ax = plt.subplots(figsize=(8, 5))

ax.fill_between(ba_xs, ba_ys, alpha=0.25, color="#3A7EBF")
ax.plot(ba_xs, ba_ys, color="#3A7EBF", linewidth=2.0,
        label=f"BeatAML  (N={beataml_sig.shape[1]:,} samples, median={np.median(ba_flat):.2f})")

ax.fill_between(ei_xs, ei_ys, alpha=0.25, color="#E07B39")
ax.plot(ei_xs, ei_ys, color="#E07B39", linewidth=2.0,
        label=f"Eisfeld  (N={eisfeld_sig.shape[1]:,} samples, median={np.median(ei_flat):.2f})")

# Vertical lines at medians
ax.axvline(np.median(ba_flat), color="#3A7EBF", linestyle=":", linewidth=1.3)
ax.axvline(np.median(ei_flat), color="#E07B39", linestyle=":", linewidth=1.3)

ax.set_xlabel("Expression value (all signature genes, all samples pooled)", fontsize=11)
ax.set_ylabel("Density", fontsize=11)
ax.set_title(
    f"Expression distribution: BeatAML vs Eisfeld\n"
    f"Restricted to {len(common_genes):,} signature genes in common",
    fontsize=11
)
ax.legend(fontsize=9, framealpha=0.8)

# Annotate global shift
shift = np.median(ei_flat) - np.median(ba_flat)
ax.text(
    0.97, 0.97,
    f"Median shift (Eisfeld - BeatAML) = {shift:+.2f}",
    transform=ax.transAxes, fontsize=8.5, ha="right", va="top",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9)
)

plt.tight_layout()
out_density = FIGURES_DIR / "norm_diag_distributions.png"
fig.savefig(out_density, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_density}")

print("\nDone.")
