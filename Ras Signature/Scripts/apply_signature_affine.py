"""
apply_signature_affine.py
=========================
Apply the pre-trained BeatAML RAS pathway signature to the Eisfeld expression
dataset using gene-wise affine correction to resolve the normalization mismatch
between the two cohorts.

MOTIVATION
----------
The original apply_signature.py uses either raw Eisfeld values (biased low
because BeatAML scaler expects log-normalized CPM ~N(4, ...) while Eisfeld uses
a different normalization space) or per-gene z-scores within Eisfeld alone (which
is limited by the 25-sample Eisfeld cohort giving noisy per-gene SD estimates).

This script uses a more principled approach: gene-wise affine correction.
For each gene shared between BeatAML (707 samples) and Eisfeld (25 samples):

    corrected_eisfeld = (eisfeld_value - eisfeld_gene_mean) / eisfeld_gene_sd
                        * beataml_gene_sd + beataml_gene_mean

This re-centers and rescales each Eisfeld gene's distribution to match the
BeatAML reference mean and variance — essentially a gene-level z-score into
BeatAML space. Because BeatAML's per-gene statistics are estimated from 707
samples, they are far more stable than the within-Eisfeld 25-sample estimates.

EDGE CASES
----------
- eisfeld_gene_sd == 0 (gene uniformly expressed across 25 samples):
    Cannot z-score; set corrected value = beataml_gene_mean for all samples
    (imputes the gene at the BeatAML population mean, equivalent to z=0).
- Gene in BeatAML but missing from Eisfeld:
    Filled with 0.0 (same zero-fill as apply_signature.py) AFTER affine
    correction of the present genes. The scaler then handles these as usual.
- Gene in Eisfeld but not in BeatAML (can't compute reference stats):
    Treated as zero-fill after alignment to the feature gene list.

OUTPUTS
-------
  applied_ras_scores_affine.csv    — per-sample RAS scores after affine correction
  Figures/norm_diag_affine_comparison.png — before/after diagnostic figure
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
matplotlib.use("Agg")   # headless — no display needed
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import spearmanr

BASE_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
BEATAML_FILE = Path(os.environ.get("BEATAML_DIR", str(Path(os.environ.get("ARHG_BASE_DIR", "..")).parent / "BeatAML"))) / "beataml_waves1to4_norm_exp_dbgap.txt"
EISFELD_FILE  = BASE_DIR / "eisfeld_expression.tsv"
MODEL_PKL     = BASE_DIR / "ras_lasso_model.pkl"
SCALER_PKL    = BASE_DIR / "ras_scaler.pkl"
FEAT_GENES    = BASE_DIR / "signature_feature_genes.txt"
ORIG_SCORES   = BASE_DIR / "applied_ras_scores.csv"
OUTPUT_FILE   = BASE_DIR / "applied_ras_scores_affine.csv"
FIGURES_DIR   = BASE_DIR / "Figures"
FIGURES_DIR.mkdir(exist_ok=True)

print("=" * 65)
print("apply_signature_affine.py — gene-wise affine correction")
print("=" * 65)

# ---------------------------------------------------------------------------
# 1.  LOAD FEATURE GENE LIST AND MODEL ARTIFACTS
# ---------------------------------------------------------------------------
print("\nLoading model artifacts ...")

with open(FEAT_GENES) as fh:
    feature_genes = [line.strip() for line in fh if line.strip()]
feature_set = set(feature_genes)
print(f"  Feature genes      : {len(feature_genes):,}")

with open(SCALER_PKL, "rb") as fh:
    scaler = pickle.load(fh)
print(f"  Scaler loaded      : {SCALER_PKL.name}")

with open(MODEL_PKL, "rb") as fh:
    model = pickle.load(fh)
print(f"  Model loaded       : {MODEL_PKL.name}")

# ---------------------------------------------------------------------------
# 2.  LOAD BEATAML — compute per-gene reference statistics
#     (replicates normalization_diagnostics.py section 3 exactly)
#     Columns: stable_id (index) | display_label | description | biotype | samples...
# ---------------------------------------------------------------------------
print("\nLoading BeatAML expression (reference) ...")

beataml_raw = pd.read_csv(
    BEATAML_FILE,
    sep="\t",
    index_col=0,
    low_memory=False
)
print(f"  Full BeatAML shape : {beataml_raw.shape}  (genes x samples+annotation)")

annot_cols  = list(beataml_raw.columns[:3])    # display_label, description, biotype
sample_cols = list(beataml_raw.columns[3:])

# Filter to signature genes using the display_label (gene symbol) column
beataml_raw["gene_symbol"] = beataml_raw[annot_cols[0]]
beataml_sig = beataml_raw[beataml_raw["gene_symbol"].isin(feature_set)].copy()
beataml_sig = beataml_sig.set_index("gene_symbol")[sample_cols]
beataml_sig = beataml_sig[~beataml_sig.index.duplicated(keep="first")]
beataml_sig = beataml_sig.apply(pd.to_numeric, errors="coerce")

print(f"  Signature genes found in BeatAML : {beataml_sig.shape[0]:,}")
print(f"  BeatAML sample count             : {beataml_sig.shape[1]:,}")

# Compute reference statistics (per-gene, across all BeatAML samples)
ba_mean = beataml_sig.mean(axis=1)    # Series indexed by gene symbol
ba_std  = beataml_sig.std(axis=1, ddof=1)

# ---------------------------------------------------------------------------
# 3.  LOAD EISFELD — orient to samples x genes
#     (replicates normalization_diagnostics.py section 4)
# ---------------------------------------------------------------------------
print("\nLoading Eisfeld expression ...")

eisfeld_raw = pd.read_csv(
    EISFELD_FILE,
    sep="\t",
    index_col=0,
    low_memory=False
)
print(f"  Full Eisfeld shape : {eisfeld_raw.shape}  (genes x samples)")

# Auto-detect orientation using gene matches (same logic as apply_signature.py)
row_matches = len(feature_set & set(eisfeld_raw.index.astype(str)))
col_matches = len(feature_set & set(eisfeld_raw.columns.astype(str)))

if row_matches >= col_matches and row_matches > 0:
    # Genes are rows — transpose to samples x genes
    expr_mat = eisfeld_raw.T.copy()
    print(f"  Orientation: genes as rows — transposed to samples x genes")
else:
    expr_mat = eisfeld_raw.copy()
    print(f"  Orientation: samples as rows — used as-is")

expr_mat.index.name   = "sample_id"
expr_mat.columns.name = None
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce")

print(f"  Oriented shape     : {expr_mat.shape}  (samples x genes)")
print(f"  Eisfeld sample IDs : {list(expr_mat.index)}")

# ---------------------------------------------------------------------------
# 4.  IDENTIFY GENE OVERLAP
# ---------------------------------------------------------------------------
eisfeld_genes = set(expr_mat.columns)
common_genes  = sorted(feature_set & eisfeld_genes & set(ba_mean.index))
missing_genes = feature_set - eisfeld_genes  # in signature but not in Eisfeld

print(f"\nGene coverage:")
print(f"  Common (signature + BeatAML + Eisfeld) : {len(common_genes):,}")
print(f"  In signature but missing from Eisfeld  : {len(missing_genes):,} (zero-filled)")
print(f"  In signature but missing from BeatAML  : {len(feature_set - set(ba_mean.index)):,} (no ref stats available)")

# ---------------------------------------------------------------------------
# 5.  APPLY GENE-WISE AFFINE CORRECTION
#     For each common gene g:
#       corrected[g] = (eisfeld[g] - ei_mean[g]) / ei_std[g] * ba_std[g] + ba_mean[g]
#
#     Edge case: ei_std[g] == 0 (gene flat across 25 Eisfeld samples)
#       We cannot compute a meaningful z-score.  Resolution: set corrected[g]
#       = ba_mean[g] for all patients (equivalent to z=0 — imputes at the
#       BeatAML population mean). This is conservative and avoids divide-by-zero.
#
#     Note on SD scaling: BeatAML's ba_std is estimated from 707 samples and
#     is reliable. Eisfeld's ei_std from 25 samples is noisy, but we only use
#     it to divide OUT the Eisfeld spread; the ba_std then imposes the correct
#     variance. The correction degrades gracefully for high-variance ei_std
#     outliers (they get rescaled to ba_std regardless).
# ---------------------------------------------------------------------------
print("\nApplying gene-wise affine correction ...")

# Work on the common-gene subset only
expr_common = expr_mat[common_genes].copy()    # shape: (25 samples x n_common_genes)

# Per-Eisfeld-sample gene statistics (axis=0 = across samples per gene)
ei_mean = expr_common.mean(axis=0)    # Series indexed by gene
ei_std  = expr_common.std(axis=0, ddof=1)

# Identify zero-variance genes in Eisfeld
zero_var_mask = ei_std == 0.0
n_zero_var    = zero_var_mask.sum()
print(f"  Genes with ei_std == 0 (uniform across 25 samples) : {n_zero_var:,}")
print(f"  These will be imputed at BeatAML mean (z=0)")

# Compute corrected values gene by gene using numpy broadcasting
# expr_common : (n_samples, n_genes) numpy array
expr_np  = expr_common.values.astype(float)           # (25, n_common)
ei_m_np  = ei_mean.values                             # (n_common,)
ei_s_np  = ei_std.values                              # (n_common,)
ba_m_np  = ba_mean.loc[common_genes].values           # (n_common,)
ba_s_np  = ba_std.loc[common_genes].values            # (n_common,)

# Replace zero SD with 1.0 temporarily to avoid divide-by-zero
# (result will be overwritten for those genes below)
ei_s_safe = np.where(ei_s_np == 0.0, 1.0, ei_s_np)

# Affine transform: vectorized across all samples and genes simultaneously
expr_corrected_np = (expr_np - ei_m_np) / ei_s_safe * ba_s_np + ba_m_np  # (25, n_common)

# For zero-variance Eisfeld genes: replace every sample's corrected value with
# the BeatAML mean (broadcasting: ba_m_np[mask] is a scalar per gene)
if n_zero_var > 0:
    zero_idx = np.where(zero_var_mask.values)[0]
    for idx in zero_idx:
        expr_corrected_np[:, idx] = ba_m_np[idx]

expr_corrected = pd.DataFrame(
    expr_corrected_np,
    index=expr_common.index,
    columns=common_genes
)

print(f"\n  Before correction (Eisfeld common genes):")
print(f"    Min={expr_np.min():.3f}  Mean={expr_np.mean():.3f}  Max={expr_np.max():.3f}")
print(f"  After correction (mapped to BeatAML space):")
print(f"    Min={expr_corrected_np.min():.3f}  Mean={expr_corrected_np.mean():.3f}  Max={expr_corrected_np.max():.3f}")
print(f"  BeatAML reference (same genes):")
ba_sub = beataml_sig.loc[common_genes]
print(f"    Min={ba_sub.values.min():.3f}  Mean={ba_sub.values.mean():.3f}  Max={ba_sub.values.max():.3f}")

# ---------------------------------------------------------------------------
# 6.  RECONSTRUCT FULL FEATURE MATRIX (corrected common + zero-fill missing)
#     Must align to the exact feature_genes ordering expected by scaler/model.
# ---------------------------------------------------------------------------
print("\nAligning corrected features to model input space ...")

# Start with corrected common genes, then reindex to full feature list
# Genes in feature_genes but not in common_genes will be NaN → filled with 0
expr_aligned = expr_corrected.reindex(columns=feature_genes, fill_value=0.0)
n_nan = expr_aligned.isna().sum().sum()
if n_nan > 0:
    expr_aligned = expr_aligned.fillna(0.0)

n_zero_filled = len(missing_genes & feature_set)
print(f"  Feature matrix shape (samples x genes) : {expr_aligned.shape}")
print(f"  Genes zero-filled (not in Eisfeld)      : {n_zero_filled:,}")

# ---------------------------------------------------------------------------
# 7.  SCORE WITH SCALER + MODEL
# ---------------------------------------------------------------------------
print("\nScoring samples with BeatAML scaler and LASSO model ...")

X_corrected = expr_aligned.values.astype(float)
X_scaled    = scaler.transform(X_corrected)              # apply BeatAML train scaler
ras_scores_corrected = model.predict_proba(X_scaled)[:, 1]   # P(RAS-mutant)

print(f"  Samples scored : {len(ras_scores_corrected):,}")
print(f"  Score range    : {ras_scores_corrected.min():.4f} - {ras_scores_corrected.max():.4f}")
print(f"  Score mean     : {ras_scores_corrected.mean():.4f}  (median: {np.median(ras_scores_corrected):.4f})")

# ---------------------------------------------------------------------------
# 8.  SAVE CORRECTED SCORES
# ---------------------------------------------------------------------------
results_df = pd.DataFrame({
    "sample_id": expr_aligned.index.astype(str),
    "RAS_score_affine": ras_scores_corrected,
})
results_df.to_csv(OUTPUT_FILE, index=False)
print(f"\nSaved corrected scores: {OUTPUT_FILE}")

# ---------------------------------------------------------------------------
# 9.  LOAD ORIGINAL SCORES AND COMPUTE COMPARISON STATS
# ---------------------------------------------------------------------------
print("\nLoading original scores for comparison ...")

orig_df = pd.read_csv(ORIG_SCORES)
orig_df = orig_df.set_index("sample_id")

# Align on sample order (use the scored samples from this run as canonical order)
sample_order    = list(expr_aligned.index.astype(str))
scores_original = orig_df.loc[sample_order, "RAS_score"].values
scores_corrected = ras_scores_corrected   # already in same order

print("\n" + "=" * 65)
print("COMPARISON STATISTICS")
print("=" * 65)

mean_orig = scores_original.mean()
mean_corr = scores_corrected.mean()
print(f"  Mean score — original  : {mean_orig:.4f}")
print(f"  Mean score — corrected : {mean_corr:.4f}")
print(f"  Mean shift (corr - orig): {mean_corr - mean_orig:+.4f}")

median_orig = np.median(scores_original)
median_corr = np.median(scores_corrected)
print(f"\n  Median score — original  : {median_orig:.4f}")
print(f"  Median score — corrected : {median_corr:.4f}")

rho, pval = spearmanr(scores_original, scores_corrected)
print(f"\n  Rank-order Spearman r (orig vs corrected): {rho:.4f}  (p = {pval:.3e})")
print(f"  Interpretation: {'high — rank ordering preserved' if rho > 0.8 else 'moderate — some rank reordering' if rho > 0.5 else 'low — substantial rank reordering'}")

THRESHOLD = 0.5
class_orig = scores_original  > THRESHOLD
class_corr = scores_corrected > THRESHOLD
n_changed  = (class_orig != class_corr).sum()
print(f"\n  Classification threshold : {THRESHOLD}")
print(f"  Original  high-RAS (>0.5) : {class_orig.sum():,} of {len(class_orig):,}  ({100*class_orig.mean():.1f}%)")
print(f"  Corrected high-RAS (>0.5) : {class_corr.sum():,} of {len(class_corr):,}  ({100*class_corr.mean():.1f}%)")
print(f"  Patients changing class    : {n_changed:,}")
if n_changed > 0:
    changed_ids = [s for s, o, c in zip(sample_order, class_orig, class_corr) if o != c]
    changed_dir = ["low->high" if not o and c else "high->low"
                   for o, c in zip(class_orig, class_corr) if o != c]
    for sid, d in zip(changed_ids, changed_dir):
        orig_s = scores_original[sample_order.index(sid)]
        corr_s = scores_corrected[sample_order.index(sid)]
        print(f"      {sid:20s}  {d}  (orig={orig_s:.3f} -> corr={corr_s:.3f})")

print("=" * 65)

# ---------------------------------------------------------------------------
# 10. LOAD BEATAML RAS-MUTANT SCORES FOR REFERENCE OVERLAY
#     ras_activation_scores.csv was produced by ras_signature.py and contains
#     per-sample BeatAML scores with a "RAS_mutant" label column.
# ---------------------------------------------------------------------------
beataml_scores_file = BASE_DIR / "ras_activation_scores.csv"
beataml_scores_df   = None
ras_mut_scores      = None
ras_wt_scores       = None

if beataml_scores_file.exists():
    beataml_scores_df = pd.read_csv(beataml_scores_file)
    print(f"\nBeatAML score file loaded: {beataml_scores_file.name}")
    print(f"  Columns: {list(beataml_scores_df.columns)}")

    # Identify the score column and the mutant label column
    score_col  = next((c for c in beataml_scores_df.columns
                       if "score" in c.lower()), None)
    label_col  = next((c for c in beataml_scores_df.columns
                       if "ras_mutant" in c.lower() or "label" in c.lower()
                          or "mutant" in c.lower()), None)

    if score_col and label_col:
        ras_mut_scores = beataml_scores_df.loc[
            beataml_scores_df[label_col] == 1, score_col].dropna().values
        ras_wt_scores  = beataml_scores_df.loc[
            beataml_scores_df[label_col] == 0, score_col].dropna().values
        print(f"  Score col: '{score_col}', Label col: '{label_col}'")
        print(f"  RAS-mutant samples : {len(ras_mut_scores):,}  (mean={ras_mut_scores.mean():.3f})")
        print(f"  RAS-WT samples     : {len(ras_wt_scores):,}  (mean={ras_wt_scores.mean():.3f})")
    else:
        print(f"  WARNING: Could not identify score or label columns — skipping BeatAML overlay")
else:
    print(f"\nNOTE: {beataml_scores_file.name} not found — BeatAML reference overlay will be omitted")

# ---------------------------------------------------------------------------
# 11. DIAGNOSTIC FIGURE
#
#     Left panel:  KDE / stripplot of score distributions
#                  - Original Eisfeld scores (orange)
#                  - Affine-corrected Eisfeld scores (teal)
#                  - BeatAML RAS-mutant scores (blue, if available)
#
#     Right panel: Scatter of original vs corrected score per patient
#                  - Points colored by classification change
#                  - Labeled with sample ID
#                  - Spearman r annotated
# ---------------------------------------------------------------------------
print("\nGenerating diagnostic figure ...")

from scipy.stats import gaussian_kde

def _kde_curve(vals, n_points=300):
    """Return x, y for a KDE over the provided array."""
    if len(vals) < 3:
        return np.array([]), np.array([])
    lo, hi = max(0.0, vals.min() - 0.05), min(1.0, vals.max() + 0.05)
    xs = np.linspace(lo, hi, n_points)
    kde = gaussian_kde(vals, bw_method="silverman")
    return xs, kde(xs)

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
fig.subplots_adjust(wspace=0.35)

# ---- Left panel: score distributions ----
ax = axes[0]

# Determine which groups to show and their colors/labels
groups = [
    (scores_original,  "#E07B39", "Eisfeld original\n(within-cohort z-score)", 2.0),
    (scores_corrected, "#2A9D8F", "Eisfeld affine-corrected\n(mapped to BeatAML space)", 2.0),
]
if ras_mut_scores is not None:
    groups.append((ras_mut_scores, "#3A7EBF",
                   f"BeatAML RAS-mutant\n(N={len(ras_mut_scores):,})", 1.3))
if ras_wt_scores is not None:
    groups.append((ras_wt_scores, "#AAAAAA",
                   f"BeatAML RAS-WT\n(N={len(ras_wt_scores):,})", 1.0))

for vals, color, label, lw in groups:
    xs, ys = _kde_curve(vals)
    if len(xs) == 0:
        continue
    ax.fill_between(xs, ys, alpha=0.15, color=color)
    ax.plot(xs, ys, color=color, linewidth=lw, label=f"{label}\n  mean={vals.mean():.3f}")

    # Rug marks at the bottom for Eisfeld groups (not BeatAML — too many)
    if "Eisfeld" in label or len(vals) <= 30:
        ax.plot(vals, np.full_like(vals, -0.03 * ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else -0.1),
                "|", color=color, alpha=0.7, markersize=8, markeredgewidth=1.5)

# Threshold line
ax.axvline(THRESHOLD, color="black", linestyle="--", linewidth=0.9,
           label=f"Threshold = {THRESHOLD}", zorder=5)

ax.set_xlabel("RAS activation score", fontsize=11)
ax.set_ylabel("Density", fontsize=11)
ax.set_title("Score distributions:\nbefore and after affine correction", fontsize=11, fontweight="bold")
ax.legend(fontsize=7.5, framealpha=0.8, loc="upper left")
ax.set_xlim(-0.05, 1.05)

# ---- Right panel: original vs corrected scatter ----
ax2 = axes[1]

# Color by classification status
no_change_same_low  = (~class_orig) & (~class_corr)   # stayed low
no_change_same_high =  class_orig   &  class_corr     # stayed high
changed_lo_to_hi    = (~class_orig) &  class_corr     # gained high
changed_hi_to_lo    =  class_orig   & (~class_corr)   # lost high

scatter_groups = [
    (no_change_same_low,  "#AAAAAA", "Stable low (both <0.5)",  50, 0.8),
    (no_change_same_high, "#2A9D8F", "Stable high (both >0.5)", 70, 0.9),
    (changed_lo_to_hi,    "#E63946", "Low -> High after correction", 90, 1.0),
    (changed_hi_to_lo,    "#457B9D", "High -> Low after correction", 90, 1.0),
]

for mask, color, label, size, alpha in scatter_groups:
    if mask.sum() == 0:
        continue
    ax2.scatter(
        scores_original[mask], scores_corrected[mask],
        s=size, color=color, alpha=alpha, edgecolors="white",
        linewidths=0.5, label=f"{label} (n={mask.sum()})", zorder=4
    )

# Label every point with sample ID
for i, sid in enumerate(sample_order):
    ax2.annotate(
        sid,
        xy=(scores_original[i], scores_corrected[i]),
        xytext=(3, 3),
        textcoords="offset points",
        fontsize=5.5,
        color="#333333",
        alpha=0.85,
    )

# Identity line (no change)
lim_min = min(scores_original.min(), scores_corrected.min()) - 0.03
lim_max = max(scores_original.max(), scores_corrected.max()) + 0.03
ax2.plot([lim_min, lim_max], [lim_min, lim_max],
         "k--", linewidth=0.8, alpha=0.5, zorder=2, label="Identity (no change)")

# Threshold crosshairs
ax2.axvline(THRESHOLD, color="grey", linestyle=":", linewidth=0.7, alpha=0.6)
ax2.axhline(THRESHOLD, color="grey", linestyle=":", linewidth=0.7, alpha=0.6)

ax2.set_xlabel("Original RAS score (within-Eisfeld z-score)", fontsize=11)
ax2.set_ylabel("Affine-corrected RAS score", fontsize=11)
ax2.set_title(
    f"Original vs affine-corrected scores\n"
    f"Spearman r = {rho:.3f}  |  {n_changed} patient(s) changed class",
    fontsize=11, fontweight="bold"
)
ax2.legend(fontsize=8, framealpha=0.8, loc="upper left")
ax2.set_xlim(lim_min, lim_max)
ax2.set_ylim(lim_min, lim_max)

# Annotate stats box
ax2.text(
    0.97, 0.05,
    f"Mean: {mean_orig:.3f} -> {mean_corr:.3f}\nSpearman r = {rho:.3f}",
    transform=ax2.transAxes, fontsize=8.5, ha="right", va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", alpha=0.9)
)

plt.suptitle(
    "RAS signature: gene-wise affine normalization correction\n"
    "(BeatAML reference statistics applied to Eisfeld expression)",
    fontsize=12, fontweight="bold", y=1.01
)

out_fig = FIGURES_DIR / "norm_diag_affine_comparison.png"
fig.savefig(out_fig, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_fig}")

# ---------------------------------------------------------------------------
# 12. FINAL SUMMARY
# ---------------------------------------------------------------------------
print("\n" + "=" * 65)
print("FINAL SUMMARY")
print("=" * 65)
print(f"  Method           : Gene-wise affine correction")
print(f"  Reference cohort : BeatAML ({beataml_sig.shape[1]:,} samples)")
print(f"  Target cohort    : Eisfeld ({expr_mat.shape[0]:,} samples)")
print(f"  Genes corrected  : {len(common_genes):,} (of {len(feature_genes):,} signature genes)")
print(f"  Genes zero-filled: {len(feature_genes) - len(common_genes):,}")
print(f"  Zero-var Eisfeld : {n_zero_var:,} genes imputed at BeatAML mean")
print()
print(f"  Score mean  : {mean_orig:.4f} (original) -> {mean_corr:.4f} (corrected)")
print(f"  Score median: {median_orig:.4f} (original) -> {median_corr:.4f} (corrected)")
print(f"  Spearman r  : {rho:.4f}  (rank-order preservation)")
print(f"  Class changes at threshold {THRESHOLD}: {n_changed}")
print()
print(f"  Output scores  : {OUTPUT_FILE}")
print(f"  Diagnostic fig : {out_fig}")
print("=" * 65)
print("Done.")
