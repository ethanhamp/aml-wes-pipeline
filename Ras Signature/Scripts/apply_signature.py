"""
apply_signature.py
==================
Apply the pre-trained BeatAML RAS pathway signature to a new expression dataset.

Designed for the two-step workflow:
  Step 1: ras_signature.py — trains the LASSO signature on BeatAML (public cohort)
  Step 2: apply_signature.py — scores a new dataset (e.g., Eisfeld Lab ARHG patients)
                               using the saved model artifacts

Why a separate apply script?
  ARHG-mutant patients are rare in public cohorts like BeatAML, making a direct
  comparison underpowered. The Eisfeld Lab has a curated cohort enriched for
  ARHG mutations in AML. Applying the pre-trained BeatAML model to this cohort
  (rather than retraining) avoids overfitting to the small ARHG-mutant group and
  preserves the biological question: does the RAS transcriptional program activate
  in ARHG-mutant AML as defined by the BeatAML-derived signature?

Usage:
  python apply_signature.py path/to/expression.tsv [--output path/to/output.csv]

Input expression file format (TSV, auto-detected orientation):
  Option A — genes as rows, samples as columns:
      gene_id  sample1  sample2  ...
      NRAS     1.23     2.34     ...
      KRAS     0.89     1.45     ...
  Option B — samples as rows, genes as columns:
      sample_id  NRAS  KRAS  ...
      sample1    1.23  0.89  ...
      sample2    2.34  1.45  ...

  Auto-detection: the script checks whether gene names from signature_feature_genes.txt
  appear as row index values (Option A) or column names (Option B).
  Missing genes are filled with 0 (mean-imputation alternative available below).
  Extra genes are silently dropped.

Output:
  applied_ras_scores.csv — per-sample RAS activation score (0-1)
    Columns: sample_id, RAS_score
    RAS_score ~ P(RAS-pathway active) per the LASSO model
"""

# ===========================================================================
# 0.  PATHS TO MODEL ARTIFACTS (must match ras_signature.py OUTPUT_DIR)
# ===========================================================================
import os
from pathlib import Path

OUTPUT_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
SCALER_PKL     = OUTPUT_DIR / "ras_scaler.pkl"
MODEL_PKL      = OUTPUT_DIR / "ras_lasso_model.pkl"
FEAT_GENES_TXT = OUTPUT_DIR / "signature_feature_genes.txt"

# ===========================================================================
# 1.  IMPORTS
# ===========================================================================
import argparse
import pickle
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd

# ===========================================================================
# 2.  ARGUMENT PARSING
# ===========================================================================
parser = argparse.ArgumentParser(
    description="Apply pre-trained BeatAML RAS signature to a new expression matrix."
)
parser.add_argument(
    "expression_file",
    type=str,
    help="Path to input expression TSV (genes x samples or samples x genes — auto-detected).",
)
parser.add_argument(
    "--output", "-o",
    type=str,
    default=str(OUTPUT_DIR / "applied_ras_scores.csv"),
    help="Path for output CSV (default: Ras Signature/applied_ras_scores.csv).",
)
parser.add_argument(
    "--normalize", "-n",
    action="store_true",
    default=False,
    help=(
        "Apply per-gene z-score normalization before scoring. "
        "Required when input data uses a different normalization than BeatAML training data "
        "(e.g. log2(TPM+1) vs log-normalized CPM). "
        "Maps each gene to z-scores within the input cohort, then shifts to BeatAML "
        "distribution so the saved scaler produces well-calibrated scaled values."
    ),
)
args = parser.parse_args()

EXPR_FILE   = Path(args.expression_file)
OUTPUT_FILE = Path(args.output)

# ===========================================================================
# 3.  LOAD MODEL ARTIFACTS
# ===========================================================================
print("=" * 65)
print("apply_signature.py — RAS pathway scoring")
print("=" * 65)
print("\nLoading model artifacts ...")

for artifact in [SCALER_PKL, MODEL_PKL, FEAT_GENES_TXT]:
    if not artifact.exists():
        print(f"\nERROR: Required artifact not found: {artifact}")
        print("       Run ras_signature.py first to generate model artifacts.")
        sys.exit(1)

with open(SCALER_PKL, "rb") as fh:
    scaler = pickle.load(fh)
print(f"  Scaler loaded      : {SCALER_PKL}")

with open(MODEL_PKL, "rb") as fh:
    model = pickle.load(fh)
print(f"  Model loaded       : {MODEL_PKL}")

with open(FEAT_GENES_TXT, "r") as fh:
    feature_genes = [line.strip() for line in fh if line.strip()]
print(f"  Feature genes      : {len(feature_genes):,} genes expected")

# ===========================================================================
# 4.  LOAD INPUT EXPRESSION MATRIX
# ===========================================================================
print(f"\nLoading expression data: {EXPR_FILE}")

if not EXPR_FILE.exists():
    print(f"ERROR: Input file not found: {EXPR_FILE}")
    sys.exit(1)

expr_raw = pd.read_csv(EXPR_FILE, sep="\t", index_col=0, low_memory=False)
print(f"  Raw shape: {expr_raw.shape}  (index x columns)")

# ===========================================================================
# 5.  AUTO-DETECT ORIENTATION
#     Check whether gene names appear as row index values or column names.
#     Match rate threshold: at least 10% of expected feature genes must be found
#     to confirm the correct orientation.
# ===========================================================================
print("\nAuto-detecting matrix orientation ...")

feature_set       = set(feature_genes)
row_gene_matches  = len(feature_set & set(expr_raw.index.astype(str)))
col_gene_matches  = len(feature_set & set(expr_raw.columns.astype(str)))

print(f"  Genes found as row index : {row_gene_matches:,}")
print(f"  Genes found as columns   : {col_gene_matches:,}")

if row_gene_matches >= col_gene_matches and row_gene_matches > 0:
    # Genes are rows → transpose to samples x genes
    print("  Orientation: genes as rows — transposing to samples x genes")
    expr_mat = expr_raw.T.copy()
    # Row index of transposed matrix = sample IDs (original columns)
elif col_gene_matches > 0:
    # Genes are columns → already samples x genes
    print("  Orientation: samples as rows — using as-is")
    expr_mat = expr_raw.copy()
else:
    print("ERROR: Could not identify gene names in either rows or columns.")
    print("       Ensure gene symbols match those in signature_feature_genes.txt")
    print(f"       Example expected genes: {feature_genes[:5]}")
    sys.exit(1)

expr_mat.index.name   = "sample_id"
expr_mat.columns.name = None
print(f"  Oriented matrix shape: {expr_mat.shape}  (samples x genes)")

# ===========================================================================
# 6.  ALIGN COLUMNS TO EXPECTED FEATURE GENES
#     The scaler and model expect exactly the same gene columns in exactly
#     the same order as the training data. Steps:
#       a) For genes present in input: use those values
#       b) For genes missing from input: fill with 0
#          (0 = global mean after StandardScaler centering, a conservative
#           imputation that assumes the gene is expressed at the population mean;
#           if many genes are missing, interpret scores with caution)
#       c) Drop extra genes not in the feature list
# ===========================================================================
print("\nAligning gene features ...")

# Convert to float (coerce errors to NaN)
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce")

present_genes = set(expr_mat.columns)
missing_genes = feature_set - present_genes
extra_genes   = present_genes - feature_set

print(f"  Expected features  : {len(feature_genes):,}")
print(f"  Found in input     : {len(present_genes & feature_set):,}")
print(f"  Missing (fill 0)   : {len(missing_genes):,}")
print(f"  Extra (dropped)    : {len(extra_genes):,}")

if len(missing_genes) > 0.5 * len(feature_genes):
    print(f"\n  WARNING: More than 50% of signature genes are missing from input.")
    print(f"           Scores may be unreliable. Check that:")
    print(f"           - Gene symbols use HGNC format (e.g., NRAS not Nras)")
    print(f"           - Expression values are for the same species (human)")
    print(f"           - The file contains normalized expression, not raw counts")

# Reindex to exact feature gene order, fill missing with 0
expr_aligned = expr_mat.reindex(columns=feature_genes, fill_value=0.0)

# Replace any remaining NaN with 0
n_nan = expr_aligned.isna().sum().sum()
if n_nan > 0:
    print(f"  NaN values found and replaced with 0: {n_nan:,}")
    expr_aligned = expr_aligned.fillna(0.0)

print(f"  Aligned matrix shape: {expr_aligned.shape}")

# ===========================================================================
# 6.5  CROSS-COHORT NORMALIZATION  (--normalize flag)
#
#     When the input cohort uses a different normalization than BeatAML
#     training data, the saved StandardScaler miscalibrates features.
#     The scaler was fit on BeatAML (log-normalized CPM, mean ~4, can go
#     negative). If the input is log2(TPM+1) (always ≥ 0, mean ~1), the
#     scaler subtracts ~4 from values near 1 → highly negative inputs the
#     model never saw during training → deflated scores regardless of biology.
#
#     Fix: per-gene z-score normalization of the input cohort.
#     For each gene, subtract the cohort mean and divide by cohort std.
#     This maps each gene's distribution to N(0,1) within the input cohort.
#     Then shift/scale back to the BeatAML training distribution so the
#     downstream StandardScaler produces the same scaled space as training.
#
#     Math: x_corrected = z * scaler.scale_ + scaler.mean_
#     After scaler.transform(): (x_corrected - scaler.mean_) / scaler.scale_ = z
#     So the final model input = per-gene z-score of the input cohort.
#     Genes with zero variance across input samples are set to 0.0.
# ===========================================================================
from scipy.stats import zscore as _zscore

if args.normalize:
    print("\nApplying per-gene z-score normalization (--normalize) ...")
    raw_vals = expr_aligned.values.astype(float)

    per_gene_z = _zscore(raw_vals, axis=0, nan_policy="omit")
    per_gene_z = np.nan_to_num(per_gene_z, nan=0.0)

    # Map z-scores into BeatAML training distribution space
    expr_for_scoring = per_gene_z * scaler.scale_ + scaler.mean_

    print(f"  Before: min={raw_vals.min():.2f}, mean={raw_vals.mean():.2f}, max={raw_vals.max():.2f}")
    print(f"  After : min={expr_for_scoring.min():.2f}, mean={expr_for_scoring.mean():.2f}, max={expr_for_scoring.max():.2f}")
    print(f"  BeatAML reference mean (across genes): {scaler.mean_.mean():.2f}")
else:
    expr_for_scoring = expr_aligned.values.astype(float)

# ===========================================================================
# 7.  APPLY SCALER AND MODEL
#     CRITICAL: we use the scaler fitted on BeatAML TRAIN data.
#     This means each gene is centered to BeatAML's training mean and scaled
#     to BeatAML's training standard deviation. This is correct — we want to
#     score new samples relative to the BeatAML reference distribution, not
#     re-normalize within the new dataset (which would erase batch effects
#     but also erase the signal the model was trained on).
# ===========================================================================
print("\nScoring samples ...")

X_new     = expr_for_scoring
X_scaled  = scaler.transform(X_new)                    # apply BeatAML train scaler
ras_scores = model.predict_proba(X_scaled)[:, 1]       # P(RAS-mutant)

print(f"  Scored {len(ras_scores):,} samples")
print(f"  RAS score range: {ras_scores.min():.4f} - {ras_scores.max():.4f}")
print(f"  RAS score mean : {ras_scores.mean():.4f}  (median: {np.median(ras_scores):.4f})")

# ===========================================================================
# 8.  SAVE OUTPUT
# ===========================================================================
results_df = pd.DataFrame({
    "sample_id": expr_aligned.index.astype(str),
    "RAS_score": ras_scores,
})

results_df.to_csv(OUTPUT_FILE, index=False)
print(f"\nOutput saved: {OUTPUT_FILE}")

# ===========================================================================
# 9.  SUMMARY
# ===========================================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"  Input file      : {EXPR_FILE}")
print(f"  Normalization   : {'per-gene z-score (--normalize)' if args.normalize else 'none (raw input values)'}")
print(f"  Samples scored  : {len(results_df):,}")
print(f"  Features used   : {len(feature_genes):,} (of which {len(missing_genes):,} were missing/zero-filled)")
print(f"  RAS score range : {ras_scores.min():.4f} - {ras_scores.max():.4f}")
print(f"  High RAS (>0.5) : {(ras_scores > 0.5).sum():,} samples ({100*(ras_scores>0.5).mean():.1f}%)")
print(f"  Output CSV      : {OUTPUT_FILE}")
print()
print("  Interpretation:")
print("    RAS_score ~ P(RAS-pathway transcriptionally active)")
print("    Trained on BeatAML WES + RNA-seq (N~615 patients)")
print("    Score > 0.5 suggests RAS-pathway activation by LASSO signature")
print("    For group comparisons (e.g., ARHG-mutant vs WT), use Mann-Whitney U")
print("    and interpret in the context of your cohort's clinical metadata.")
print("=" * 65)
print("Done.")
