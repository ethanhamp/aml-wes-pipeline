"""
apply_signature_logcpm.py
=========================
Re-score Eisfeld ARHG-mutant patients using log2(CPM+1) normalization of raw
counts, matching BeatAML's training normalization more closely than the z-score
workaround.

Why this is better than z-score normalization:
  Z-scoring removes all cross-cohort absolute level information — it only
  preserves relative ordering within the Eisfeld cohort. Log2(CPM+1) puts
  the Eisfeld data into the same distributional space as BeatAML training
  data, so the saved scaler and model see values similar to what they were
  trained on. No information is discarded.

Why log2(CPM+1) not log2(CPM):
  BeatAML uses log2(CPM) which can go negative for lowly expressed genes.
  For zero-count genes, log2(0) is undefined. Adding +1 to CPM avoids -inf
  for zero-count genes while staying in the same range as BeatAML for all
  well-expressed genes (which make up the signature).

Steps:
  1. Load raw_counts.tsv (Ensembl IDs × 25 samples)
  2. Map Ensembl ID → HGNC symbol via gene_id_to_name.tsv
     - Duplicate symbols: keep the row with highest mean count
  3. Compute log2(CPM+1) per sample
  4. Align to signature_feature_genes.txt (fill missing with 0)
  5. Apply saved BeatAML scaler + LASSO model directly (no z-score step)
  6. Save applied_ras_scores_logcpm.csv
  7. Print comparison with previous z-score scores
"""

import os
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
COUNTS     = BASE / "raw_counts.tsv"
ID_MAP     = BASE / "gene_id_to_name.tsv"
FEAT_GENES = BASE / "signature_feature_genes.txt"
SCALER_PKL = BASE / "ras_scaler.pkl"
MODEL_PKL  = BASE / "ras_lasso_model.pkl"
PREV_SCORES = BASE / "applied_ras_scores.csv"
OUT        = BASE / "applied_ras_scores_logcpm.csv"

# ---------------------------------------------------------------------------
# Load model artifacts
# ---------------------------------------------------------------------------
with open(SCALER_PKL, "rb") as f:
    scaler = pickle.load(f)
with open(MODEL_PKL, "rb") as f:
    model = pickle.load(f)
with open(FEAT_GENES, "r") as f:
    feature_genes = [l.strip() for l in f if l.strip()]

print(f"Model loaded. Feature genes: {len(feature_genes):,}")

# ---------------------------------------------------------------------------
# Load raw counts (Ensembl IDs as rows, samples as columns)
# ---------------------------------------------------------------------------
# Header has literal \t text instead of real tabs — parse manually
with open(COUNTS) as f:
    header_line = f.readline().rstrip("\n")
    col_names = [c.strip() for c in header_line.split("\\t")]  # literal \t

counts = pd.read_csv(COUNTS, sep="\t", skiprows=1, header=None)
counts.columns = col_names
counts = counts.set_index(col_names[0])
print(f"Raw counts shape: {counts.shape}  (genes × samples)")

# ---------------------------------------------------------------------------
# Map Ensembl IDs → HGNC symbols
# ---------------------------------------------------------------------------
id_map = pd.read_csv(ID_MAP, sep="\t", header=None, names=["ensembl", "symbol"])
id_map = id_map.set_index("ensembl")["symbol"]

counts.index = counts.index.map(lambda x: id_map.get(x, None))
counts = counts[counts.index.notna()]                  # drop unmapped
counts.index = counts.index.astype(str)

# Duplicate gene symbols: keep the row with the highest mean count
counts["_mean"] = counts.mean(axis=1)
counts = counts.sort_values("_mean", ascending=False)
counts = counts[~counts.index.duplicated(keep="first")]
counts = counts.drop(columns=["_mean"])
print(f"After ID mapping + dedup: {counts.shape}  (genes × samples)")

# ---------------------------------------------------------------------------
# log2(CPM + 1) normalization
# Library size per sample = total raw counts in that sample
# CPM = count / library_size * 1e6
# log2_cpm1 = log2(CPM + 1)
# ---------------------------------------------------------------------------
library_sizes = counts.sum(axis=0)                     # per-sample totals
cpm = counts.divide(library_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

print(f"\nLog2(CPM+1) stats across all genes/samples:")
vals = log_cpm.values.flatten()
print(f"  min={vals.min():.2f}  mean={vals.mean():.2f}  max={vals.max():.2f}")

# Compare to BeatAML training space (scaler mean = BeatAML per-gene means)
print(f"\nBeatAML training space (scaler.mean_):")
print(f"  min={scaler.mean_.min():.2f}  mean={scaler.mean_.mean():.2f}  max={scaler.mean_.max():.2f}")

# ---------------------------------------------------------------------------
# Align to feature genes (samples as rows)
# ---------------------------------------------------------------------------
expr = log_cpm.T                                       # samples × genes
expr.index = expr.index.str.replace("_blast", "", regex=False)

feature_set = set(feature_genes)
present = len(feature_set & set(expr.columns))
missing = len(feature_set - set(expr.columns))
print(f"\nFeature gene alignment:")
print(f"  Found : {present:,} / {len(feature_genes):,}")
print(f"  Missing (fill 0): {missing:,}")

# Fill missing genes with the BeatAML training mean for that gene (scaler.mean_).
# This makes missing genes contribute 0 after scaler.transform — neutral, not deflating.
# Filling with 0 would map to (0 - mean)/std ≈ -3 to -4, which severely deflates scores.
gene_to_idx  = {g: i for i, g in enumerate(feature_genes)}
missing_fill = {g: scaler.mean_[gene_to_idx[g]]
                for g in feature_genes if g not in expr.columns}

expr_aligned = expr.reindex(columns=feature_genes)
for gene, val in missing_fill.items():
    expr_aligned[gene] = val

# ---------------------------------------------------------------------------
# Score: apply BeatAML scaler + LASSO model directly (no z-score step)
# ---------------------------------------------------------------------------
X_scaled = scaler.transform(expr_aligned.values.astype(float))
scores   = model.predict_proba(X_scaled)[:, 1]

results = pd.DataFrame({"sample_id": expr_aligned.index, "RAS_score_logcpm": scores})
results.to_csv(OUT, index=False)
print(f"\nScores saved: {OUT}")

# ---------------------------------------------------------------------------
# Compare to previous z-score scores
# ---------------------------------------------------------------------------
prev = pd.read_csv(PREV_SCORES).rename(
    columns={"RAS_score": "RAS_score_zscore", "sample_id": "sample_id"}
)
comparison = results.merge(prev, on="sample_id").sort_values("RAS_score_logcpm", ascending=False)
comparison["delta"] = comparison["RAS_score_logcpm"] - comparison["RAS_score_zscore"]

print("\n" + "=" * 65)
print("COMPARISON: log-CPM vs z-score normalization")
print("=" * 65)
print(f"{'Sample':<15} {'log-CPM':>10} {'z-score':>10} {'delta':>8}")
print("-" * 45)
for _, row in comparison.iterrows():
    print(f"{row['sample_id']:<15} {row['RAS_score_logcpm']:>10.3f} {row['RAS_score_zscore']:>10.3f} {row['delta']:>+8.3f}")

print(f"\nMean log-CPM score : {scores.mean():.3f}")
print(f"Mean z-score score : {prev['RAS_score_zscore'].mean():.3f}")
print(f"Samples >0.5 (log-CPM) : {(scores > 0.5).sum()} / {len(scores)}")
print(f"Samples >0.5 (z-score) : {(prev['RAS_score_zscore'] > 0.5).sum()} / {len(prev)}")
print("=" * 65)
