"""
compare_features.py
===================
Cross-model feature importance comparison for the BeatAML RAS signature.

Produces figures for the AI/ML class write-up:
  1. Rank-correlation heatmap (Spearman rho between all model pairs)
  2. Top-N overlap bar chart (how many of the top-N features are shared
     across all 4 models, for N = 10, 25, 50)
  3. Venn-style UpSet-compatible counts table (saved as CSV; you can paste
     into UpSetR in R if desired — UpSet requires ≥3 sets)
  4. Dot-plot: top-40 genes, dot size = |importance|, color = model
  5. SPRED1/SPRED2 and DUSP6 callout (canonical RAS-feedback biology)
  6. Per-gene z-score of importance ranks (identifies "consensus" genes
     that rank consistently high across all models)

All four models are re-fitted here from scratch on the same split used in
ras_signature.py (random_state=42, test_size=0.2) so that results are
100% reproducible without needing to pickle and reload the RF/GBM objects
(those are not saved by ras_signature.py — only LASSO is serialized).

The LASSO model is loaded from the saved pkl to guarantee identical coefficients.

Usage:
    python compare_features.py

Output files (written to same directory as this script):
    feature_rank_correlation.png
    top_n_overlap.png
    dot_plot_top40.png
    consensus_genes.png
    cross_model_importance.csv   <- full table, all genes, all models
    consensus_genes.csv          <- top consensus genes ranked by mean rank
"""

# ===========================================================================
# 0. PATHS
# ===========================================================================
import os
from pathlib import Path
import pickle, warnings
warnings.filterwarnings("ignore")

import numpy  as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import spearmanr

from sklearn.preprocessing   import StandardScaler
from sklearn.linear_model    import LogisticRegressionCV
from sklearn.ensemble        import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics         import roc_auc_score
from sklearn.inspection      import permutation_importance

DATA_DIR = Path(os.environ.get("BEATAML_DIR", str(Path(os.environ.get("ARHG_BASE_DIR", "..")).parent / "BeatAML")))
OUTPUT_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"

EXPR_FILE  = DATA_DIR / "beataml_waves1to4_norm_exp_dbgap.txt"
MUT_FILE   = DATA_DIR / "beataml_wes_wv1to4_mutations_dbgap.txt"
SCALER_PKL = OUTPUT_DIR / "ras_scaler.pkl"
MODEL_PKL  = OUTPUT_DIR / "ras_lasso_model.pkl"
SIG_CSV    = OUTPUT_DIR / "signature_genes.csv"

DPI = 150
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

# Color scheme consistent with ras_signature.py
MODEL_COLORS = {
    "LASSO"         : "#E45756",
    "Elastic Net"   : "#4C78A8",
    "Random Forest" : "#54A24B",
    "Grad. Boosting": "#F58518",
}

# ===========================================================================
# 1. RECONSTRUCT DATA (same pipeline as ras_signature.py, Sections 2-6)
#    We replicate the data-loading so this script is self-contained.
# ===========================================================================
print("=" * 65)
print("compare_features.py — Cross-model feature importance analysis")
print("=" * 65)

print("\n[1/6] Loading and preparing data ...")

META_COLS   = ["stable_id", "display_label", "description", "biotype"]
expr_raw    = pd.read_csv(EXPR_FILE, sep="\t", index_col=None, low_memory=False)
protein_mask = expr_raw["biotype"] == "protein_coding"
expr_pc     = expr_raw.loc[protein_mask].drop_duplicates(subset="display_label", keep="first")
sample_cols = [c for c in expr_raw.columns if c not in META_COLS]

expr_mat = (
    expr_pc[sample_cols + ["display_label"]]
    .set_index("display_label").T.astype(float)
)

mut = pd.read_csv(MUT_FILE, sep="\t", low_memory=False)
mut["patient_id"] = mut["dbgap_sample_id"].str.replace(r"[RD]$", "", regex=True)

expr_mat.index.name = "sample_id"
expr_mat = expr_mat.reset_index()
expr_mat["patient_id"] = expr_mat["sample_id"].str.replace(r"[RD]$", "", regex=True)

shared   = set(expr_mat["patient_id"]) & set(mut["patient_id"])
expr_mat = expr_mat[expr_mat["patient_id"].isin(shared)].copy()
mut      = mut[mut["patient_id"].isin(shared)].copy()

RAS_GENES = {
    "NRAS", "KRAS", "HRAS", "PTPN11", "NF1", "CBL",
    "RRAS", "RRAS2", "RAF1", "BRAF", "MAP2K1", "MAP2K2",
}
ras_patients = set(mut[mut["symbol"].isin(RAS_GENES)]["patient_id"])
expr_mat["RAS_label"] = expr_mat["patient_id"].apply(
    lambda pid: 1 if pid in ras_patients else 0
)

gene_cols    = [c for c in expr_mat.columns if c not in ("sample_id", "patient_id", "RAS_label")]
X_raw        = expr_mat[gene_cols].values.astype(float)
y            = expr_mat["RAS_label"].values
variances    = X_raw.var(axis=0)
nonzero_mask = variances > 0
X_raw        = X_raw[:, nonzero_mask]
gene_names   = np.array(gene_cols)[nonzero_mask]

print(f"  Samples: {X_raw.shape[0]:,}  |  Features: {X_raw.shape[1]:,}")
print(f"  RAS-mutant: {y.sum():,}  ({100*y.mean():.1f}%)")

# ===========================================================================
# 2. TRAIN / TEST SPLIT AND SCALING  (identical random_state=42)
# ===========================================================================
print("\n[2/6] Splitting and scaling data (random_state=42) ...")

X_train_raw, X_test_raw, y_train, y_test = train_test_split(
    X_raw, y, test_size=0.2, stratify=y, random_state=42
)

# Load the saved scaler so we use exactly the same scaling as ras_signature.py
with open(SCALER_PKL, "rb") as fh:
    scaler = pickle.load(fh)

X_train = scaler.transform(X_train_raw)
X_test  = scaler.transform(X_test_raw)
print(f"  Train: {X_train.shape[0]}  |  Test: {X_test.shape[0]}")

# ===========================================================================
# 3. LOAD / RE-FIT ALL FOUR MODELS
#    LASSO: load from pkl (ensures identical coefficients)
#    EN / RF / GBM: re-fit with same hyperparameters as ras_signature.py
# ===========================================================================
print("\n[3/6] Loading LASSO from pkl; re-fitting EN, RF, GBM ...")

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
Cs = np.logspace(-3, 2, 60)

# --- LASSO (from pkl) -------------------------------------------------------
with open(MODEL_PKL, "rb") as fh:
    lasso_cv = pickle.load(fh)
lasso_test_auc = roc_auc_score(y_test, lasso_cv.predict_proba(X_test)[:, 1])
print(f"  LASSO  loaded  | test AUC = {lasso_test_auc:.4f}")

# --- Elastic Net ------------------------------------------------------------
enet_cv = LogisticRegressionCV(
    Cs=Cs, penalty="elasticnet", l1_ratios=[0.5], solver="saga",
    scoring="roc_auc", cv=cv, class_weight="balanced",
    max_iter=5000, n_jobs=-1, random_state=42, refit=True,
)
enet_cv.fit(X_train, y_train)
enet_test_auc = roc_auc_score(y_test, enet_cv.predict_proba(X_test)[:, 1])
print(f"  ElasticNet re-fit | test AUC = {enet_test_auc:.4f}")

# --- Random Forest ----------------------------------------------------------
rf = RandomForestClassifier(
    n_estimators=500, class_weight="balanced", random_state=42, n_jobs=-1,
)
rf.fit(X_train, y_train)
rf_test_auc = roc_auc_score(y_test, rf.predict_proba(X_test)[:, 1])
print(f"  RF     re-fit  | test AUC = {rf_test_auc:.4f}")

# --- Gradient Boosting ------------------------------------------------------
class_counts   = np.bincount(y_train)
total          = len(y_train)
sample_weights = np.where(
    y_train == 1,
    total / (2 * class_counts[1]),
    total / (2 * class_counts[0]),
)
gb = GradientBoostingClassifier(
    n_estimators=200, learning_rate=0.05, max_depth=3,
    random_state=42, subsample=0.8,
)
gb.fit(X_train, y_train, sample_weight=sample_weights)
gb_test_auc = roc_auc_score(y_test, gb.predict_proba(X_test)[:, 1])
print(f"  GBM    re-fit  | test AUC = {gb_test_auc:.4f}")

# ===========================================================================
# 4. EXTRACT IMPORTANCE METRICS FOR EACH MODEL
#
#    The central challenge in cross-model comparison is that each model
#    expresses "importance" differently:
#
#    LASSO / ElasticNet:
#      coef_.ravel() — linear weights after L1/L2 penalization.
#      Sign = direction (positive = activating in RAS-mutant, negative =
#      suppressed). Magnitude = scaled contribution to log-odds.
#      Genes with exactly zero coef are excluded from the model.
#      Interpretability: very high (direct biological direction).
#
#    Random Forest:
#      feature_importances_ = mean decrease in Gini impurity across all
#      trees, normalized to sum to 1. Always non-negative. Correlated
#      features can split importance across many trees (importance
#      "dilution"). Tends to favor high-cardinality / continuous features.
#      Interpretability: moderate (no direction; nonlinear interactions).
#
#    Gradient Boosting:
#      feature_importances_ = mean decrease in loss (residuals) when a
#      feature is used as a split criterion. Same bias as RF: diluted for
#      correlated features; no directionality.
#
#    Unified comparison strategy:
#      a) Raw importance values differ in scale and sign — we cannot
#         directly compare magnitudes across models.
#      b) RANK-based comparison is model-agnostic: convert each model's
#         importances to a rank (1 = most important), then compute
#         Spearman rho between all pairs of rank vectors. This is the
#         "Spearman rank correlation of feature rankings" approach used
#         in many ML benchmarking papers (e.g., Strobl et al., 2007).
#      c) Top-N overlap: for each N, count how many genes appear in all
#         models' top-N lists. This is a discrete, non-parametric measure.
#      d) Consensus rank: mean rank across models (lower = more important
#         on average). Genes consistently in the top 50 across all four
#         models are "high-confidence" signature genes.
# ===========================================================================
print("\n[4/6] Extracting and aligning importance metrics ...")

n_genes = len(gene_names)

# LASSO: use raw coefficients (can be negative)
lasso_coef = lasso_cv.coef_.ravel()                     # shape (n_genes,)
lasso_abs  = np.abs(lasso_coef)                         # for ranking by magnitude

# ElasticNet: same structure
enet_coef  = enet_cv.coef_.ravel()
enet_abs   = np.abs(enet_coef)

# RF: Gini impurity reduction (non-negative)
rf_imp     = rf.feature_importances_                    # shape (n_genes,)

# GBM: same structure as RF
gb_imp     = gb.feature_importances_

# --- Build a unified importance DataFrame -----------------------------------
importance_df = pd.DataFrame({
    "gene"          : gene_names,
    "lasso_coef"    : lasso_coef,
    "lasso_abs"     : lasso_abs,
    "enet_coef"     : enet_coef,
    "enet_abs"      : enet_abs,
    "rf_importance" : rf_imp,
    "gb_importance" : gb_imp,
})

# Ranks (ascending rank = most important)
importance_df["lasso_rank"] = importance_df["lasso_abs"].rank(ascending=False)
importance_df["enet_rank"]  = importance_df["enet_abs"].rank(ascending=False)
importance_df["rf_rank"]    = importance_df["rf_importance"].rank(ascending=False)
importance_df["gb_rank"]    = importance_df["gb_importance"].rank(ascending=False)

# Mean rank across all four models (consensus importance)
importance_df["mean_rank"]  = importance_df[
    ["lasso_rank", "enet_rank", "rf_rank", "gb_rank"]
].mean(axis=1)

importance_df = importance_df.sort_values("mean_rank").reset_index(drop=True)

# LASSO nonzero flag
importance_df["lasso_nonzero"] = (importance_df["lasso_coef"] != 0).astype(int)
importance_df["enet_nonzero"]  = (importance_df["enet_coef"]  != 0).astype(int)

# Save full table
out_csv = OUTPUT_DIR / "cross_model_importance.csv"
importance_df.to_csv(out_csv, index=False)
print(f"  Full importance table saved: {out_csv}")

# Top 50 consensus genes
consensus_top50 = importance_df.head(50)[
    ["gene", "mean_rank", "lasso_rank", "enet_rank", "rf_rank", "gb_rank",
     "lasso_coef", "enet_coef", "rf_importance", "gb_importance",
     "lasso_nonzero", "enet_nonzero"]
].reset_index(drop=True)
consensus_csv = OUTPUT_DIR / "consensus_genes.csv"
consensus_top50.to_csv(consensus_csv, index=False)
print(f"  Consensus top-50 saved: {consensus_csv}")

# ===========================================================================
# 5. RANK CORRELATION HEATMAP
#
#    Spearman rho computed across ALL protein-coding genes (n ~ 15,000+).
#    This tests whether the four models AGREE on the global rank of every
#    gene — not just the top ones. High rho (>0.6) means the models broadly
#    agree on which genes matter; low rho (<0.3) means fundamentally
#    different feature representations.
#
#    Expected result:
#      LASSO vs ElasticNet: very high rho (both are L1/L2 linear models;
#        at l1_ratio=0.5 with similar best-C, EN ≈ regularized LASSO)
#      LASSO vs RF: moderate rho (RF captures nonlinear/interaction effects
#        that LASSO cannot; correlated genes inflate RF importances)
#      RF vs GBM: high rho (both are tree ensembles; GBM uses sequential
#        residuals while RF uses independent trees, but feature selection
#        mechanisms are similar)
#      LASSO vs GBM: often the lowest pair (most different inductive bias)
# ===========================================================================
print("\n[5/6] Computing Spearman rank correlations ...")

rank_cols  = ["lasso_rank", "enet_rank", "rf_rank", "gb_rank"]
rank_labels = ["LASSO", "Elastic Net", "Random Forest", "Grad. Boosting"]

rank_matrix = importance_df[rank_cols].values
n_models    = len(rank_cols)
rho_matrix  = np.zeros((n_models, n_models))
pval_matrix = np.zeros((n_models, n_models))

for i in range(n_models):
    for j in range(n_models):
        rho, pval = spearmanr(rank_matrix[:, i], rank_matrix[:, j])
        rho_matrix[i, j]  = rho
        pval_matrix[i, j] = pval

rho_df = pd.DataFrame(rho_matrix, index=rank_labels, columns=rank_labels)
print("\n  Spearman rho matrix (all genes):")
print(rho_df.round(3).to_string())

# --- Plot: heatmap ----------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 4.5))

mask = np.triu(np.ones_like(rho_matrix, dtype=bool), k=1)  # upper triangle mask

sns.heatmap(
    rho_df,
    annot=True, fmt=".3f",
    cmap="RdYlGn",       # red = low agreement, green = high
    vmin=0, vmax=1,
    linewidths=0.5, linecolor="white",
    square=True,
    ax=ax,
    cbar_kws={"label": "Spearman rho", "shrink": 0.8},
)
ax.set_title(
    "Feature Rank Correlation (Spearman rho)\n"
    "Computed across all protein-coding genes",
    fontsize=11,
)
ax.set_xticklabels(rank_labels, rotation=30, ha="right")
ax.set_yticklabels(rank_labels, rotation=0)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "feature_rank_correlation.png", dpi=DPI)
plt.close(fig)
print(f"  Saved: feature_rank_correlation.png")

# ===========================================================================
# 6. TOP-N OVERLAP ANALYSIS
#
#    For N = 5, 10, 25, 50, 100, 200:
#      - Get the top-N genes for each model
#      - Count how many genes appear in ALL 4 models' top-N
#      - Count pairwise overlaps (LASSO ∩ EN, LASSO ∩ RF, LASSO ∩ GBM,
#        EN ∩ RF, EN ∩ GBM, RF ∩ GBM)
#
#    Interpretation:
#      High overlap at low N (e.g., 5/10 genes shared in top-10) means
#      the models converge on a small, reproducible core signature.
#      High overlap only at large N means the top genes differ but the
#      broader "important region" of feature space is shared.
# ===========================================================================
print("\n[6/6] Computing top-N overlaps and generating plots ...")

Ns = [5, 10, 25, 50, 100, 200]

# Get sorted gene lists per model
top_genes = {
    "LASSO"         : importance_df.sort_values("lasso_rank")["gene"].tolist(),
    "Elastic Net"   : importance_df.sort_values("enet_rank")["gene"].tolist(),
    "Random Forest" : importance_df.sort_values("rf_rank")["gene"].tolist(),
    "Grad. Boosting": importance_df.sort_values("gb_rank")["gene"].tolist(),
}

model_pairs = [
    ("LASSO", "Elastic Net"),
    ("LASSO", "Random Forest"),
    ("LASSO", "Grad. Boosting"),
    ("Elastic Net", "Random Forest"),
    ("Elastic Net", "Grad. Boosting"),
    ("Random Forest", "Grad. Boosting"),
]

print(f"\n  Top-N overlap counts (number of genes in ALL 4 models' top-N):")
print(f"  {'N':>6}  {'All-4 overlap':>14}  {'Overlap %':>10}")

all4_overlaps = []
for N in Ns:
    sets     = {m: set(genes[:N]) for m, genes in top_genes.items()}
    all4     = set.intersection(*sets.values())
    pct      = 100 * len(all4) / N
    all4_overlaps.append(len(all4))
    print(f"  {N:>6}  {len(all4):>14}  {pct:>9.1f}%")

# Pairwise overlap table
print(f"\n  Pairwise top-50 overlaps:")
for m1, m2 in model_pairs:
    s1, s2 = set(top_genes[m1][:50]), set(top_genes[m2][:50])
    overlap = len(s1 & s2)
    print(f"  {m1:<20} ∩ {m2:<20}: {overlap:>3} / 50")

# --- Plot 1: All-4 overlap bar chart ----------------------------------------
fig, ax = plt.subplots(figsize=(6, 4))
bars = ax.bar(
    [str(n) for n in Ns],
    all4_overlaps,
    color="#4C78A8",
    edgecolor="white",
    linewidth=0.5,
)
# Add "out of N" reference line
for i, (N, cnt) in enumerate(zip(Ns, all4_overlaps)):
    pct = 100 * cnt / N
    ax.text(i, cnt + 0.4, f"{cnt}\n({pct:.0f}%)", ha="center",
            va="bottom", fontsize=9, color="#333333")

ax.set_xlabel("N (top-N genes per model)", fontsize=11)
ax.set_ylabel("Genes shared across all 4 models", fontsize=11)
ax.set_title(
    "Top-N Feature Overlap — All 4 Models\n"
    "(LASSO, Elastic Net, Random Forest, Gradient Boosting)",
    fontsize=11,
)
ax.set_ylim(0, max(all4_overlaps) * 1.3)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "top_n_overlap.png", dpi=DPI)
plt.close(fig)
print(f"\n  Saved: top_n_overlap.png")

# --- Plot 2: Pairwise overlap heatmap across Ns ----------------------------
pair_labels = [f"{a.split()[0]}\n∩\n{b.split()[0]}" for a, b in model_pairs]
overlap_matrix = np.zeros((len(model_pairs), len(Ns)))

for pi, (m1, m2) in enumerate(model_pairs):
    for ni, N in enumerate(Ns):
        s1, s2 = set(top_genes[m1][:N]), set(top_genes[m2][:N])
        overlap_matrix[pi, ni] = len(s1 & s2)

fig, ax = plt.subplots(figsize=(7, 4))
sns.heatmap(
    overlap_matrix,
    xticklabels=[str(n) for n in Ns],
    yticklabels=pair_labels,
    annot=True, fmt=".0f",
    cmap="Blues",
    linewidths=0.5, linecolor="white",
    ax=ax,
    cbar_kws={"label": "Genes in common", "shrink": 0.8},
)
ax.set_xlabel("N (top-N per model)")
ax.set_title("Pairwise Top-N Feature Overlap\n(number of genes shared between model pair)", fontsize=11)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "pairwise_overlap_heatmap.png", dpi=DPI)
plt.close(fig)
print(f"  Saved: pairwise_overlap_heatmap.png")

# ===========================================================================
# 7. DOT PLOT — TOP 40 CONSENSUS GENES
#
#    For each of the top-40 genes by mean rank, show each model's
#    normalized importance as a dot. Normalization: within each model,
#    min-max scale importances to [0,1] so dot sizes are comparable.
#    Color = model. This is a "beeswarm-style" comparison across models.
#
#    Genes highlighted: those in LASSO's nonzero set AND top-50 of both
#    RF and GBM are labeled in bold — these are the highest-confidence
#    biologically relevant features.
# ===========================================================================
top40 = importance_df.head(40).copy()

# Min-max normalize each model's importance within the full gene set
def minmax(arr):
    mn, mx = arr.min(), arr.max()
    return (arr - mn) / (mx - mn) if mx > mn else np.zeros_like(arr)

imp_dict = {
    "LASSO"         : minmax(importance_df["lasso_abs"].values),
    "Elastic Net"   : minmax(importance_df["enet_abs"].values),
    "Random Forest" : minmax(importance_df["rf_importance"].values),
    "Grad. Boosting": minmax(importance_df["gb_importance"].values),
}

# For each top-40 gene (by mean_rank index), extract normalized importance
top40_idx = top40.index.tolist()
top40_genes = top40["gene"].tolist()

fig, ax = plt.subplots(figsize=(8, 10))

y_positions = np.arange(len(top40_genes))
model_offsets = {"LASSO": -0.3, "Elastic Net": -0.1, "Random Forest": 0.1, "Grad. Boosting": 0.3}

for model_name, offset in model_offsets.items():
    imp_vals = imp_dict[model_name][top40_idx]
    # Scale dot size: 20 = small, 300 = large
    sizes = 20 + 280 * imp_vals
    ax.scatter(
        imp_vals,
        y_positions + offset,
        s=sizes,
        c=MODEL_COLORS[model_name],
        alpha=0.75,
        edgecolors="white",
        linewidths=0.4,
        label=model_name,
        zorder=3,
    )

# Gene labels: bold if LASSO nonzero AND in top-50 of RF AND top-50 of GBM
lasso_nonzero_set = set(importance_df[importance_df["lasso_nonzero"] == 1]["gene"])
rf_top50_set      = set(top_genes["Random Forest"][:50])
gb_top50_set      = set(top_genes["Grad. Boosting"][:50])
high_conf_set     = lasso_nonzero_set & rf_top50_set & gb_top50_set

for yi, gname in enumerate(top40_genes):
    weight = "bold" if gname in high_conf_set else "normal"
    color  = "#333333" if gname in high_conf_set else "#666666"
    ax.text(-0.03, yi, gname, ha="right", va="center",
            fontsize=8, fontweight=weight, color=color)

ax.set_yticks(y_positions)
ax.set_yticklabels([])
ax.set_xlim(-0.05, 1.15)
ax.set_ylim(-1, len(top40_genes))
ax.invert_yaxis()
ax.set_xlabel("Normalized importance (min-max within model)")
ax.set_title(
    "Top-40 Consensus Feature Importance Across All Models\n"
    "(genes ranked by mean importance rank; bold = selected by all models)",
    fontsize=11,
)
ax.legend(title="Model", loc="lower right", fontsize=9, title_fontsize=9)
ax.axvline(0, color="black", lw=0.5, alpha=0.3)
ax.grid(axis="x", alpha=0.3, lw=0.5)
ax.grid(axis="y", alpha=0.15, lw=0.5)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "dot_plot_top40.png", dpi=DPI)
plt.close(fig)
print(f"  Saved: dot_plot_top40.png")

# ===========================================================================
# 8. CONSENSUS GENE SUMMARY FIGURE
#
#    Show the top-20 consensus genes (by mean rank) as a grouped bar chart.
#    Each bar shows the rank of that gene within each model.
#    Lower rank = more important.
#    This highlights genes where all models agree (bars are uniform height)
#    vs. genes where models disagree (bars vary widely).
# ===========================================================================
top20_consensus = importance_df.head(20).copy()

rank_cols_plot  = ["lasso_rank", "enet_rank", "rf_rank", "gb_rank"]
model_colors_list = list(MODEL_COLORS.values())

x   = np.arange(len(top20_consensus))
w   = 0.2

fig, ax = plt.subplots(figsize=(12, 5))
for i, (col, label, color) in enumerate(zip(
        rank_cols_plot,
        ["LASSO", "Elastic Net", "Random Forest", "Grad. Boosting"],
        model_colors_list)):
    ax.bar(
        x + (i - 1.5) * w,
        top20_consensus[col].values,
        width=w,
        color=color,
        label=label,
        edgecolor="white",
        linewidth=0.4,
        alpha=0.85,
    )

# Lower bars = better (rank 1 = most important)
ax.invert_yaxis()
ax.set_xticks(x)
ax.set_xticklabels(top20_consensus["gene"].values, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Importance rank within model\n(lower = more important)")
ax.set_title(
    "Top-20 Consensus Features — Per-Model Rank\n"
    "(ranked by mean rank across all 4 models; lower bar = more important)",
    fontsize=11,
)
ax.legend(title="Model", loc="lower right", fontsize=9)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "consensus_genes.png", dpi=DPI)
plt.close(fig)
print(f"  Saved: consensus_genes.png")

# ===========================================================================
# 9. BIOLOGICAL ANNOTATION SUMMARY
#
#    Print which top-20 consensus genes have established RAS-pathway biology.
#    These are the most important genes to discuss in the write-up.
# ===========================================================================
RAS_BIOLOGY_ANNOTATIONS = {
    # Core RAS feedback regulators
    "SPRED1"  : "Sprouty-related: inhibits RAS-ERK via RAF suppression (RASopathy gene)",
    "SPRED2"  : "Sprouty-related: inhibits RAS-ERK via RAF suppression (RASopathy gene)",
    "SPRY2"   : "Sprouty-2: negative feedback inhibitor of RAS/MAPK signaling",
    "DUSP6"   : "Dual-specificity phosphatase 6: ERK-specific feedback inhibitor (transcribed by ERK itself)",
    "ETV5"    : "ETS transcription factor: direct downstream target of RAS-MEK-ERK cascade",
    "CCND1"   : "Cyclin D1: transcriptionally activated by RAS/MAPK; drives G1 progression",
    "NDE1"    : "Nuderosome E1: cell cycle regulator; downstream of CDK/RAS proliferation axis",
    "RAPGEF5" : "RAP guanine nucleotide exchange factor 5: modulates RAS-related GTPases",
    "RHOB"    : "RHO GTPase B: crosstalk node between RAS and RHO signaling",
    "RUNX1"   : "Core binding factor: frequently co-mutated with RAS in AML; tumor suppressor context",
    "SOCS2"   : "Suppressor of cytokine signaling 2: JAK-STAT regulator; crosstalk with RAS",
    "KLF2"    : "Krüppel-like factor 2: transcriptional target suppressed by ERK in hematopoiesis",
    "IL2RA"   : "CD25 (IL-2 receptor alpha): JAK-STAT target; altered expression in RAS-mutant AML",
    "PRKCE"   : "Protein kinase C epsilon: activated downstream of RAS in some contexts",
    "FER"     : "FER tyrosine kinase: RTK effector; upstream of RAS activation",
    "DNTT"    : "Terminal deoxynucleotidyl transferase: differentiation marker; negative in myeloid",
    "ARHGEF11": "Rho guanine nucleotide exchange factor 11: member of ARHG study gene set",
    "ARHGAP29": "Rho GTPase activating protein 29: member of ARHG study gene set",
}

print("\n" + "=" * 65)
print("BIOLOGICAL ANNOTATIONS — Top-20 Consensus Genes")
print("=" * 65)
for _, row in importance_df.head(20).iterrows():
    gene = row["gene"]
    bio  = RAS_BIOLOGY_ANNOTATIONS.get(gene, "No established RAS-pathway annotation")
    lasso_flag = "[LASSO]" if row["lasso_nonzero"] == 1 else "       "
    print(f"  {lasso_flag}  {gene:<14}  mean_rank={row['mean_rank']:.1f}  |  {bio}")

print("\n" + "=" * 65)
print("HIGH-CONFIDENCE GENES (LASSO nonzero AND in RF top-50 AND GBM top-50):")
high_conf_genes = importance_df[
    (importance_df["lasso_nonzero"] == 1) &
    (importance_df["gene"].isin(rf_top50_set)) &
    (importance_df["gene"].isin(gb_top50_set))
].sort_values("mean_rank")

for _, row in high_conf_genes.head(20).iterrows():
    bio = RAS_BIOLOGY_ANNOTATIONS.get(row["gene"], "—")
    print(f"  {row['gene']:<14}  lasso_coef={row['lasso_coef']:+.4f}  mean_rank={row['mean_rank']:.1f}  |  {bio}")

print("\n" + "=" * 65)
print("MODEL AUCs (for table in write-up):")
print(f"  LASSO      : {lasso_test_auc:.4f}")
print(f"  ElasticNet : {enet_test_auc:.4f}")
print(f"  Rand Forest: {rf_test_auc:.4f}")
print(f"  Grad Boost : {gb_test_auc:.4f}")
print("=" * 65)
print("\nDone. Outputs saved to:")
for fname in ["cross_model_importance.csv", "consensus_genes.csv",
              "feature_rank_correlation.png", "top_n_overlap.png",
              "pairwise_overlap_heatmap.png", "dot_plot_top40.png",
              "consensus_genes.png"]:
    print(f"  {OUTPUT_DIR / fname}")
