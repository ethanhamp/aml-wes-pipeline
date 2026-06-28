"""
ras_signature.py
================
Build a RAS pathway activation signature in AML using LASSO logistic regression
on BeatAML Waves 1-4 RNA-seq data (dbGaP release).

Pipeline:
  1.  Load & match RNA / DNA data (615 patients with both modalities)
  2.  Define RAS-mutant vs. wild-type labels
  3.  Protein-coding gene filtering + variance filter
  4.  Outer train/test split (80/20, stratified) — fix for data leakage bug
  5.  StandardScaler fit on TRAIN only, applied to both train and test
  6.  Fit four models on TRAIN with 5-fold StratifiedKFold CV:
        - LASSO (L1 logistic regression)         [primary: sparse / interpretable]
        - Elastic Net (L1+L2 logistic regression) [comparison]
        - Random Forest (n=500)                   [comparison]
        - Gradient Boosting (sklearn)             [comparison]
  7.  Class imbalance sensitivity check: rerun LASSO on 1:1 downsampled TRAIN
  8.  Final AUC = score on held-out TEST set for all models
  9.  Extract signature genes from LASSO nonzero coefficients
  10. Compute per-sample RAS activation scores (LASSO predict_proba, TEST+TRAIN)
  11. Flag ARHG-mutant samples in BeatAML (exploratory; see NOTE in Section 10)
  12. Save model artifacts (scaler, model, feature gene list) for apply_signature.py
  13. Generate all diagnostic plots
  14. Print summary to stdout

Fixes vs. prior version:
  - BUG FIX: outer train/test split prevents in-sample ROC inflation (Section 4)
  - BUG FIX: scaler fit only on X_train (prevents data leakage into test set)
  - NEW: multi-model comparison with test-set AUCs (Section 6)
  - NEW: class imbalance sensitivity check via 1:1 downsampling (Section 7)
  - NEW: model artifacts saved for downstream apply_signature.py (Section 12)
"""

# ===========================================================================
# 0.  PATHS — change only these if data locations change
# ===========================================================================
import os
from pathlib import Path

DATA_DIR = Path(os.environ.get("BEATAML_DIR", str(Path(os.environ.get("ARHG_BASE_DIR", "..")).parent / "BeatAML")))
OUTPUT_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"

EXPR_FILE  = DATA_DIR / "beataml_waves1to4_norm_exp_dbgap.txt"
MUT_FILE   = DATA_DIR / "beataml_wes_wv1to4_mutations_dbgap.txt"

# Output CSVs
SCORES_CSV        = OUTPUT_DIR / "ras_activation_scores.csv"
SIG_CSV           = OUTPUT_DIR / "signature_genes.csv"
FEAT_GENES_TXT    = OUTPUT_DIR / "signature_feature_genes.txt"

# Output plots
ROC_PNG           = OUTPUT_DIR / "roc_curve.png"
CV_PNG            = OUTPUT_DIR / "cv_regularization.png"
BOX_PNG           = OUTPUT_DIR / "ras_score_by_group.png"
BAR_PNG           = OUTPUT_DIR / "top_signature_genes.png"
COMPARE_PNG       = OUTPUT_DIR / "model_comparison.png"

# Serialized model artifacts (for apply_signature.py)
SCALER_PKL        = OUTPUT_DIR / "ras_scaler.pkl"
MODEL_PKL         = OUTPUT_DIR / "ras_lasso_model.pkl"

# ===========================================================================
# 1.  IMPORTS (standard library, then third-party)
# ===========================================================================
import pickle
import warnings
warnings.filterwarnings("ignore")   # suppress sklearn convergence noise

import numpy  as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")               # non-interactive backend for saving PNGs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from sklearn.preprocessing      import StandardScaler
from sklearn.linear_model       import LogisticRegressionCV, LogisticRegression
from sklearn.ensemble           import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection    import (StratifiedKFold, train_test_split,
                                        cross_val_score)
from sklearn.metrics            import roc_auc_score, roc_curve
from sklearn.utils               import resample
from scipy.stats                import mannwhitneyu

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)
DPI = 150

# ===========================================================================
# 2.  LOAD EXPRESSION DATA
#     Rows = genes (22,843), Cols = samples (+ 4 metadata cols)
#     Transpose so rows = samples, cols = gene symbols
# ===========================================================================
print("=" * 65)
print("SECTION 2: Loading expression data ...")

expr_raw = pd.read_csv(EXPR_FILE, sep="\t", index_col=None, low_memory=False)

META_COLS   = ["stable_id", "display_label", "description", "biotype"]
sample_cols = [c for c in expr_raw.columns if c not in META_COLS]

# Protein-coding filter:
# Removes lncRNAs, pseudogenes, miRNAs, etc. — they add dimensionality without
# a clear mechanistic link to RAS signaling and inflate the feature space.
protein_mask = expr_raw["biotype"] == "protein_coding"
print(f"  Total genes         : {len(expr_raw):,}")
print(f"  Protein-coding genes: {protein_mask.sum():,}")

expr_pc = expr_raw.loc[protein_mask].copy()
expr_pc = expr_pc.drop_duplicates(subset="display_label", keep="first")

# Build samples × genes matrix (rows = RNA sample IDs, cols = gene symbols)
expr_mat = (
    expr_pc[sample_cols + ["display_label"]]
    .set_index("display_label")
    .T
    .astype(float)
)
print(f"  Expression matrix shape (samples x genes): {expr_mat.shape}")

# ===========================================================================
# 3.  LOAD MUTATION DATA & INTERSECT RNA + DNA PATIENTS
#     Strip trailing R/D suffix to create shared patient_id join key
# ===========================================================================
print("\nSECTION 3: Loading mutation data and matching patients ...")

mut = pd.read_csv(MUT_FILE, sep="\t", low_memory=False)
mut["patient_id"] = mut["dbgap_sample_id"].str.replace(r"[RD]$", "", regex=True)

expr_mat.index.name = "sample_id"
expr_mat = expr_mat.reset_index()
expr_mat["patient_id"] = expr_mat["sample_id"].str.replace(r"[RD]$", "", regex=True)

shared = set(expr_mat["patient_id"]) & set(mut["patient_id"])
print(f"  Shared (RNA + DNA) patients: {len(shared):,}")

expr_mat = expr_mat[expr_mat["patient_id"].isin(shared)].copy()
mut      = mut[mut["patient_id"].isin(shared)].copy()

# ===========================================================================
# 4.  DEFINE RAS-PATHWAY LABELS
# ===========================================================================
print("\nSECTION 4: Defining RAS-pathway labels ...")

RAS_GENES = {
    "NRAS", "KRAS", "HRAS",   # canonical RAS oncogenes (hotspot Gly12/Gly13/Gln61)
    "PTPN11",                  # SHP2: upstream RAS activator
    "NF1",                     # RAS-GAP tumor suppressor; LOF = RAS hyperactivation
    "CBL",                     # E3 ubiquitin ligase; LOF = RTK/RAS activation
    "RRAS", "RRAS2",           # R-RAS subfamily
    "RAF1", "BRAF",            # RAF kinases: directly downstream of RAS
    "MAP2K1", "MAP2K2",        # MEK1/2: downstream effectors
}

ras_patients = set(mut[mut["symbol"].isin(RAS_GENES)]["patient_id"])
expr_mat["RAS_label"] = expr_mat["patient_id"].apply(
    lambda pid: 1 if pid in ras_patients else 0
)

n_ras = expr_mat["RAS_label"].sum()
n_wt  = (expr_mat["RAS_label"] == 0).sum()
print(f"  RAS-mutant samples: {n_ras:,}")
print(f"  RAS-WT samples    : {n_wt:,}")
print(f"  RAS prevalence    : {100 * n_ras / len(expr_mat):.1f}%")

# ===========================================================================
# 5.  BUILD FULL FEATURE MATRIX (gene selection applied before split)
#     Zero-variance filter only — scaling happens AFTER split to prevent leakage
# ===========================================================================
print("\nSECTION 5: Building feature matrix ...")

gene_cols = [c for c in expr_mat.columns
             if c not in ("sample_id", "patient_id", "RAS_label")]

X_raw = expr_mat[gene_cols].values.astype(float)
y     = expr_mat["RAS_label"].values

# Zero-variance filter: constant genes carry no discriminatory signal
variances    = X_raw.var(axis=0)
nonzero_mask = variances > 0
X_raw        = X_raw[:, nonzero_mask]
gene_names   = np.array(gene_cols)[nonzero_mask]

print(f"  Genes after zero-variance filter: {X_raw.shape[1]:,}")

# Save the full protein-coding feature gene list for alignment in apply_signature.py
with open(FEAT_GENES_TXT, "w") as fh:
    fh.write("\n".join(gene_names.tolist()))
print(f"  Feature gene list saved: {FEAT_GENES_TXT}")

# ===========================================================================
# 6.  OUTER TRAIN / TEST SPLIT  [BUG FIX vs. prior version]
#
#     The prior script had NO holdout test set. LogisticRegressionCV uses
#     internal CV only for tuning C, then refits on ALL data. The ROC curve
#     shown was in-sample (optimistic / overfit).
#
#     Fix: reserve 20% of samples as a held-out test set BEFORE any fitting.
#     ALL models are trained on X_train / y_train only. Final reported AUC
#     is computed on X_test (the model has never seen these samples).
#
#     Stratify=y ensures both splits maintain the ~25% RAS-mutant proportion.
#     StandardScaler must be fit ONLY on X_train to prevent data leakage
#     (if we fit on all data, test-set mean/variance information bleeds into
#     the scaler, giving an artificially favorable normalization).
# ===========================================================================
print("\nSECTION 6: Train/test split and scaling ...")

X_train_raw, X_test_raw, y_train, y_test = train_test_split(
    X_raw, y, test_size=0.2, stratify=y, random_state=42
)
print(f"  Train: {X_train_raw.shape[0]} samples  "
      f"({y_train.sum()} RAS-mutant, {(y_train==0).sum()} WT)")
print(f"  Test : {X_test_raw.shape[0]} samples  "
      f"({y_test.sum()} RAS-mutant, {(y_test==0).sum()} WT)")

# Fit scaler on TRAIN only, transform both
scaler  = StandardScaler()
X_train = scaler.fit_transform(X_train_raw)   # fit + transform train
X_test  = scaler.transform(X_test_raw)        # transform test with TRAIN params

# ===========================================================================
# 7.  FIT FOUR MODELS (all trained on X_train, evaluated on X_test)
#
#     Model rationale:
#       LASSO    — sparse L1 penalty → interpretable signature (PRIMARY model)
#       Elastic Net — L1+L2 hybrid → handles correlated gene groups better
#       Random Forest — nonlinear ensemble; robust baseline; no scaling needed
#                       (we still apply scaled features for consistency, but RF
#                        is invariant to monotone transforms of features)
#       Gradient Boosting — sequential boosting; handles imbalance via
#                           sample_weight; different inductive bias from RF
#
#     5-fold StratifiedKFold CV is used within training data for hyperparameter
#     tuning (LASSO/ElasticNet: best C; RF/GBM: fixed hyperparams evaluated by CV).
# ===========================================================================
print("\nSECTION 7: Fitting models ...")

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
Cs = np.logspace(-3, 2, 60)   # 60 log-spaced C values from 0.001 to 100

# ---------------------------------------------------------------------------
# 7a. LASSO (L1 logistic regression)  — PRIMARY SIGNATURE MODEL
# ---------------------------------------------------------------------------
print("  [1/4] Fitting LASSO LogisticRegressionCV ...")
print("        (60 Cs x 5 folds x saga solver — may take a few minutes)")

lasso_cv = LogisticRegressionCV(
    Cs           = Cs,
    penalty      = "l1",
    solver       = "saga",        # only solver supporting L1 + warm starts
    scoring      = "roc_auc",
    cv           = cv,
    class_weight = "balanced",    # upweights RAS-mutant minority class
    max_iter     = 5000,
    n_jobs       = -1,
    random_state = 42,
    refit        = True,
)
lasso_cv.fit(X_train, y_train)

best_C     = lasso_cv.C_[0]
cv_scores  = lasso_cv.scores_[1]              # shape: (n_folds, n_Cs)
mean_auc_per_C = cv_scores.mean(axis=0)
std_auc_per_C  = cv_scores.std(axis=0)
best_C_idx     = np.argmax(mean_auc_per_C)
lasso_cv_auc   = mean_auc_per_C[best_C_idx]
lasso_cv_std   = std_auc_per_C[best_C_idx]

y_score_train_lasso = lasso_cv.predict_proba(X_train)[:, 1]
y_score_test_lasso  = lasso_cv.predict_proba(X_test)[:, 1]
lasso_test_auc      = roc_auc_score(y_test, y_score_test_lasso)

print(f"        Best C      : {best_C:.5f}")
print(f"        CV AUC      : {lasso_cv_auc:.4f} +/- {lasso_cv_std:.4f}")
print(f"        Test AUC    : {lasso_test_auc:.4f}")

# ---------------------------------------------------------------------------
# 7b. Elastic Net (L1+L2 logistic regression)
# ---------------------------------------------------------------------------
print("  [2/4] Fitting Elastic Net LogisticRegressionCV ...")

enet_cv = LogisticRegressionCV(
    Cs           = Cs,
    penalty      = "elasticnet",
    l1_ratios    = [0.5],   # mix of L1 and L2 regularization
    solver       = "saga",
    scoring      = "roc_auc",
    cv           = cv,
    class_weight = "balanced",
    max_iter     = 5000,
    n_jobs       = -1,
    random_state = 42,
    refit        = True,
)
enet_cv.fit(X_train, y_train)

# CV AUC: scores_ key is (l1_ratio, class) for elasticnet
# Best l1_ratio is stored in l1_ratio_; get its CV scores
best_l1_ratio  = enet_cv.l1_ratio_[0]
enet_cv_scores = enet_cv.scores_[1]   # shape: (n_folds, n_Cs) at best l1_ratio
enet_cv_auc    = float(enet_cv_scores.mean(axis=0).max())
enet_cv_std    = float(enet_cv_scores.std(axis=0)[enet_cv_scores.mean(axis=0).argmax()])

y_score_test_enet = enet_cv.predict_proba(X_test)[:, 1]
enet_test_auc     = roc_auc_score(y_test, y_score_test_enet)

print(f"        Best l1_ratio: {best_l1_ratio:.2f}")
print(f"        CV AUC       : {enet_cv_auc:.4f} +/- {enet_cv_std:.4f}")
print(f"        Test AUC     : {enet_test_auc:.4f}")

# ---------------------------------------------------------------------------
# 7c. Random Forest
# ---------------------------------------------------------------------------
print("  [3/4] Fitting Random Forest (n_estimators=500) ...")

rf = RandomForestClassifier(
    n_estimators = 500,
    class_weight = "balanced",
    random_state = 42,
    n_jobs       = -1,
)
rf_cv_aucs = cross_val_score(rf, X_train, y_train, cv=cv, scoring="roc_auc",
                              n_jobs=-1)
rf_cv_auc  = rf_cv_aucs.mean()
rf_cv_std  = rf_cv_aucs.std()

# Refit on full train set for test-set evaluation
rf.fit(X_train, y_train)
y_score_test_rf = rf.predict_proba(X_test)[:, 1]
rf_test_auc     = roc_auc_score(y_test, y_score_test_rf)

print(f"        CV AUC  : {rf_cv_auc:.4f} +/- {rf_cv_std:.4f}")
print(f"        Test AUC: {rf_test_auc:.4f}")

# ---------------------------------------------------------------------------
# 7d. Gradient Boosting (sklearn GradientBoostingClassifier)
#     Class imbalance handled via sample_weight (equivalent to class_weight
#     for boosting, which does not natively support class_weight= kwarg)
# ---------------------------------------------------------------------------
print("  [4/4] Fitting Gradient Boosting ...")

# Compute per-sample weights that upweight the minority class
class_counts    = np.bincount(y_train)
total           = len(y_train)
sample_weights  = np.where(
    y_train == 1,
    total / (2 * class_counts[1]),   # weight for RAS-mutant
    total / (2 * class_counts[0]),   # weight for RAS-WT
)

gb = GradientBoostingClassifier(
    n_estimators = 200,
    learning_rate = 0.05,
    max_depth    = 3,
    random_state = 42,
    subsample    = 0.8,              # stochastic boosting reduces overfitting
)

# Cross-val with sample weights requires manual loop
gb_cv_aucs = []
for train_idx, val_idx in cv.split(X_train, y_train):
    Xtr, Xval   = X_train[train_idx], X_train[val_idx]
    ytr, yval   = y_train[train_idx], y_train[val_idx]
    sw_tr       = sample_weights[train_idx]
    gb_fold     = GradientBoostingClassifier(
        n_estimators=200, learning_rate=0.05, max_depth=3,
        random_state=42, subsample=0.8,
    )
    gb_fold.fit(Xtr, ytr, sample_weight=sw_tr)
    fold_score  = roc_auc_score(yval, gb_fold.predict_proba(Xval)[:, 1])
    gb_cv_aucs.append(fold_score)

gb_cv_auc = np.mean(gb_cv_aucs)
gb_cv_std = np.std(gb_cv_aucs)

# Refit on full train set for test-set evaluation
gb.fit(X_train, y_train, sample_weight=sample_weights)
y_score_test_gb = gb.predict_proba(X_test)[:, 1]
gb_test_auc     = roc_auc_score(y_test, y_score_test_gb)

print(f"        CV AUC  : {gb_cv_auc:.4f} +/- {gb_cv_std:.4f}")
print(f"        Test AUC: {gb_test_auc:.4f}")

# ===========================================================================
# 8.  CLASS IMBALANCE SENSITIVITY CHECK
#
#     The main model uses class_weight='balanced', which is algorithmically
#     correct for L1 logistic regression under class imbalance. However,
#     to confirm results are not an artifact of the ~25/75 imbalance, we
#     run a sensitivity check: downsample the majority class (RAS-WT) to
#     match the minority class (RAS-mutant) 1:1, refit LASSO, and compare AUC.
#
#     IMPORTANT: this is ONLY a sensitivity check. The main model (trained on
#     the full training set with class_weight='balanced') is the reported model.
#     Downsampling reduces effective sample size and is not used for the
#     primary signature.
# ===========================================================================
print("\nSECTION 8: Class imbalance sensitivity check (1:1 downsampling) ...")

ras_idx = np.where(y_train == 1)[0]
wt_idx  = np.where(y_train == 0)[0]

# Downsample WT to match RAS-mutant count
wt_idx_down = resample(wt_idx, replace=False, n_samples=len(ras_idx),
                        random_state=42)
down_idx    = np.concatenate([ras_idx, wt_idx_down])
X_train_down = X_train[down_idx]
y_train_down = y_train[down_idx]

print(f"  Downsampled train: {len(y_train_down)} samples  "
      f"({y_train_down.sum()} RAS-mutant, {(y_train_down==0).sum()} WT)")

# Refit LASSO on downsampled data (class_weight='balanced' still applied,
# but with 1:1 ratio it has negligible effect — mainly serves as control)
lasso_down = LogisticRegressionCV(
    Cs           = Cs,
    penalty      = "l1",
    solver       = "saga",
    scoring      = "roc_auc",
    cv           = StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
    class_weight = "balanced",
    max_iter     = 5000,
    n_jobs       = -1,
    random_state = 42,
    refit        = True,
)
lasso_down.fit(X_train_down, y_train_down)

# Evaluate on the SAME held-out test set (not downsampled)
y_score_test_down  = lasso_down.predict_proba(X_test)[:, 1]
lasso_down_test_auc = roc_auc_score(y_test, y_score_test_down)

print(f"  LASSO (full train, balanced weights) test AUC : {lasso_test_auc:.4f}")
print(f"  LASSO (1:1 downsample sensitivity)  test AUC : {lasso_down_test_auc:.4f}")
if abs(lasso_test_auc - lasso_down_test_auc) < 0.03:
    print("  --> Results consistent: imbalance is NOT driving performance.")
else:
    print("  --> NOTE: >3% AUC difference — review class imbalance handling.")

# ===========================================================================
# 9.  EXTRACT LASSO SIGNATURE GENES (nonzero coefficients)
#
#     Why use LASSO for the signature even if RF/GBM has higher AUC?
#     LASSO's L1 penalty drives coefficients exactly to zero, producing a
#     sparse linear model. Surviving gene coefficients have a clear biological
#     interpretation: positive coefficients = genes upregulated in RAS-mutant
#     AML; negative = downregulated. This directional, gene-level interpretability
#     is essential for our biological question and for applying the signature
#     to new datasets without re-training.
# ===========================================================================
print("\nSECTION 9: Extracting LASSO signature genes ...")

coef        = lasso_cv.coef_.ravel()
nonzero_idx = np.where(coef != 0)[0]

sig_df = pd.DataFrame({
    "gene"       : gene_names[nonzero_idx],
    "coefficient": coef[nonzero_idx],
    "abs_coef"   : np.abs(coef[nonzero_idx]),
    "direction"  : np.where(coef[nonzero_idx] > 0, "activating", "suppressive"),
}).sort_values("abs_coef", ascending=False).reset_index(drop=True)

sig_df.to_csv(SIG_CSV, index=False)
print(f"  Signature genes selected: {len(sig_df):,}")
print(f"  Saved: {SIG_CSV}")

# ===========================================================================
# 10. RAS ACTIVATION SCORES (LASSO predict_proba, all samples)
#     We score ALL samples (train + test) for visualization and downstream use.
#     Scores are posterior probabilities from the LASSO model = P(RAS-mutant).
# ===========================================================================
print("\nSECTION 10: Computing RAS activation scores ...")

# Scale all samples using the train-fit scaler, then score
X_all = scaler.transform(X_raw)
y_score_all = lasso_cv.predict_proba(X_all)[:, 1]

all_indices     = np.arange(len(X_raw))
_, test_indices = train_test_split(all_indices, test_size=0.2,
                                   stratify=y, random_state=42)
split_col = ["test" if i in set(test_indices) else "train"
             for i in range(len(X_raw))]

scores_df = pd.DataFrame({
    "sample_id" : expr_mat["sample_id"].values,
    "patient_id": expr_mat["patient_id"].values,
    "RAS_score" : y_score_all,
    "RAS_label" : y,
    "split"     : split_col,
})

# ===========================================================================
# 11. FLAG ARHG-MUTANT SAMPLES IN BeatAML
#
#     NOTE: This exploratory comparison (ARHG-mutant vs. Double-WT in BeatAML)
#     is likely UNDERPOWERED because ARHGEF/ARHGAP mutations are rare somatic
#     events in public AML cohorts. A significant result here is informative,
#     but a null result does NOT rule out ARHG-RAS convergence.
#
#     The PROPER analysis is to apply the trained model to the Eisfeld Lab
#     ARHG-enriched patient cohort using apply_signature.py, which has
#     substantially more ARHG-mutant samples and targeted patient selection.
# ===========================================================================
print("Flagging ARHG-mutant samples in BeatAML (exploratory) ...")

arhg_patients = set(
    mut[mut["symbol"].str.startswith(("ARHGEF", "ARHGAP"), na=False)]["patient_id"]
)
print(f"  ARHG-mutant patients in BeatAML: {len(arhg_patients):,}")

def assign_group(row):
    if row["RAS_label"] == 1:
        return "RAS-mutant"
    elif row["patient_id"] in arhg_patients:
        return "ARHG-mutant"
    else:
        return "Double-WT"

scores_df["group"] = scores_df.apply(assign_group, axis=1)

group_counts = scores_df["group"].value_counts()
print("  Group breakdown:")
for grp, cnt in group_counts.items():
    print(f"    {grp}: {cnt}")

scores_df.to_csv(SCORES_CSV, index=False)
print(f"  Saved: {SCORES_CSV}")

# ===========================================================================
# 12. SAVE MODEL ARTIFACTS
#     These allow apply_signature.py to score new expression datasets without
#     re-training. The scaler and model are serialized with pickle (stdlib).
#     signature_feature_genes.txt (already saved in Section 5) defines the
#     gene alignment order expected by the scaler.
# ===========================================================================
print("\nSECTION 12: Saving model artifacts ...")

with open(SCALER_PKL, "wb") as fh:
    pickle.dump(scaler, fh)
print(f"  Saved scaler : {SCALER_PKL}")

with open(MODEL_PKL, "wb") as fh:
    pickle.dump(lasso_cv, fh)
print(f"  Saved model  : {MODEL_PKL}")
print(f"  Feature genes: {FEAT_GENES_TXT}")

# ===========================================================================
# 13. PLOTS
# ===========================================================================
print("\nSECTION 13: Generating plots ...")

# ---- 13a. ROC Curves (all 4 models, test set) --------------------------------
# All ROC curves are on the held-out 20% test set — no in-sample inflation.
fpr_lasso, tpr_lasso, _ = roc_curve(y_test, y_score_test_lasso)
fpr_enet,  tpr_enet,  _ = roc_curve(y_test, y_score_test_enet)
fpr_rf,    tpr_rf,    _ = roc_curve(y_test, y_score_test_rf)
fpr_gb,    tpr_gb,    _ = roc_curve(y_test, y_score_test_gb)

fig, ax = plt.subplots(figsize=(5.5, 5))
ax.plot(fpr_lasso, tpr_lasso, color="#E45756", lw=2,
        label=f"LASSO       (AUC = {lasso_test_auc:.3f})")
ax.plot(fpr_enet,  tpr_enet,  color="#4C78A8", lw=2,
        label=f"Elastic Net (AUC = {enet_test_auc:.3f})")
ax.plot(fpr_rf,    tpr_rf,    color="#54A24B", lw=2,
        label=f"Rand. Forest(AUC = {rf_test_auc:.3f})")
ax.plot(fpr_gb,    tpr_gb,    color="#F58518", lw=2,
        label=f"Grad. Boost (AUC = {gb_test_auc:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1, label="Random classifier")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC Curve — RAS Pathway Signature\n(BeatAML; held-out 20% test set)")
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(ROC_PNG, dpi=DPI)
plt.close(fig)
print(f"  Saved: {ROC_PNG}")

# ---- 13b. LASSO CV Regularization Path (trained on 80% holdout) -----------
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(np.log10(Cs), mean_auc_per_C, color="#E45756", lw=2, label="Mean CV AUC")
ax.fill_between(
    np.log10(Cs),
    mean_auc_per_C - std_auc_per_C,
    mean_auc_per_C + std_auc_per_C,
    alpha=0.25, color="#E45756", label="+/-1 SD"
)
ax.axvline(np.log10(best_C), color="black", linestyle="--",
           label=f"Best C = {best_C:.4f}")
ax.set_xlabel("log10(C)   [higher C = less regularization]")
ax.set_ylabel("CV ROC-AUC")
ax.set_title("LASSO Regularization Path — 5-Fold CV\n(trained on 80% of data; holdout test set excluded)")
ax.legend()
fig.tight_layout()
fig.savefig(CV_PNG, dpi=DPI)
plt.close(fig)
print(f"  Saved: {CV_PNG}")

# ---- 13c. RAS Score by Group (Violin + Strip) — BeatAML exploratory --------
# NOTE: ARHG-mutant group in BeatAML is typically small (underpowered).
# Use apply_signature.py on Eisfeld Lab data for the primary ARHG analysis.
GROUP_ORDER  = ["RAS-mutant", "ARHG-mutant", "Double-WT"]
GROUP_COLORS = {
    "RAS-mutant" : "#E45756",
    "ARHG-mutant": "#4C78A8",
    "Double-WT"  : "#72B7B2",
}

fig, ax = plt.subplots(figsize=(6, 5))
sns.violinplot(
    data=scores_df, x="group", y="RAS_score",
    order=GROUP_ORDER, palette=GROUP_COLORS,
    inner=None, linewidth=0.8, alpha=0.6, ax=ax,
)
sns.stripplot(
    data=scores_df, x="group", y="RAS_score",
    order=GROUP_ORDER, palette=GROUP_COLORS,
    size=2.5, alpha=0.55, jitter=True, dodge=False, ax=ax,
)

arhg_scores = scores_df.loc[scores_df["group"] == "ARHG-mutant", "RAS_score"]
dwt_scores  = scores_df.loc[scores_df["group"] == "Double-WT",   "RAS_score"]
stat, pval  = mannwhitneyu(arhg_scores, dwt_scores, alternative="two-sided")

if pval < 0.001:
    pval_str = f"p = {pval:.2e}"
elif pval < 0.05:
    pval_str = f"p = {pval:.4f}"
else:
    pval_str = f"p = {pval:.4f} (ns)"

y_max   = scores_df["RAS_score"].max()
bracket = y_max + 0.05
ax.plot([1, 1, 2, 2], [bracket, bracket + 0.03, bracket + 0.03, bracket],
        lw=1.2, color="black")
ax.text(1.5, bracket + 0.04, pval_str, ha="center", va="bottom", fontsize=9)

n_arhg = (scores_df["group"] == "ARHG-mutant").sum()
ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score (LASSO predict_proba)")
ax.set_title(
    f"RAS Activation Score by Mutation Group\n"
    f"(BeatAML exploratory; ARHG-mutant n={n_arhg} — see apply_signature.py for primary analysis)"
)
ax.set_ylim(bottom=-0.05, top=y_max + 0.18)
fig.tight_layout()
fig.savefig(BOX_PNG, dpi=DPI)
plt.close(fig)
print(f"  Saved: {BOX_PNG}")

# ---- 13d. Top 20 Signature Genes -------------------------------------------
top20  = sig_df.head(20).copy()
colors = ["#E45756" if d == "activating" else "#4C78A8"
          for d in top20["direction"]]

fig, ax = plt.subplots(figsize=(7, 6))
ax.barh(
    y=top20["gene"][::-1],
    width=top20["abs_coef"][::-1],
    color=colors[::-1],
    edgecolor="white",
    linewidth=0.5,
)
ax.set_xlabel("Absolute LASSO Coefficient")
ax.set_title("Top 20 RAS Signature Genes\n(Red = positive coef / activating, Blue = negative coef / suppressive)")
legend_elements = [
    mpatches.Patch(facecolor="#E45756", label="Positive coef (activating)"),
    mpatches.Patch(facecolor="#4C78A8", label="Negative coef (suppressive)"),
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(BAR_PNG, dpi=DPI)
plt.close(fig)
print(f"  Saved: {BAR_PNG}")

# ---- 13e. Model Comparison: AUC +/- std (CV and test) ----------------------
model_labels  = ["LASSO\n(primary)", "Elastic Net", "Random\nForest", "Grad.\nBoosting"]
cv_aucs       = [lasso_cv_auc, enet_cv_auc, rf_cv_auc, gb_cv_auc]
cv_stds       = [lasso_cv_std, enet_cv_std, rf_cv_std, gb_cv_std]
test_aucs     = [lasso_test_auc, enet_test_auc, rf_test_auc, gb_test_auc]
bar_colors    = ["#E45756", "#4C78A8", "#54A24B", "#F58518"]

x     = np.arange(len(model_labels))
width = 0.35

fig, ax = plt.subplots(figsize=(7, 5))
bars_cv   = ax.bar(x - width/2, cv_aucs, width, yerr=cv_stds, capsize=4,
                   color=bar_colors, alpha=0.65, label="5-Fold CV AUC (+/-SD)",
                   error_kw={"elinewidth": 1.5, "ecolor": "black"})
bars_test = ax.bar(x + width/2, test_aucs, width, color=bar_colors, alpha=1.0,
                   label="Held-out Test AUC",
                   edgecolor="black", linewidth=0.8)

# Add LASSO downsampling sensitivity result as a dashed line annotation
ax.axhline(lasso_down_test_auc, color="#E45756", linestyle="--", lw=1.5,
           label=f"LASSO (1:1 downsample) test AUC = {lasso_down_test_auc:.3f}")

ax.set_xticks(x)
ax.set_xticklabels(model_labels)
ax.set_ylabel("ROC-AUC")
ax.set_ylim(0.5, 1.02)
ax.set_title("Model Comparison — RAS Pathway Signature\n(BeatAML; 5-fold CV on train + test-set AUC)")
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout()
fig.savefig(COMPARE_PNG, dpi=DPI)
plt.close(fig)
print(f"  Saved: {COMPARE_PNG}")

# ===========================================================================
# 14. SUMMARY PRINTOUT
# ===========================================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"  Total samples           : {len(scores_df):,}")
print(f"  Train / Test split      : {len(X_train)} / {len(X_test)} (80/20, stratified)")
print(f"  RAS-mutant              : {n_ras:,}  ({100*n_ras/len(scores_df):.1f}%)")
print(f"  RAS-WT                  : {n_wt:,}  ({100*n_wt/len(scores_df):.1f}%)")
print(f"  ARHG-mutant (RAS-WT)    : {(scores_df['group']=='ARHG-mutant').sum():,}  "
      f"[exploratory only — underpowered]")
print(f"  Double-WT               : {(scores_df['group']=='Double-WT').sum():,}")
print()
print(f"  Feature genes (protein-coding, non-zero-variance): {len(gene_names):,}")
print(f"  Signature genes (nonzero LASSO coef)             : {len(sig_df):,}")
print()
print("  Model Performance (held-out TEST set):")
print(f"    LASSO (primary)            : CV AUC = {lasso_cv_auc:.4f} +/- {lasso_cv_std:.4f}  |  Test AUC = {lasso_test_auc:.4f}")
print(f"    Elastic Net                : CV AUC = {enet_cv_auc:.4f} +/- {enet_cv_std:.4f}  |  Test AUC = {enet_test_auc:.4f}")
print(f"    Random Forest              : CV AUC = {rf_cv_auc:.4f} +/- {rf_cv_std:.4f}  |  Test AUC = {rf_test_auc:.4f}")
print(f"    Gradient Boosting          : CV AUC = {gb_cv_auc:.4f} +/- {gb_cv_std:.4f}  |  Test AUC = {gb_test_auc:.4f}")
print()
print(f"  Class imbalance sensitivity check (1:1 downsample):")
print(f"    LASSO (full, balanced)     : Test AUC = {lasso_test_auc:.4f}")
print(f"    LASSO (1:1 downsampled)    : Test AUC = {lasso_down_test_auc:.4f}")
delta = abs(lasso_test_auc - lasso_down_test_auc)
print(f"    Delta AUC = {delta:.4f}  ({'consistent — imbalance not driving result' if delta < 0.03 else 'NOTE: review imbalance handling'})")
print()
print(f"  Best LASSO regularization C : {best_C:.5f}")
print()
print("  Top 10 signature genes:")
for _, row in sig_df.head(10).iterrows():
    sign_str = "+" if row["direction"] == "activating" else "-"
    print(f"    {sign_str}  {row['gene']:<14}  coef = {row['coefficient']:+.4f}")
print()
print(f"  Wilcoxon ARHG-mutant vs Double-WT (BeatAML exploratory):")
print(f"    U = {stat:.1f},  {pval_str}")

# Performance note
best_test_auc = max(lasso_test_auc, enet_test_auc, rf_test_auc, gb_test_auc)
best_model    = ["LASSO", "Elastic Net", "Random Forest", "Gradient Boosting"][
    [lasso_test_auc, enet_test_auc, rf_test_auc, gb_test_auc].index(best_test_auc)
]
if best_test_auc - lasso_test_auc > 0.05:
    print(f"\n  NOTE: {best_model} substantially outperforms LASSO "
          f"(delta AUC = {best_test_auc - lasso_test_auc:.3f}). "
          f"Consider whether the LASSO signature requires revision, "
          f"though LASSO is retained for biological interpretability.")

print()
print("  Model artifacts saved for apply_signature.py:")
print(f"    Scaler     : {SCALER_PKL}")
print(f"    Model      : {MODEL_PKL}")
print(f"    Feat. genes: {FEAT_GENES_TXT}")
print("=" * 65)
print("Done. Run apply_signature.py to score Eisfeld Lab ARHG-mutant patients.")
