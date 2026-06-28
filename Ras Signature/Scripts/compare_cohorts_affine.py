"""
compare_cohorts_affine.py
=========================
Re-runs the RAS pathway activation group comparison using affine-corrected
Eisfeld ARHG-mutant scores (applied_ras_scores_affine.csv).

Affine correction context:
  The affine transformation shifts and scales Eisfeld scores to match the
  BeatAML training distribution (matched on mean and variance), reducing the
  systematic low-score bias that arose from distributional mismatch between
  the Eisfeld RNA-seq pre-processing and the BeatAML training data.  This is
  the preferred analysis for publication because it places both cohorts on a
  common score axis.

Biological rationale (unchanged from compare_cohorts_ras_clean.py):
  We ask whether ARHG mutations can activate RAS transcriptional programs
  INDEPENDENTLY of direct RAS pathway mutations.  Patients who carry both an
  ARHG mutation AND a RAS pathway mutation (NRAS/KRAS/HRAS/PTPN11/NF1/CBL/
  RRAS/RRAS2/RAF1/BRAF/MAP2K1/MAP2K2) are excluded so that elevated scores
  cannot be attributed to the co-occurring RAS mutation.

Exclusion logic (two independent filters):
  1. RAS co-mutation (ras_pathway_mutants.csv, in_score_csv == True):
       Patients with confirmed somatic RAS pathway mutations — their elevated
       scores cannot be attributed to the ARHG mutation alone.
  2. No WES data (unknown RAS status):
       Patients who have an RNA-seq score but whose exome data is absent from
       arhg_mutant_unified_long.rds.  Because we lack WES for these patients
       we cannot confirm they are RAS-pathway-clean.  Including them as
       "RAS-clean" would be assumption-based, not evidence-based.  They are
       classified as unknown status and excluded.
       Patient IDs loaded at runtime from config/wes_patients.txt (gitignored).

Inputs:
  ras_activation_scores.csv         — BeatAML cohort (from ras_signature.py)
  applied_ras_scores_affine.csv     — Eisfeld ARHG-mutant, affine-corrected
  ras_pathway_mutants.csv           — exclusion list (from ras_pathway_mutants.R)

Outputs (all written to Figures/ with _affine suffix):
  Figures/cohort_comparison_violin_rasclean_affine.png
  Figures/cohort_comparison_boxplot_rasclean_affine.png
  Figures/cohort_comparison_before_after_affine.png
  Printed: exclusion details + Mann-Whitney U summary table
"""

import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# ===========================================================================
# PATHS
# ===========================================================================
BASE_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
OUTPUT_DIR    = BASE_DIR / "Figures"
BEATAML_CSV   = BASE_DIR / "ras_activation_scores.csv"
EISFELD_CSV   = BASE_DIR / "applied_ras_scores_affine.csv"   # affine-corrected
EXCLUSION_CSV = BASE_DIR / "ras_pathway_mutants.csv"

# Patients confirmed to have WES data — used to identify no-WES unknowns.
# Derived from unique patient_id values in arhg_mutant_unified_long.rds.
# Store in config/wes_patients.txt (one ID per line; gitignored).
_wes_file = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "config/wes_patients.txt"
PATIENTS_WITH_WES = set()
if _wes_file.exists():
    with open(_wes_file) as _f:
        PATIENTS_WITH_WES = {l.strip() for l in _f if l.strip() and not l.startswith("#")}

OUTPUT_DIR.mkdir(exist_ok=True)

# ===========================================================================
# LOAD DATA
# ===========================================================================
beataml    = pd.read_csv(BEATAML_CSV)
eisfeld    = pd.read_csv(EISFELD_CSV)
exclusions = pd.read_csv(EXCLUSION_CSV)

# Normalise column name: affine CSV uses 'RAS_score_affine'; rename to 'RAS_score'
# so all downstream logic is identical to compare_cohorts_ras_clean.py.
if "RAS_score_affine" in eisfeld.columns and "RAS_score" not in eisfeld.columns:
    eisfeld = eisfeld.rename(columns={"RAS_score_affine": "RAS_score"})

print("=" * 65)
print("compare_cohorts_affine.py")
print("=" * 65)
print(f"\nLoaded {len(eisfeld)} Eisfeld ARHG-mutant scored samples (affine-corrected)")
print(f"Loaded {len(beataml)} BeatAML scored samples")
print(f"Loaded {len(exclusions)} rows from exclusion table "
      f"({exclusions['patient_id'].nunique()} patients)")

# ===========================================================================
# BUILD EXCLUSION SET — Stage 1: confirmed RAS pathway co-mutations
# Patients to exclude from the Eisfeld group:
#   in_score_csv == True → they have both a RAS mutation AND an RNA-seq score.
#   Patients with in_score_csv == False had no RNA-seq data and are already
#   absent from applied_ras_scores_affine.csv — no action needed for those.
# ===========================================================================
to_exclude_ras = (
    exclusions[exclusions["in_score_csv"] == True]          # noqa: E712
    [["patient_id", "ras_score_id", "gene", "hgvs_mane", "vaf", "timepoint"]]
    .drop_duplicates(subset="ras_score_id")                  # one row per patient
    .sort_values("ras_score_id")
)

exclude_ids_ras = set(to_exclude_ras["ras_score_id"].dropna())

print(f"\n{'─'*65}")
print("EXCLUSIONS (Stage 1) — Eisfeld patients with somatic RAS pathway mutations")
print(f"{'─'*65}")
print(f"  Patients to exclude: {len(exclude_ids_ras)}")
print()

for _, row in to_exclude_ras.iterrows():
    patient_variants = exclusions[exclusions["ras_score_id"] == row["ras_score_id"]]
    var_summary = "; ".join(
        f"{r['gene']} {r['hgvs_mane']} (VAF={r['vaf']:.3f}, {r['timepoint']})"
        for _, r in patient_variants.iterrows()
    )
    print(f"  EXCLUDE  {row['ras_score_id']}")
    print(f"           RAS mutations: {var_summary}")
    score = eisfeld.loc[eisfeld["sample_id"] == row["ras_score_id"], "RAS_score"]
    if len(score) > 0:
        print(f"           RAS score (affine): {score.values[0]:.4f}")
    print()

# ===========================================================================
# BUILD EXCLUSION SET — Stage 2: no WES data → unknown RAS pathway status
# A patient can only be confirmed RAS-clean if we have exome data for them.
# Patients with an RNA-seq score but no WES entry in arhg_mutant_unified_long
# are NOT confirmed clean — they are unknown.  Treating them as RAS-clean
# would silently assume absence of evidence equals evidence of absence.
# These patients are excluded under the "unknown status" category.
# ===========================================================================

def strip_blast(sample_id: str) -> str:
    """Convert e.g. '4472_blast' → '4472'; pass-through for plain IDs."""
    return sample_id.replace("_blast", "")

# Among scored patients not already excluded by Stage 1, identify those
# whose patient_id (after stripping _blast suffix) is absent from WES set.
remaining_after_stage1 = eisfeld[~eisfeld["sample_id"].isin(exclude_ids_ras)].copy()
no_wes_mask = remaining_after_stage1["sample_id"].apply(
    lambda sid: strip_blast(sid) not in PATIENTS_WITH_WES
)
to_exclude_no_wes = remaining_after_stage1.loc[no_wes_mask, "sample_id"].tolist()
exclude_ids_no_wes = set(to_exclude_no_wes)

print(f"{'─'*65}")
print("EXCLUSIONS (Stage 2) — Scored patients with NO WES data (unknown RAS status)")
print(f"{'─'*65}")
print(f"  Patients to exclude: {len(exclude_ids_no_wes)}")
print()
for sid in sorted(exclude_ids_no_wes):
    score = eisfeld.loc[eisfeld["sample_id"] == sid, "RAS_score"]
    score_str = f"{score.values[0]:.4f}" if len(score) > 0 else "N/A"
    print(f"  EXCLUDE  {sid}")
    print(f"           Reason: absent from arhg_mutant_unified_long.rds "
          f"(no WES → RAS status unknown)")
    print(f"           RAS score (affine): {score_str}")
    print()

# Combined exclusion set
exclude_ids_all = exclude_ids_ras | exclude_ids_no_wes

# ===========================================================================
# FILTER EISFELD TO RAS-CLEAN SUBSET
# ===========================================================================
eisfeld_clean = eisfeld[~eisfeld["sample_id"].isin(exclude_ids_all)].copy()
n_original        = len(eisfeld)
n_excluded_ras    = len(exclude_ids_ras)
n_excluded_no_wes = len(exclude_ids_no_wes)
n_excluded        = len(exclude_ids_all)
n_clean           = len(eisfeld_clean)

print(f"{'─'*65}")
print(f"  Original Eisfeld ARHG-mutant     : {n_original}")
print(f"  Excluded (RAS co-mutated)         : {n_excluded_ras}")
print(f"  Excluded (no WES / unknown status): {n_excluded_no_wes}")
print(f"  Total excluded                    : {n_excluded}")
print(f"  RAS-clean subset (confirmed)      : {n_clean}")
print(f"{'─'*65}\n")

no_score = exclusions[exclusions["in_score_csv"] == False]["patient_id"].unique()  # noqa: E712
if len(no_score) > 0:
    print(f"Note: {len(no_score)} additional RAS-mutant patients "
          f"({', '.join(sorted(no_score))})")
    print(f"      had no RNA-seq data and were not in applied_ras_scores_affine.csv.")
    print(f"      These are already absent from the analysis (no action needed).\n")

# ===========================================================================
# COMBINE FOR PLOTTING
# ===========================================================================
# Label uses explicit '(affine)' annotation so plots are self-documenting.
EISFELD_LABEL = "Eisfeld ARHG-mutant\n(RAS-clean, WES-confirmed, affine)"
EISFELD_ALL_LABEL = f"Eisfeld ARHG-mutant\n(all, n={n_original}, affine)"

beataml_plot = beataml[["sample_id", "RAS_score", "group"]].copy()
eisfeld_plot = eisfeld_clean[["sample_id", "RAS_score"]].copy()
eisfeld_plot["group"] = EISFELD_LABEL

combined = pd.concat([beataml_plot, eisfeld_plot], ignore_index=True)

# Pre-exclusion version for before/after panel
eisfeld_orig = eisfeld[["sample_id", "RAS_score"]].copy()
eisfeld_orig["group"] = EISFELD_ALL_LABEL

# ===========================================================================
# GROUP ORDER & COLORS
# ===========================================================================
GROUP_ORDER = [
    "RAS-mutant",
    EISFELD_LABEL,
    "ARHG-mutant",   # BeatAML ARHG — underpowered reference
    "Double-WT",
]
GROUP_ORDER = [g for g in GROUP_ORDER if g in combined["group"].unique()]

PALETTE = {
    "RAS-mutant"    : "#E45756",   # red
    EISFELD_LABEL   : "#F28E2B",   # orange — primary result
    "ARHG-mutant"   : "#76B7B2",   # muted teal — underpowered BeatAML reference
    "Double-WT"     : "#BAB0AC",   # grey
    EISFELD_ALL_LABEL: "#AEC7E8",  # muted blue — pre-exclusion
}

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)

# ===========================================================================
# HELPER: draw significance bracket
# ===========================================================================
def draw_bracket(ax, x1, x2, y, p, tip=0.02):
    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    ax.plot([x1, x1, x2, x2],
            [y, y + tip, y + tip, y],
            lw=1.2, color="black")
    ax.text((x1 + x2) / 2, y + tip + 0.005, stars,
            ha="center", va="bottom", fontsize=11)


# ===========================================================================
# PLOT 1 — Violin + strip (all groups)
# ===========================================================================
fig, ax = plt.subplots(figsize=(7, 5))

sns.violinplot(
    data=combined, x="group", y="RAS_score", hue="group",
    order=GROUP_ORDER, palette=PALETTE,
    inner=None, linewidth=0.8, alpha=0.55, legend=False, ax=ax,
)
sns.stripplot(
    data=combined, x="group", y="RAS_score", hue="group",
    order=GROUP_ORDER, palette=PALETTE,
    size=3, alpha=0.6, jitter=True, legend=False, ax=ax,
)

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5,
           label="Score = 0.5 threshold")
ax.legend(loc="upper right", fontsize=8)
ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score (affine-corrected)")
ax.set_title(
    "RAS Pathway Activation — Eisfeld ARHG-mutant (RAS-clean, WES-confirmed, affine) vs BeatAML\n"
    f"(n={n_excluded_ras} RAS co-mutated + n={n_excluded_no_wes} no-WES excluded; "
    f"BeatAML-trained LASSO; affine score correction)"
)

for i, grp in enumerate(GROUP_ORDER):
    grp_scores = combined.loc[combined["group"] == grp, "RAS_score"]
    n   = len(grp_scores)
    pct = 100 * (grp_scores > 0.5).mean()
    ax.text(i, -0.07, f"n={n}\n{pct:.0f}% >0.5", ha="center", va="top",
            fontsize=8, transform=ax.get_xaxis_transform())

if "ARHG-mutant" in GROUP_ORDER:
    idx = GROUP_ORDER.index("ARHG-mutant")
    ax.text(idx, 1.02, "ref\n(BeatAML)", ha="center", va="bottom",
            fontsize=7, color="grey", transform=ax.get_xaxis_transform())

fig.tight_layout()
out1 = OUTPUT_DIR / "cohort_comparison_violin_rasclean_affine.png"
fig.savefig(out1, dpi=150)
plt.close(fig)
print(f"Saved: {out1}")


# ===========================================================================
# PLOT 2 — Boxplot (3-group: RAS-mut, Eisfeld RAS-clean, Double-WT)
#           with significance brackets
# ===========================================================================
plot_groups = ["RAS-mutant", EISFELD_LABEL, "Double-WT"]
plot_groups = [g for g in plot_groups if g in combined["group"].unique()]
plot_data   = combined[combined["group"].isin(plot_groups)].copy()

eisfeld_clean_scores = combined.loc[combined["group"] == EISFELD_LABEL, "RAS_score"]

fig, ax = plt.subplots(figsize=(6, 5))

sns.boxplot(
    data=plot_data, x="group", y="RAS_score", hue="group",
    order=plot_groups, palette=PALETTE,
    width=0.5, linewidth=1.2, fliersize=0, legend=False, ax=ax,
)
sns.stripplot(
    data=plot_data, x="group", y="RAS_score", hue="group",
    order=plot_groups, palette=PALETTE,
    size=3.5, alpha=0.55, jitter=True, legend=False, ax=ax,
)

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5,
           label="Score = 0.5 threshold")
ax.legend(loc="upper right", fontsize=8)

y_max = plot_data["RAS_score"].max()
eisfeld_idx = plot_groups.index(EISFELD_LABEL)

# Eisfeld RAS-clean (affine) vs RAS-mutant
if "RAS-mutant" in plot_groups:
    ras_idx = plot_groups.index("RAS-mutant")
    other   = combined.loc[combined["group"] == "RAS-mutant", "RAS_score"]
    _, p    = mannwhitneyu(eisfeld_clean_scores, other, alternative="two-sided")
    draw_bracket(ax, ras_idx, eisfeld_idx, y_max + 0.03, p)

# Eisfeld RAS-clean (affine) vs Double-WT
if "Double-WT" in plot_groups:
    wt_idx = plot_groups.index("Double-WT")
    other  = combined.loc[combined["group"] == "Double-WT", "RAS_score"]
    _, p   = mannwhitneyu(eisfeld_clean_scores, other, alternative="two-sided")
    draw_bracket(ax, eisfeld_idx, wt_idx, y_max + 0.12, p)

for i, grp in enumerate(plot_groups):
    n = (plot_data["group"] == grp).sum()
    ax.text(i, -0.07, f"n={n}", ha="center", va="top",
            fontsize=9, transform=ax.get_xaxis_transform())

n_high  = (eisfeld_clean_scores > 0.5).sum()
n_total = len(eisfeld_clean_scores)
ax.text(
    eisfeld_idx, y_max + 0.20,
    f"{n_high}/{n_total}\n>0.5",
    ha="center", va="bottom", fontsize=8,
    color=PALETTE[EISFELD_LABEL],
    fontweight="bold",
)

ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score (affine-corrected)")
ax.set_title(
    f"RAS Activation: ARHG-mutant (RAS-clean, WES-confirmed, n={n_clean}) vs BeatAML Groups\n"
    f"(n={n_excluded_ras} RAS co-mut + n={n_excluded_no_wes} no-WES excluded; "
    f"affine correction; dashed = 0.5)"
)
ax.set_ylim(bottom=-0.05, top=y_max + 0.35)

fig.tight_layout()
out2 = OUTPUT_DIR / "cohort_comparison_boxplot_rasclean_affine.png"
fig.savefig(out2, dpi=150)
plt.close(fig)
print(f"Saved: {out2}")


# ===========================================================================
# PLOT 3 — Before/after panel (original vs RAS-clean, affine-corrected)
# ===========================================================================
combined_ba = pd.concat([
    beataml[["sample_id", "RAS_score", "group"]],
    eisfeld_orig[["sample_id", "RAS_score", "group"]],
    eisfeld_plot[["sample_id", "RAS_score", "group"]],
], ignore_index=True)

ba_order = [
    "RAS-mutant",
    "Double-WT",
    EISFELD_ALL_LABEL,
    EISFELD_LABEL,
]
ba_order = [g for g in ba_order if g in combined_ba["group"].unique()]

ba_palette = {**PALETTE}  # inherits all keys defined above

fig, ax = plt.subplots(figsize=(8, 5))

sns.boxplot(
    data=combined_ba, x="group", y="RAS_score", hue="group",
    order=ba_order, palette=ba_palette,
    width=0.5, linewidth=1.2, fliersize=0, legend=False, ax=ax,
)
sns.stripplot(
    data=combined_ba, x="group", y="RAS_score", hue="group",
    order=ba_order, palette=ba_palette,
    size=3.5, alpha=0.6, jitter=True, legend=False, ax=ax,
)

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.4)

for i, grp in enumerate(ba_order):
    scores = combined_ba.loc[combined_ba["group"] == grp, "RAS_score"]
    n   = len(scores)
    med = scores.median()
    ax.text(i, -0.07, f"n={n}\nmed={med:.2f}", ha="center", va="top",
            fontsize=7.5, transform=ax.get_xaxis_transform())

wt_scores = combined_ba.loc[combined_ba["group"] == "Double-WT", "RAS_score"]
y_top = combined_ba["RAS_score"].max()

# All (pre-exclusion) vs Double-WT
if EISFELD_ALL_LABEL in ba_order and "Double-WT" in ba_order:
    all_scores = combined_ba.loc[combined_ba["group"] == EISFELD_ALL_LABEL, "RAS_score"]
    _, p = mannwhitneyu(all_scores, wt_scores, alternative="two-sided")
    all_idx = ba_order.index(EISFELD_ALL_LABEL)
    wt_idx  = ba_order.index("Double-WT")
    draw_bracket(ax, wt_idx, all_idx, y_top + 0.03, p)

# RAS-clean (affine) vs Double-WT — separate bracket to show effect of exclusion
if EISFELD_LABEL in ba_order and "Double-WT" in ba_order:
    clean_scores = combined_ba.loc[combined_ba["group"] == EISFELD_LABEL, "RAS_score"]
    _, p = mannwhitneyu(clean_scores, wt_scores, alternative="two-sided")
    all_idx   = ba_order.index(EISFELD_ALL_LABEL)
    clean_idx = ba_order.index(EISFELD_LABEL)
    draw_bracket(ax, all_idx, clean_idx, y_top + 0.03, p)

ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score (affine-corrected)")
ax.set_title(
    "Before/After Exclusions (affine-corrected scores)\n"
    f"(n={n_excluded_ras} RAS co-mutated + n={n_excluded_no_wes} no-WES/unknown removed)"
)
ax.set_ylim(bottom=-0.12, top=y_top + 0.35)

fig.tight_layout()
out3 = OUTPUT_DIR / "cohort_comparison_before_after_affine.png"
fig.savefig(out3, dpi=150)
plt.close(fig)
print(f"Saved: {out3}")


# ===========================================================================
# WILCOXON SUMMARY TABLE
# ===========================================================================
print("\n" + "=" * 65)
print("WILCOXON TESTS — Eisfeld ARHG-mutant (RAS-clean, WES-confirmed, affine-corrected)")
print(f"  n = {n_clean}")
print("=" * 65)

comparisons = [
    ("RAS-mutant", "vs RAS-mutant"),
    ("Double-WT",  "vs Double-WT"),
]
if "ARHG-mutant" in combined["group"].unique():
    comparisons.append(("ARHG-mutant", "vs BeatAML ARHG-mutant"))

for grp, label in comparisons:
    other_scores = combined.loc[combined["group"] == grp, "RAS_score"]
    if len(other_scores) == 0:
        continue
    stat, p = mannwhitneyu(eisfeld_clean_scores, other_scores, alternative="two-sided")
    stars   = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {label:<30} U={stat:.0f}, p={p:.4f} {stars}")


# ===========================================================================
# FULL GROUP SUMMARY
# ===========================================================================
print("\n" + "=" * 65)
print("GROUP SUMMARY (affine-corrected Eisfeld scores)")
print("=" * 65)
all_groups = list(combined["group"].unique()) + [EISFELD_ALL_LABEL]
for grp in all_groups:
    src    = combined if EISFELD_ALL_LABEL not in grp else combined_ba
    scores = src.loc[src["group"] == grp, "RAS_score"]
    if len(scores) == 0:
        continue
    marker = " <- RAS-clean (primary result)" if "RAS-clean" in grp else \
             " <- pre-exclusion comparison"   if "all"       in grp else ""
    print(f"  {grp.replace(chr(10), ' '):<50}  n={len(scores):3d}  "
          f"median={scores.median():.3f}  mean={scores.mean():.3f}{marker}")

print("=" * 65)
print(f"\nExcluded IDs (RAS co-mutated): {sorted(exclude_ids_ras)}")
print(f"Excluded IDs (no WES / unknown): {sorted(exclude_ids_no_wes)}")
print("\nDone. Plots saved to:", OUTPUT_DIR)
