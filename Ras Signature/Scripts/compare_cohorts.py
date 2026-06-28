"""
compare_cohorts.py
==================
Compare RAS activation scores across:
  - BeatAML RAS-mutant
  - BeatAML Double-WT
  - BeatAML ARHG-mutant (underpowered — kept for reference)
  - Eisfeld Lab ARHG-mutant (primary comparison group)

Inputs:
  ras_activation_scores.csv  — BeatAML cohort (from ras_signature.py)
  applied_ras_scores.csv     — Eisfeld Lab ARHG-mutant (from apply_signature.py --normalize)

NOTE: applied_ras_scores.csv must be generated with the --normalize flag:
  python apply_signature.py eisfeld_expression.tsv --normalize
"""

import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import itertools

# ===========================================================================
# PATHS
# ===========================================================================
OUTPUT_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
BEATAML_CSV  = OUTPUT_DIR / "ras_activation_scores.csv"
EISFELD_CSV  = OUTPUT_DIR / "applied_ras_scores.csv"

# ===========================================================================
# LOAD & COMBINE
# ===========================================================================
beataml = pd.read_csv(BEATAML_CSV)
eisfeld = pd.read_csv(EISFELD_CSV)

# Label Eisfeld samples as their own group
eisfeld["group"] = "Eisfeld ARHG-mutant"

# Combine into one dataframe with consistent columns
beataml_plot = beataml[["sample_id", "RAS_score", "group"]].copy()
eisfeld_plot  = eisfeld[["sample_id", "RAS_score", "group"]].copy()

combined = pd.concat([beataml_plot, eisfeld_plot], ignore_index=True)

# ===========================================================================
# GROUP ORDER & COLORS
# ===========================================================================
GROUP_ORDER = [
    "RAS-mutant",
    "Eisfeld ARHG-mutant",
    "ARHG-mutant",       # BeatAML ARHG (underpowered, reference only)
    "Double-WT",
]

# Filter to only groups present in data
GROUP_ORDER = [g for g in GROUP_ORDER if g in combined["group"].unique()]

PALETTE = {
    "RAS-mutant"          : "#E45756",   # red
    "Eisfeld ARHG-mutant" : "#F28E2B",   # orange — primary result
    "ARHG-mutant"         : "#76B7B2",   # muted teal — underpowered reference
    "Double-WT"           : "#BAB0AC",   # grey
}

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)

# ===========================================================================
# PLOT 1: Violin + strip across all groups
# ===========================================================================
fig, ax = plt.subplots(figsize=(7, 5))

sns.violinplot(
    data=combined, x="group", y="RAS_score",
    order=GROUP_ORDER, palette=PALETTE,
    inner=None, linewidth=0.8, alpha=0.55, ax=ax,
)
sns.stripplot(
    data=combined, x="group", y="RAS_score",
    order=GROUP_ORDER, palette=PALETTE,
    size=3, alpha=0.6, jitter=True, ax=ax,
)

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5,
           label="Score = 0.5 threshold")
ax.legend(loc="upper right", fontsize=8)

ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score")
ax.set_title("RAS Pathway Activation Score by Mutation Group\n"
             "(BeatAML-trained LASSO; Eisfeld scores per-gene z-score normalized)")

# Add n= and % > 0.5 labels below each group
for i, grp in enumerate(GROUP_ORDER):
    grp_scores = combined.loc[combined["group"] == grp, "RAS_score"]
    n = len(grp_scores)
    pct_high = 100 * (grp_scores > 0.5).mean()
    ax.text(i, -0.07, f"n={n}\n{pct_high:.0f}% >0.5", ha="center", va="top",
            fontsize=8, transform=ax.get_xaxis_transform())

# Annotate BeatAML ARHG as underpowered if present
if "ARHG-mutant" in GROUP_ORDER:
    idx = GROUP_ORDER.index("ARHG-mutant")
    ax.text(idx, 1.02, "ref\n(BeatAML)", ha="center", va="bottom",
            fontsize=7, color="grey", transform=ax.get_xaxis_transform())

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "cohort_comparison_violin.png", dpi=150)
plt.close(fig)
print("Saved: cohort_comparison_violin.png")

# ===========================================================================
# PLOT 2: Pairwise Wilcoxon p-values — Eisfeld ARHG vs key groups
# ===========================================================================
eisfeld_scores = combined.loc[combined["group"] == "Eisfeld ARHG-mutant", "RAS_score"]

comparisons = [
    ("RAS-mutant",  "vs RAS-mutant"),
    ("Double-WT",   "vs Double-WT"),
]
if "ARHG-mutant" in GROUP_ORDER:
    comparisons.append(("ARHG-mutant", "vs BeatAML ARHG-mutant"))

print("\nWilcoxon tests — Eisfeld ARHG-mutant:")
print(f"  n = {len(eisfeld_scores)}")
for grp, label in comparisons:
    other = combined.loc[combined["group"] == grp, "RAS_score"]
    stat, p = mannwhitneyu(eisfeld_scores, other, alternative="two-sided")
    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {label}: U={stat:.0f}, p={p:.4f} {stars}")

# ===========================================================================
# PLOT 3: Boxplot with significance brackets (Eisfeld ARHG vs RAS-mut & WT)
# ===========================================================================
plot_groups = ["RAS-mutant", "Eisfeld ARHG-mutant", "Double-WT"]
plot_groups = [g for g in plot_groups if g in combined["group"].unique()]
plot_data   = combined[combined["group"].isin(plot_groups)].copy()

fig, ax = plt.subplots(figsize=(6, 5))

sns.boxplot(
    data=plot_data, x="group", y="RAS_score",
    order=plot_groups, palette=PALETTE,
    width=0.5, linewidth=1.2, fliersize=0, ax=ax,
)
sns.stripplot(
    data=plot_data, x="group", y="RAS_score",
    order=plot_groups, palette=PALETTE,
    size=3.5, alpha=0.55, jitter=True, ax=ax,
)

# Significance brackets
y_max = plot_data["RAS_score"].max()
bracket_height = 0.06
tip_height     = 0.02

def draw_bracket(ax, x1, x2, y, p):
    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    ax.plot([x1, x1, x2, x2],
            [y, y + tip_height, y + tip_height, y],
            lw=1.2, color="black")
    ax.text((x1 + x2) / 2, y + tip_height + 0.005, stars,
            ha="center", va="bottom", fontsize=11)

eisfeld_idx = plot_groups.index("Eisfeld ARHG-mutant")

# Eisfeld ARHG vs RAS-mutant
if "RAS-mutant" in plot_groups:
    ras_idx = plot_groups.index("RAS-mutant")
    other   = combined.loc[combined["group"] == "RAS-mutant", "RAS_score"]
    _, p    = mannwhitneyu(eisfeld_scores, other, alternative="two-sided")
    draw_bracket(ax, ras_idx, eisfeld_idx, y_max + 0.03, p)

# Eisfeld ARHG vs Double-WT
if "Double-WT" in plot_groups:
    wt_idx = plot_groups.index("Double-WT")
    other  = combined.loc[combined["group"] == "Double-WT", "RAS_score"]
    _, p   = mannwhitneyu(eisfeld_scores, other, alternative="two-sided")
    draw_bracket(ax, eisfeld_idx, wt_idx, y_max + 0.12, p)

# n= labels
for i, grp in enumerate(plot_groups):
    n = (plot_data["group"] == grp).sum()
    ax.text(i, -0.07, f"n={n}", ha="center", va="top",
            fontsize=9, transform=ax.get_xaxis_transform())

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5,
           label="Score = 0.5 threshold")
ax.legend(loc="upper right", fontsize=8)

# Annotate Eisfeld high-scorers count
n_high = (eisfeld_scores > 0.5).sum()
n_total = len(eisfeld_scores)
if "Eisfeld ARHG-mutant" in plot_groups:
    eis_x = plot_groups.index("Eisfeld ARHG-mutant")
    ax.text(eis_x, y_max + 0.20,
            f"{n_high}/{n_total}\n>0.5",
            ha="center", va="bottom", fontsize=8, color=PALETTE["Eisfeld ARHG-mutant"],
            fontweight="bold")

ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score")
ax.set_title("RAS Activation: Eisfeld ARHG-mutant vs BeatAML Groups\n"
             "(per-gene z-score normalized; dashed line = 0.5 threshold)")
ax.set_ylim(bottom=-0.05, top=y_max + 0.35)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "cohort_comparison_boxplot.png", dpi=150)
plt.close(fig)
print("Saved: cohort_comparison_boxplot.png")

# ===========================================================================
# SUMMARY
# ===========================================================================
print("\n" + "=" * 50)
print("GROUP SUMMARY")
print("=" * 50)
for grp in GROUP_ORDER:
    scores = combined.loc[combined["group"] == grp, "RAS_score"]
    print(f"  {grp:<25} n={len(scores):3d}  "
          f"median={scores.median():.3f}  "
          f"mean={scores.mean():.3f}")
print("=" * 50)
print("Done. Plots saved to:", OUTPUT_DIR)
