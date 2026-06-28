"""
compare_cohorts_logcpm.py
=========================
Regenerates cohort comparison figures using log2(CPM+1) normalized Eisfeld
scores instead of the z-score approach. Old figures are preserved; new ones
have _logcpm suffix and include "log-CPM normalized" in all titles.

Figures produced (Figures/ directory):
  cohort_comparison_boxplot_logcpm.png
  ras_comutat_validation_logcpm.png
  ras_comutat_validation_logcpm.pdf
  ras_score_by_gene_logcpm.pdf
  ras_score_gap_vs_gef_logcpm.pdf
  ras_score_4groups_logcpm.pdf

Inputs:
  ras_activation_scores.csv          — BeatAML cohort scores
  applied_ras_scores_logcpm.csv      — Eisfeld log-CPM scores (apply_signature_logcpm.py)
  ras_pathway_mutants.csv            — RAS co-mutation exclusion/highlight list
  ../Variants/arhg_variants_clean.txt — per-patient ARHG gene assignments
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"
FIG_DIR    = BASE / "Figures"
BEATAML    = BASE / "ras_activation_scores.csv"
EISFELD    = BASE / "applied_ras_scores_logcpm.csv"
EXCLUSIONS = BASE / "ras_pathway_mutants.csv"
VARIANTS = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Variants/arhg_variants_clean.txt"

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)

PALETTE = {
    "RAS-mutant"          : "#E45756",
    "Eisfeld ARHG-mutant" : "#F28E2B",
    "Double-WT"           : "#BAB0AC",
}
NORM_LABEL = "log-CPM normalized"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
beataml = pd.read_csv(BEATAML)
eisfeld = pd.read_csv(EISFELD).rename(columns={"RAS_score_logcpm": "RAS_score"})
eisfeld["group"] = "Eisfeld ARHG-mutant"

exclusions = pd.read_csv(EXCLUSIONS)
# Patients present in the scoring CSV who also carry RAS mutations
ras_comut_ids = set(
    exclusions.loc[exclusions["in_score_csv"] == True, "ras_score_id"]
    .str.replace("_blast", "", regex=False)
)

# RAS-clean Eisfeld subset: exclude patients with RAS pathway co-mutations
eisfeld_clean = eisfeld[~eisfeld["sample_id"].isin(ras_comut_ids)].copy()
print(f"Eisfeld RAS-clean: {len(eisfeld_clean)}/{len(eisfeld)} patients (excluded {len(eisfeld)-len(eisfeld_clean)} RAS co-mutant)")

combined = pd.concat(
    [beataml[["sample_id", "RAS_score", "group"]],
     eisfeld_clean[["sample_id", "RAS_score", "group"]]],
    ignore_index=True
)

eisfeld_scores = eisfeld_clean["RAS_score"]
GROUP_ORDER_BOX = ["RAS-mutant", "Eisfeld ARHG-mutant", "Double-WT"]

# ---------------------------------------------------------------------------
# Helper: significance bracket
# ---------------------------------------------------------------------------
def draw_bracket(ax, x1, x2, y, p, tip=0.02):
    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    ax.plot([x1, x1, x2, x2], [y, y + tip, y + tip, y], lw=1.2, color="black")
    ax.text((x1 + x2) / 2, y + tip + 0.005, stars, ha="center", va="bottom", fontsize=11)

# ===========================================================================
# FIGURE 1 — Cohort comparison boxplot
# ===========================================================================
plot_data = combined[combined["group"].isin(GROUP_ORDER_BOX)].copy()
y_max = plot_data["RAS_score"].max()

fig, ax = plt.subplots(figsize=(6, 5))

sns.boxplot(
    data=plot_data, x="group", y="RAS_score",
    order=GROUP_ORDER_BOX, palette=PALETTE,
    width=0.5, linewidth=1.2, fliersize=0, ax=ax,
)
sns.stripplot(
    data=plot_data, x="group", y="RAS_score",
    order=GROUP_ORDER_BOX, palette=PALETTE,
    size=3.5, alpha=0.55, jitter=True, ax=ax,
)

eis_idx = GROUP_ORDER_BOX.index("Eisfeld ARHG-mutant")
ras_idx = GROUP_ORDER_BOX.index("RAS-mutant")
wt_idx  = GROUP_ORDER_BOX.index("Double-WT")

_, p_ras = mannwhitneyu(eisfeld_scores,
                         combined.loc[combined["group"] == "RAS-mutant", "RAS_score"],
                         alternative="two-sided")
_, p_wt  = mannwhitneyu(eisfeld_scores,
                         combined.loc[combined["group"] == "Double-WT", "RAS_score"],
                         alternative="two-sided")

draw_bracket(ax, ras_idx, eis_idx, y_max + 0.03,  p_ras)
draw_bracket(ax, eis_idx, wt_idx,  y_max + 0.12, p_wt)

for i, grp in enumerate(GROUP_ORDER_BOX):
    n = (plot_data["group"] == grp).sum()
    ax.text(i, -0.07, f"n={n}", ha="center", va="top",
            fontsize=9, transform=ax.get_xaxis_transform())

n_high = (eisfeld_scores > 0.5).sum()
ax.text(eis_idx, y_max + 0.22,
        f"{n_high}/{len(eisfeld_scores)}\n>0.5",
        ha="center", va="bottom", fontsize=8,
        color=PALETTE["Eisfeld ARHG-mutant"], fontweight="bold")

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5)
ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score")
ax.set_title(f"RAS Activation: Eisfeld ARHG-mutant vs BeatAML Groups\n({NORM_LABEL}; dashed = 0.5 threshold)")
ax.set_ylim(bottom=-0.05, top=y_max + 0.38)

fig.tight_layout()
out = FIG_DIR / "cohort_comparison_boxplot_logcpm.png"
fig.savefig(out, dpi=150)
plt.close(fig)
print(f"Saved: {out.name}")

# ===========================================================================
# FIGURE 2 — RAS co-mutation validation violin
# ===========================================================================
# Groups:
#   RAS-mutant      — BeatAML
#   ARHG-only       — Eisfeld ARHG-mutant WITHOUT RAS co-mutation
#   Eisfeld (all)   — Eisfeld ARHG-mutant (all, with co-muts highlighted)
#   Double-WT       — BeatAML

eisfeld["ras_comut"] = eisfeld["sample_id"].isin(ras_comut_ids)

arhg_only   = eisfeld[~eisfeld["ras_comut"]].copy()
arhg_only["group"] = "ARHG-only\n(no RAS co-mut)"
arhg_all    = eisfeld.copy()
arhg_all["group"] = "ARHG-mutant"

ras_mut = beataml[beataml["group"] == "RAS-mutant"].copy()
dbl_wt  = beataml[beataml["group"] == "Double-WT"].copy()

VPAL = {
    "RAS-mutant"           : "#E45756",
    "ARHG-only\n(no RAS co-mut)": "#F28E2B",
    "ARHG-mutant"          : "#76B7B2",
    "Double-WT"            : "#BAB0AC",
}
VORDER = ["RAS-mutant", "ARHG-only\n(no RAS co-mut)", "ARHG-mutant", "Double-WT"]

vdata = pd.concat([
    ras_mut[["sample_id", "RAS_score"]].assign(group="RAS-mutant"),
    arhg_only[["sample_id", "RAS_score", "group"]],
    arhg_all[["sample_id", "RAS_score", "group"]],
    dbl_wt[["sample_id", "RAS_score"]].assign(group="Double-WT"),
], ignore_index=True)

# Bonferroni-corrected pairwise tests for ARHG-only vs others
arhg_only_scores = arhg_only["RAS_score"].values
comparisons_v = [
    ("RAS-mutant", ras_mut["RAS_score"].values),
    ("Double-WT",  dbl_wt["RAS_score"].values),
]
raw_ps = [mannwhitneyu(arhg_only_scores, other, alternative="two-sided")[1]
          for _, other in comparisons_v]
_, adj_ps, _, _ = multipletests(raw_ps, method="bonferroni")
p_vs_wt_bonf = adj_ps[1]

print(f"\nARHG-only vs Double-WT (Bonferroni): p={p_vs_wt_bonf:.4f}")

# Build highlight dataframe — co-mutant patients with their RAS gene label
highlight_rows = []
for _, row in exclusions[exclusions["in_score_csv"] == True].drop_duplicates("ras_score_id").iterrows():
    sid = str(row["ras_score_id"]).replace("_blast", "")
    score_match = eisfeld.loc[eisfeld["sample_id"] == sid, "RAS_score"]
    if len(score_match) == 0:
        continue
    highlight_rows.append({
        "sample_id": sid,
        "RAS_score": score_match.iloc[0],
        "label": row["gene"],
    })
highlights = pd.DataFrame(highlight_rows)

fig, ax = plt.subplots(figsize=(7, 5.5))

sns.violinplot(
    data=vdata, x="group", y="RAS_score",
    order=VORDER, palette=VPAL,
    inner=None, linewidth=0.8, alpha=0.55, ax=ax,
)

# Highlight RAS co-mutant dots within ARHG-mutant column
arhg_mut_x = VORDER.index("ARHG-mutant")
for _, hr in highlights.iterrows():
    ax.scatter(arhg_mut_x, hr["RAS_score"],
               color="#6B0000", s=55, zorder=5, edgecolors="white", linewidths=0.4)

# Label co-mutant dots
if len(highlights) > 0:
    for _, hr in highlights.sort_values("RAS_score", ascending=False).iterrows():
        ax.text(arhg_mut_x + 0.12, hr["RAS_score"], hr["label"],
                va="center", ha="left", fontsize=7.5, color="#6B0000")

# Bonferroni p annotation on ARHG-only column
arhg_only_x = VORDER.index("ARHG-only\n(no RAS co-mut)")
stars_bonf = "***" if p_vs_wt_bonf < 0.001 else "**" if p_vs_wt_bonf < 0.01 else "*" if p_vs_wt_bonf < 0.05 else "ns"
ax.text(arhg_only_x, vdata["RAS_score"].max() + 0.02,
        f"Bonferroni p={p_vs_wt_bonf:.4f}{stars_bonf}",
        ha="center", va="bottom", fontsize=8, color="grey")

ax.set_xlabel("")
ax.set_ylabel("RAS activation score")
ax.set_title(f"RAS pathway activation: ARHG-mutant AML\nRAS co-mutation patients highlighted ({NORM_LABEL})")
ax.set_ylim(top=vdata["RAS_score"].max() + 0.12)

fig.tight_layout()
for ext in ("png", "pdf"):
    out = FIG_DIR / f"ras_comutat_validation_logcpm.{ext}"
    fig.savefig(out, dpi=150 if ext == "png" else None, bbox_inches="tight")
    print(f"Saved: {out.name}")
plt.close(fig)

# ===========================================================================
# FIGURE 3 — RAS score by ARHG gene (dot plot, RAS-clean only)
# ===========================================================================
scores_map = eisfeld_clean.set_index("sample_id")["RAS_score"].to_dict()

variants = pd.read_csv(VARIANTS, sep="\t")
variants["sample"] = variants["Sample_ID"].str.replace("_blast", "", regex=False)
variants = variants[["Hugo_Symbol", "sample"]].rename(columns={"Hugo_Symbol": "gene"})

df = variants.merge(
    pd.Series(scores_map, name="RAS_score").reset_index().rename(columns={"index": "sample"}),
    on="sample", how="inner"
)

def classify(gene):
    if gene == "ARHGDIG": return "GDI"
    if "GEF" in gene:     return "GEF"
    return "GAP"

df["class"] = df["gene"].apply(classify)

gene_order = (
    df.groupby("gene")["RAS_score"]
    .mean()
    .sort_values(ascending=True)
    .index.tolist()
)

DPAL = {"GEF": "#E07B39", "GAP": "#4A90A4", "GDI": "#7B5EA7"}

fig, ax = plt.subplots(figsize=(6, 7))

for gene in gene_order:
    sub = df[df["gene"] == gene]
    y = gene_order.index(gene)
    for _, row in sub.iterrows():
        ax.scatter(row["RAS_score"], y,
                   color=DPAL[row["class"]], s=70, zorder=3,
                   edgecolors="white", linewidths=0.5)

ax.axvline(0.5, color="black", linestyle="--", linewidth=0.9, alpha=0.6)
ax.set_yticks(range(len(gene_order)))
ax.set_yticklabels([f"$\\it{{{g}}}$" for g in gene_order], fontsize=9)

for tick, gene in zip(ax.get_yticklabels(), gene_order):
    tick.set_color(DPAL[df[df["gene"] == gene]["class"].iloc[0]])

ax.set_xlabel("RAS Activation Score", fontsize=10)
ax.set_xlim(-0.02, 1.08)
ax.set_ylim(-0.7, len(gene_order) - 0.3)
ax.set_title(f"RAS Activation Score by ARHG Gene\n(Eisfeld cohort; {NORM_LABEL}; RAS-clean)", fontsize=11)

legend_handles = [
    mlines.Line2D([], [], color=DPAL["GEF"], marker="o", linestyle="None", markersize=7, label="GEF"),
    mlines.Line2D([], [], color=DPAL["GAP"], marker="o", linestyle="None", markersize=7, label="GAP"),
    mlines.Line2D([], [], color=DPAL["GDI"], marker="o", linestyle="None", markersize=7, label="GDI"),
]
ax.legend(handles=legend_handles, fontsize=8, loc="lower right", framealpha=0.8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.6)

fig.tight_layout()
out = FIG_DIR / "ras_score_by_gene_logcpm.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out.name}")

# ===========================================================================
# FIGURE 4 — RAS score by gene family: ARHGAP vs ARHGEF (RAS-clean)
# ===========================================================================
variants_fam = pd.read_csv(VARIANTS, sep="\t")
variants_fam["sample"] = variants_fam["Sample_ID"].str.replace("_blast", "", regex=False)

def classify_family(gene):
    if "GEF" in gene: return "GEF"
    if "GDI" in gene: return "GDI"
    return "GAP"

variants_fam["class"] = variants_fam["Hugo_Symbol"].apply(classify_family)
by_sample_fam = variants_fam.groupby("sample")["class"].apply(set)

def assign_family(classes):
    if "GAP" in classes and "GEF" in classes: return "Both"
    if "GEF" in classes: return "GEF"
    if "GAP" in classes: return "GAP"
    return "GDI"

sample_family = by_sample_fam.apply(assign_family).reset_index()
sample_family.columns = ["sample_id", "family"]

fam_df = eisfeld_clean[["sample_id", "RAS_score"]].merge(sample_family, on="sample_id", how="left")

# Restrict to GAP vs GEF only (exclude 1 Both + 1 GDI)
fam_plot = fam_df[fam_df["family"].isin(["GAP", "GEF"])].copy()
gap_scores = fam_plot.loc[fam_plot["family"] == "GAP", "RAS_score"].values
gef_scores = fam_plot.loc[fam_plot["family"] == "GEF", "RAS_score"].values

_, p_fam = mannwhitneyu(gap_scores, gef_scores, alternative="two-sided")
print(f"\nGAP vs GEF Mann-Whitney U: p={p_fam:.4f}  (n_GAP={len(gap_scores)}, n_GEF={len(gef_scores)})")

FAM_PAL = {"GAP": "#4A90A4", "GEF": "#E07B39"}
FAM_ORDER = ["GAP", "GEF"]

fig, ax = plt.subplots(figsize=(4, 5))

sns.boxplot(
    data=fam_plot, x="family", y="RAS_score",
    order=FAM_ORDER, palette=FAM_PAL,
    width=0.45, linewidth=1.2, fliersize=0, ax=ax,
)
sns.stripplot(
    data=fam_plot, x="family", y="RAS_score",
    order=FAM_ORDER, palette=FAM_PAL,
    size=6, alpha=0.75, jitter=True, ax=ax,
)

y_top = fam_plot["RAS_score"].max()
stars_fam = "***" if p_fam < 0.001 else "**" if p_fam < 0.01 else "*" if p_fam < 0.05 else "ns"
draw_bracket(ax, 0, 1, y_top + 0.04, p_fam)

for i, fam in enumerate(FAM_ORDER):
    n = (fam_plot["family"] == fam).sum()
    ax.text(i, -0.07, f"n={n}", ha="center", va="top",
            fontsize=9, transform=ax.get_xaxis_transform())

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5)
ax.set_xlabel("ARHG Gene Family")
ax.set_ylabel("RAS Activation Score")
ax.set_title(f"RAS Score: ARHGAP vs ARHGEF\n({NORM_LABEL}; RAS-clean; excludes 1 GAP+GEF, 1 GDI)")
ax.set_ylim(bottom=-0.05, top=y_top + 0.22)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
out = FIG_DIR / "ras_score_gap_vs_gef_logcpm.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out.name}")

# ===========================================================================
# FIGURE 5 — 4-group comparison: RAS-mutant, ARHGAP, ARHGEF, Double-WT
# ===========================================================================
gap_df = fam_plot[fam_plot["family"] == "GAP"][["sample_id", "RAS_score"]].assign(group="ARHGAP-mutant")
gef_df = fam_plot[fam_plot["family"] == "GEF"][["sample_id", "RAS_score"]].assign(group="ARHGEF-mutant")
ras_df  = beataml[beataml["group"] == "RAS-mutant"][["sample_id", "RAS_score"]].assign(group="RAS-mutant")
wt_df   = beataml[beataml["group"] == "Double-WT"][["sample_id", "RAS_score"]].assign(group="Double-WT")

four_groups = pd.concat([ras_df, gap_df, gef_df, wt_df], ignore_index=True)

G4_ORDER = ["RAS-mutant", "ARHGAP-mutant", "ARHGEF-mutant", "Double-WT"]
G4_PAL   = {
    "RAS-mutant"   : "#E45756",
    "ARHGAP-mutant": "#4A90A4",
    "ARHGEF-mutant": "#E07B39",
    "Double-WT"    : "#BAB0AC",
}

wt_scores = wt_df["RAS_score"].values
comparisons_4 = [
    ("RAS-mutant",    ras_df["RAS_score"].values),
    ("ARHGAP-mutant", gap_df["RAS_score"].values),
    ("ARHGEF-mutant", gef_df["RAS_score"].values),
]
raw_ps_4 = [mannwhitneyu(wt_scores, s, alternative="two-sided")[1] for _, s in comparisons_4]
_, adj_ps_4, _, _ = multipletests(raw_ps_4, method="bonferroni")

print("\n4-group Bonferroni-corrected vs Double-WT:")
for (label, _), p_raw, p_adj in zip(comparisons_4, raw_ps_4, adj_ps_4):
    n = (four_groups["group"] == label).sum()
    print(f"  {label:20s}  raw p={p_raw:.4f}  Bonf p={p_adj:.4f}  n={n}")

fig, ax = plt.subplots(figsize=(6, 5))

sns.boxplot(
    data=four_groups, x="group", y="RAS_score",
    order=G4_ORDER, palette=G4_PAL,
    width=0.5, linewidth=1.2, fliersize=0, ax=ax,
    hue="group", legend=False,
)
sns.stripplot(
    data=four_groups, x="group", y="RAS_score",
    order=G4_ORDER, palette=G4_PAL,
    size=3.5, alpha=0.55, jitter=True, ax=ax,
    hue="group", legend=False,
)

y_top = four_groups["RAS_score"].max()
wt_x = G4_ORDER.index("Double-WT")
bracket_y = y_top + 0.04
for label, _, p_adj in zip([c[0] for c in comparisons_4], comparisons_4, adj_ps_4):
    x_pos = G4_ORDER.index(label)
    draw_bracket(ax, x_pos, wt_x, bracket_y, p_adj)
    bracket_y += 0.10

for i, grp in enumerate(G4_ORDER):
    n = (four_groups["group"] == grp).sum()
    ax.text(i, -0.07, f"n={n}", ha="center", va="top",
            fontsize=9, transform=ax.get_xaxis_transform())

ax.axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.5)
ax.set_xlabel("")
ax.set_ylabel("RAS Activation Score")
ax.set_title(f"RAS Activation by Mutation Group\n({NORM_LABEL}; ARHG groups RAS-clean; brackets vs Double-WT Bonf.)")
ax.set_ylim(bottom=-0.05, top=bracket_y + 0.05)
ax.set_xticks(range(len(G4_ORDER)))
ax.set_xticklabels(G4_ORDER, rotation=15, ha="right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
out = FIG_DIR / "ras_score_4groups_logcpm.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out.name}")

print("\nDone. All log-CPM figures saved. Original figures untouched.")
