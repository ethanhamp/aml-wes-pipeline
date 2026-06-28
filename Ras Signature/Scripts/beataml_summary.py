"""
beataml_summary.py
==================
Summary visualizations of the BeatAML cohort used to train the RAS signature.

Plots generated:
  1. cohort_overview.png       — total RNA / WES / matched sample counts
  2. ras_mutant_breakdown.png  — RAS-mutant vs WT + per-gene breakdown
  3. ras_vaf_distribution.png  — VAF distribution of RAS pathway mutations
  4. variant_classes.png       — mutation types in RAS pathway genes
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ===========================================================================
# PATHS
# ===========================================================================
DATA_DIR = Path(os.environ.get("BEATAML_DIR", str(Path(os.environ.get("ARHG_BASE_DIR", "..")).parent / "BeatAML")))
OUTPUT_DIR = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Ras Signature"

EXPR_FILE = DATA_DIR / "beataml_waves1to4_norm_exp_dbgap.txt"
MUT_FILE  = DATA_DIR / "beataml_wes_wv1to4_mutations_dbgap.txt"

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25)
DPI = 150

RAS_GENES = [
    "NRAS", "KRAS", "HRAS",
    "PTPN11", "NF1", "CBL",
    "RRAS", "RRAS2",
    "RAF1", "BRAF",
    "MAP2K1", "MAP2K2",
]

# ===========================================================================
# LOAD DATA
# ===========================================================================
print("Loading data ...")

expr_raw = pd.read_csv(EXPR_FILE, sep="\t", index_col=None, low_memory=False)
META_COLS   = ["stable_id", "display_label", "description", "biotype"]
sample_cols = [c for c in expr_raw.columns if c not in META_COLS]

mut = pd.read_csv(MUT_FILE, sep="\t", low_memory=False)
mut["patient_id"] = mut["dbgap_sample_id"].str.replace(r"[RD]$", "", regex=True)

# Patient IDs from RNA and WES
rna_patients = set(s.rstrip("R") for s in sample_cols if s.endswith("R"))
wes_patients = set(mut["patient_id"].unique())
matched      = rna_patients & wes_patients

print(f"  RNA samples   : {len(rna_patients):,}")
print(f"  WES patients  : {len(wes_patients):,}")
print(f"  Matched       : {len(matched):,}")

# RAS labels on matched patients
ras_patients = set(
    mut[(mut["patient_id"].isin(matched)) & (mut["symbol"].isin(RAS_GENES))]["patient_id"]
)
wt_patients  = matched - ras_patients

print(f"  RAS-mutant    : {len(ras_patients):,}  ({100*len(ras_patients)/len(matched):.1f}%)")
print(f"  RAS-WT        : {len(wt_patients):,}  ({100*len(wt_patients)/len(matched):.1f}%)")

# RAS mutations subset
ras_mut = mut[(mut["patient_id"].isin(matched)) & (mut["symbol"].isin(RAS_GENES))].copy()

# ===========================================================================
# PLOT 1 — Cohort overview
# ===========================================================================
print("\nPlot 1: Cohort overview ...")

labels = ["RNA-seq\nsamples", "WES\nsamples", "Matched\n(RNA + WES)", "RAS-mutant", "RAS-WT"]
values = [len(rna_patients), len(wes_patients), len(matched), len(ras_patients), len(wt_patients)]
colors = ["#4C78A8", "#4C78A8", "#4C78A8", "#E45756", "#BAB0AC"]

fig, ax = plt.subplots(figsize=(7, 4))
bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.8, width=0.55)

for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 8,
            f"{val:,}", ha="center", va="bottom", fontsize=10, fontweight="bold")

# Bracket: matched → RAS-mutant + RAS-WT
y_br = max(values[:3]) + 40
ax.annotate("", xy=(3, y_br), xytext=(4, y_br),
            arrowprops=dict(arrowstyle="-", color="grey", lw=1))
ax.plot([2, 2, 3, 4, 4],
        [y_br - 15, y_br, y_br, y_br, y_br - 15],
        color="grey", lw=1)
ax.text(3, y_br + 8, "from matched", ha="center", fontsize=8, color="grey")

ax.set_ylabel("Number of samples")
ax.set_title("BeatAML Cohort Overview\n(Waves 1–4, dbGaP release)")
ax.set_ylim(0, max(values) * 1.18)
ax.yaxis.grid(True, alpha=0.5)
ax.set_axisbelow(True)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "cohort_overview.png", dpi=DPI)
plt.close(fig)
print("  Saved: cohort_overview.png")

# ===========================================================================
# PLOT 2 — RAS-mutant vs WT + per-gene breakdown
# ===========================================================================
print("Plot 2: RAS breakdown ...")

# Per-gene patient counts (unique patients per gene)
gene_counts = (
    ras_mut.groupby("symbol")["patient_id"]
    .nunique()
    .reindex(RAS_GENES, fill_value=0)
    .sort_values(ascending=True)
)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Left: donut
ax = axes[0]
sizes  = [len(ras_patients), len(wt_patients)]
clrs   = ["#E45756", "#BAB0AC"]
wedges, texts = ax.pie(
    sizes, colors=clrs, startangle=90,
    wedgeprops=dict(width=0.55, edgecolor="white", linewidth=1.5),
)
total = len(matched)
ax.text(0, 0.12, f"{len(ras_patients):,}", ha="center", va="center",
        fontsize=20, fontweight="bold", color="#E45756")
ax.text(0, -0.18, f"{100*len(ras_patients)/total:.1f}%\nRAS-mutant",
        ha="center", va="center", fontsize=10, color="#E45756")
legend_patches = [
    mpatches.Patch(color="#E45756", label=f"RAS-mutant  (n={len(ras_patients):,})"),
    mpatches.Patch(color="#BAB0AC", label=f"RAS-WT       (n={len(wt_patients):,})"),
]
ax.legend(handles=legend_patches, loc="lower center", bbox_to_anchor=(0.5, -0.12),
          fontsize=9, frameon=False)
ax.set_title(f"RAS-Mutant vs WT\n(matched cohort, n={total:,})", fontsize=11)

# Right: horizontal bar per gene
ax = axes[1]
gene_colors = ["#E45756" if g in {"NRAS","KRAS","HRAS"} else "#F28E2B"
               for g in gene_counts.index]
bars = ax.barh(gene_counts.index, gene_counts.values,
               color=gene_colors, edgecolor="white", linewidth=0.5)
for bar, val in zip(bars, gene_counts.values):
    if val > 0:
        ax.text(val + 0.8, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=9)

legend_patches2 = [
    mpatches.Patch(color="#E45756", label="Canonical RAS (NRAS/KRAS/HRAS)"),
    mpatches.Patch(color="#F28E2B", label="RAS pathway (upstream/downstream)"),
]
ax.legend(handles=legend_patches2, loc="lower right", fontsize=8, frameon=False)
ax.set_xlabel("Unique patients mutated")
ax.set_title("RAS Pathway Mutations by Gene\n(patients with ≥1 mutation)", fontsize=11)
ax.set_xlim(0, gene_counts.max() * 1.2)
ax.xaxis.grid(True, alpha=0.4)
ax.set_axisbelow(True)

fig.suptitle("BeatAML — RAS Pathway Mutation Landscape", fontsize=12, y=1.01)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "ras_mutant_breakdown.png", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  Saved: ras_mutant_breakdown.png")

# ===========================================================================
# PLOT 3 — VAF distribution of RAS mutations
# ===========================================================================
print("Plot 3: VAF distribution ...")

vaf_data = ras_mut[ras_mut["t_vaf"].notna()].copy()
vaf_data["t_vaf"] = pd.to_numeric(vaf_data["t_vaf"], errors="coerce")
vaf_data = vaf_data[vaf_data["t_vaf"].between(0, 1)]

# Top genes by count for per-gene panel
top_genes = (
    vaf_data.groupby("symbol")["patient_id"].nunique()
    .sort_values(ascending=False).head(6).index.tolist()
)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# Left: overall VAF histogram
ax = axes[0]
ax.hist(vaf_data["t_vaf"], bins=40, color="#E45756", edgecolor="white",
        linewidth=0.5, alpha=0.85)
ax.axvline(vaf_data["t_vaf"].median(), color="black", linestyle="--",
           linewidth=1.5, label=f"Median VAF = {vaf_data['t_vaf'].median():.2f}")
ax.set_xlabel("Tumor VAF")
ax.set_ylabel("Variant count")
ax.set_title("RAS Pathway Mutation VAF\n(all genes, all variants)")
ax.legend(fontsize=9)

# Right: per-gene VAF boxplot (top 6 genes)
ax = axes[1]
vaf_top = vaf_data[vaf_data["symbol"].isin(top_genes)].copy()
gene_order = (
    vaf_top.groupby("symbol")["patient_id"].nunique()
    .sort_values(ascending=False).index.tolist()
)
sns.boxplot(
    data=vaf_top, x="symbol", y="t_vaf",
    order=gene_order, color="#E45756",
    width=0.5, linewidth=1, fliersize=2, ax=ax,
)
ax.set_xlabel("")
ax.set_ylabel("Tumor VAF")
ax.set_title("VAF by Gene (top 6 by patient count)")

# Add n= labels
for i, gene in enumerate(gene_order):
    n = vaf_top[vaf_top["symbol"] == gene]["patient_id"].nunique()
    ax.text(i, -0.06, f"n={n}", ha="center", va="top",
            fontsize=8, transform=ax.get_xaxis_transform())

fig.suptitle("BeatAML — RAS Mutation VAF Distribution", fontsize=12)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "ras_vaf_distribution.png", dpi=DPI)
plt.close(fig)
print("  Saved: ras_vaf_distribution.png")

# ===========================================================================
# PLOT 4 — Variant classification breakdown
# ===========================================================================
print("Plot 4: Variant classification ...")

vc_counts = (
    ras_mut["variant_classification"]
    .value_counts()
    .head(8)
)

VARIANT_COLORS = {
    "Missense_Mutation"      : "#E45756",
    "Nonsense_Mutation"      : "#4C78A8",
    "Frame_Shift_Del"        : "#F28E2B",
    "Frame_Shift_Ins"        : "#FFBF00",
    "Splice_Site"            : "#54A24B",
    "In_Frame_Del"           : "#B279A2",
    "In_Frame_Ins"           : "#9D755D",
    "Translation_Start_Site" : "#BAB0AC",
}
bar_colors = [VARIANT_COLORS.get(v, "#BAB0AC") for v in vc_counts.index]

fig, ax = plt.subplots(figsize=(7, 4))
bars = ax.barh(vc_counts.index[::-1], vc_counts.values[::-1],
               color=bar_colors[::-1], edgecolor="white", linewidth=0.5)
for bar, val in zip(bars, vc_counts.values[::-1]):
    ax.text(val + 0.5, bar.get_y() + bar.get_height() / 2,
            str(val), va="center", fontsize=9)

ax.set_xlabel("Variant count")
ax.set_title("RAS Pathway Mutation Types\n(BeatAML; all RAS pathway genes)")
ax.set_xlim(0, vc_counts.max() * 1.18)
ax.xaxis.grid(True, alpha=0.4)
ax.set_axisbelow(True)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "variant_classes.png", dpi=DPI)
plt.close(fig)
print("  Saved: variant_classes.png")

# ===========================================================================
# SUMMARY
# ===========================================================================
print("\n" + "=" * 55)
print("BEATAML COHORT SUMMARY")
print("=" * 55)
print(f"  RNA-seq samples          : {len(rna_patients):,}")
print(f"  WES samples             : {len(wes_patients):,}")
print(f"  Matched (RNA + WES)      : {len(matched):,}")
print(f"  RAS-mutant               : {len(ras_patients):,}  ({100*len(ras_patients)/len(matched):.1f}%)")
print(f"  RAS-WT                   : {len(wt_patients):,}  ({100*len(wt_patients)/len(matched):.1f}%)")
print()
print("  RAS gene breakdown (unique patients):")
for gene, cnt in gene_counts.sort_values(ascending=False).items():
    if cnt > 0:
        print(f"    {gene:<10} : {cnt:3d}")
print()
print("  Variant classifications (RAS pathway mutations):")
for vc, cnt in vc_counts.items():
    print(f"    {vc:<30} : {cnt:4d}")
print("=" * 55)
print("Done. Plots saved to:", OUTPUT_DIR)
