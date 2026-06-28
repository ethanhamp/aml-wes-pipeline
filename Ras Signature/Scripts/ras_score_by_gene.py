"""
ras_score_by_gene.py
====================
Dot plot of RAS activation score per ARHG gene, colored by functional class
(GEF / GAP / GDI). Genes sorted by mean score descending. One dot per patient.
RAS co-mutant patients excluded (ARHG-only, RAS-clean subset).

Inputs:
  applied_ras_scores.csv            — Eisfeld z-score normalized scores
  ras_pathway_mutants.csv           — RAS co-mutation exclusion list
  ../Variants/arhg_variants_clean.txt — per-patient ARHG gene assignments

Output:
  Figures/ras_score_by_gene.pdf
"""

import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent)))
SCORES    = BASE / "Ras Signature/applied_ras_scores.csv"
EXCLUSIONS = BASE / "Ras Signature/ras_pathway_mutants.csv"
VARIANTS  = BASE / "Variants/arhg_variants_clean.txt"
OUT       = BASE / "Ras Signature/Figures/ras_score_by_gene.pdf"

# ---------------------------------------------------------------------------
# Load & join (RAS-clean only)
# ---------------------------------------------------------------------------
scores = pd.read_csv(SCORES).rename(columns={"sample_id": "sample"})
scores["sample"] = scores["sample"].str.replace("_blast", "", regex=False)

exclusions = pd.read_csv(EXCLUSIONS)
ras_comut_ids = set(
    exclusions.loc[exclusions["in_score_csv"] == True, "ras_score_id"]
    .str.replace("_blast", "", regex=False)
)
scores = scores[~scores["sample"].isin(ras_comut_ids)]
print(f"RAS-clean patients: {len(scores)} (excluded {len(ras_comut_ids)} RAS co-mutant)")

variants = pd.read_csv(VARIANTS, sep="\t")
variants["sample"] = variants["Sample_ID"].str.replace("_blast", "", regex=False)
variants = variants[["Hugo_Symbol", "sample"]].rename(columns={"Hugo_Symbol": "gene"})

df = variants.merge(scores[["sample", "RAS_score"]], on="sample", how="inner")

# Functional class
def classify(gene):
    if gene == "ARHGDIG":
        return "GDI"
    if "GEF" in gene:
        return "GEF"
    return "GAP"

df["class"] = df["gene"].apply(classify)

# Sort genes by mean score descending
gene_order = (
    df.groupby("gene")["RAS_score"]
    .mean()
    .sort_values(ascending=True)  # ascending=True → highest at top in horizontal plot
    .index.tolist()
)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
PALETTE = {"GEF": "#E07B39", "GAP": "#4A90A4", "GDI": "#7B5EA7"}

fig, ax = plt.subplots(figsize=(6, 7))

for gene in gene_order:
    sub = df[df["gene"] == gene]
    y = gene_order.index(gene)
    for _, row in sub.iterrows():
        ax.scatter(
            row["RAS_score"], y,
            color=PALETTE[row["class"]],
            s=70, zorder=3, edgecolors="white", linewidths=0.5,
        )

ax.axvline(0.5, color="black", linestyle="--", linewidth=0.9, alpha=0.6,
           label="Score = 0.5")

ax.set_yticks(range(len(gene_order)))
ax.set_yticklabels(
    [f"$\\it{{{g}}}$" for g in gene_order],
    fontsize=9,
)

# Color gene labels by class
for tick, gene in zip(ax.get_yticklabels(), gene_order):
    cls = df[df["gene"] == gene]["class"].iloc[0]
    tick.set_color(PALETTE[cls])

ax.set_xlabel("RAS Activation Score", fontsize=10)
ax.set_xlim(-0.02, 1.08)
ax.set_ylim(-0.7, len(gene_order) - 0.3)
ax.set_title("RAS Activation Score by ARHG Gene\n(Eisfeld cohort; z-score normalized; RAS-clean)", fontsize=11)

# Legend
legend_handles = [
    mlines.Line2D([], [], color=PALETTE["GEF"], marker="o", linestyle="None",
                  markersize=7, label="GEF"),
    mlines.Line2D([], [], color=PALETTE["GAP"], marker="o", linestyle="None",
                  markersize=7, label="GAP"),
    mlines.Line2D([], [], color=PALETTE["GDI"], marker="o", linestyle="None",
                  markersize=7, label="GDI"),
]
ax.legend(handles=legend_handles, fontsize=8, loc="lower right", framealpha=0.8)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.6)

fig.tight_layout()
fig.savefig(OUT, bbox_inches="tight")
print(f"Saved: {OUT}")
