library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# ── Configuration ──────────────────────────────────────────────────────────────
# VAF_THRESHOLD: minimum T-N VAF to count a variant as "present" at a timepoint.
# A variant with VAF < threshold at one timepoint is treated as undetected (VAF
# filled to 0 from pivot_wider), which drives the Gained/Lost classification.
VAF_THRESHOLD <- 0.02

# Status color palette — canonical source: oncoprint_dx_relapse.R
COLORS <- c(
  Persistent = "#2166AC",   # blue   — detected at both Dx and Relapse
  Lost       = "#D6604D",   # coral  — detected at Dx only
  Gained     = "#4DAC26"    # green  — detected at Relapse only
)

# Artifact / low-complexity genes excluded from all Eisfeld Lab analyses.
# List sourced from combine_pass_variants.R via oncoprint_dx_relapse.R.
artifact_genes <- c(
  "ABCA13", "ABCA4", "AHNAK", "AHNAK2", "ANKRD30A", "ANKRD30B", "BSN", "CNTNAP2",
  "COL1A2", "COL6A3", "DNAH1", "DNAH10", "DNAH11", "DNAH17", "DNAH2", "DNAH5",
  "DNAH8", "DNAH9", "FAT1", "FAT2", "FLG", "FRAS1", "FREM2", "HRNR", "KLRC1",
  "KMT2C", "KRTAP10-1", "KRTAP10-11", "KRTAP10-4", "KRTAP1-1", "KRTAP20-4",
  "LAMA1", "LAMA5", "LRP1", "LRP1B", "MAGEB1", "MKI67", "MST1", "MUC12", "MUC16",
  "MUC17", "MUC19", "MUC2", "MUC21", "MUC22", "MUC3A", "MUC4", "MUC5B", "MUC7",
  "MYH2", "MYH3", "MYH4", "MYO15A", "NBPF10", "NBPF11", "NBPF12", "NBPF8", "NBPF9",
  "NEB", "OR10H1", "OR10H5", "OR10Z1", "OR11G2", "OR13C9", "OR13G1", "OR1A1", "OR1L4",
  "OR1Q1", "OR2G3", "OR2H1", "OR2L8", "OR2T6", "OR2T7", "OR2T8", "OR4C11", "OR4C13",
  "OR4C15", "OR4D1", "OR4D11", "OR4F5", "OR51A7", "OR51G2", "OR52A1", "OR52D1",
  "OR52E6", "OR52E8", "OR56B1", "OR5D16", "OR5T2", "OR6C2", "OR6C65", "OR7D4",
  "OR8B12", "OR8G5", "OR9A4", "PCDHA1", "PCDHA10", "PCDHA12", "PCDHA4", "PCDHA6",
  "PCDHB16", "PCDHB6", "PCDHGA12", "PCDHGA6", "PCDHGB6", "PCDHGC4", "PCDHGC5",
  "PKHD1L1", "PRAMEF17", "PSG3", "PSG4", "PSG6", "PSG9", "RBMX", "RBMXL1", "RBMXL3",
  "RNF213", "RYR1", "RYR3", "SPTBN5", "SSX5", "STAB2", "SYNE2", "TACC2", "TDRD3",
  "TNXB", "TPTE", "TRIM51", "TTN", "UBR4", "USH2A", "UTRN", "VPS13B", "VPS13D",
  "VWF", "ZFHX3"
)

# ── 1. Load combined variant file ─────────────────────────────────────────────
message("Loading relapse_pass_variants_combined.xlsx ...")
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))
message(sprintf("  %d rows, %d columns", nrow(df), ncol(df)))

# ── 2. Filter: timepoints, artifact genes, hypermutated patients ───────────────
# Patient_23 and Patient_25 are excluded as hypermutators (Groups 11 and 4).
# artifact_genes are recurrent sequencing artifacts unlikely to represent true somatic
# mutations; they inflate burden counts and confound clonal evolution analyses.
df_somatic <- df |>
  filter(
    Timepoint      %in% c("diagnosis", "relapse"),
    !Gene          %in% artifact_genes,
    !Patient_Group %in% c("Patient_23", "Patient_25")
  ) |>
  # Collapse all ARHG family members (ARHGEF*, ARHGAP*) into a single label so
  # the scatter treats the family as one gene entity, consistent with other scripts.
  mutate(Gene = if_else(startsWith(Gene, "ARHG"), "ARHG", Gene))

message(sprintf("  After filters: %d rows", nrow(df_somatic)))

# ── 3. Deduplicate: one row per Pair_ID × Gene × HGVS × Timepoint ─────────────
# When a variant appears multiple times (e.g., from overlapping capture regions or
# duplicate sample entries), take the maximum observed VAF as the most conservative
# estimate of clonal prevalence.
variants_dedup <- df_somatic |>
  group_by(Pair_ID, Gene, HGVS, Timepoint) |>
  summarise(vaf = max(`T-N VAF`, na.rm = TRUE), .groups = "drop")

# ── 4. Pivot wide: one row per Pair_ID × Gene × HGVS ─────────────────────────
# variants absent at a timepoint receive VAF = 0 (undetected, not truly zero),
# which propagates correctly into the threshold-based classification below.
variants_wide <- variants_dedup |>
  pivot_wider(
    names_from  = Timepoint,
    values_from = vaf,
    values_fill = 0
  ) |>
  rename(vaf_dx = diagnosis, vaf_rel = relapse)

# Guard: ensure both columns exist even when all variants are one-sided
if (!"vaf_dx"  %in% names(variants_wide)) variants_wide$vaf_dx  <- 0
if (!"vaf_rel" %in% names(variants_wide)) variants_wide$vaf_rel <- 0

# ── 5. Retain only variants detectable at >= 1 timepoint ──────────────────────
# A row where both VAFs are 0 (or sub-threshold) after pivoting would represent a
# variant that passed other filters but was never robustly detected — exclude these.
variants_wide <- variants_wide |>
  filter(vaf_dx >= VAF_THRESHOLD | vaf_rel >= VAF_THRESHOLD)

message(sprintf("  Variants after threshold filter: %d", nrow(variants_wide)))

# ── 6. Classify clonal trajectory ─────────────────────────────────────────────
# Persistent: VAF >= threshold at both timepoints — clone survived into relapse
# Lost:       VAF >= threshold at Dx only        — clone eliminated / below detection at relapse
# Gained:     VAF >= threshold at Relapse only   — de novo or sub-threshold clone that expanded
#
# Rows where neither timepoint meets the threshold after this step are impossible
# (filtered above), but the NA guard is retained for safety.
variants_classified <- variants_wide |>
  mutate(
    present_dx  = (vaf_dx  >= VAF_THRESHOLD),
    present_rel = (vaf_rel >= VAF_THRESHOLD),
    status = case_when(
      present_dx  & present_rel  ~ "Persistent",
      present_dx  & !present_rel ~ "Lost",
      !present_dx & present_rel  ~ "Gained",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(status)) |>
  mutate(status = factor(status, levels = c("Persistent", "Lost", "Gained")))

n_total     <- nrow(variants_classified)
n_persist   <- sum(variants_classified$status == "Persistent")
n_lost      <- sum(variants_classified$status == "Lost")
n_gained    <- sum(variants_classified$status == "Gained")

message(sprintf(
  "  Classified — Persistent: %d | Lost: %d | Gained: %d  (total: %d)",
  n_persist, n_lost, n_gained, n_total
))

# ── 7. Build annotation label for plot ────────────────────────────────────────
# Shown in the upper-left quadrant (Gained region): total variant count plus
# per-status breakdown so the figure is self-contained without a data table.
count_label <- sprintf(
  "n = %d\n  Persistent: %d\n  Lost: %d\n  Gained: %d",
  n_total, n_persist, n_lost, n_gained
)

# ── 7b. Top 10 most changed variants (by absolute VAF shift) ─────────────────
# Label = Gene name. For ties in the same gene, the HGVS suffix distinguishes them.
top_changed <- variants_classified |>
  mutate(abs_delta = abs(vaf_rel - vaf_dx)) |>
  slice_max(abs_delta, n = 20, with_ties = FALSE) |>
  mutate(label = Gene)

# ── 8. Scatter plot ───────────────────────────────────────────────────────────
# Diagonal identity line (y = x): points above the line have higher Relapse VAF
# (clone expanded), points below have lower Relapse VAF (clone contracted or lost).
# geom_rug adds marginal tick marks along both axes to show 1-D density per status
# without obscuring the scatter; kept light (alpha = 0.3) to avoid clutter.
p <- ggplot(variants_classified, aes(x = vaf_dx, y = vaf_rel, color = status)) +

  # Identity line drawn first so points render on top
  geom_abline(
    intercept = 0, slope = 1,
    linetype  = "dashed", color = "grey60", linewidth = 0.5
  ) +

  # Marginal rug: 1-D projection of VAF distributions along each axis
  geom_rug(
    aes(color = status),
    alpha       = 0.3,
    linewidth   = 0.4,
    outside     = FALSE,
    show.legend = FALSE
  ) +

  # Main scatter
  geom_point(size = 1.8, alpha = 0.6) +

  # Label top 10 most changed variants; repel prevents overlap
  geom_text_repel(
    data        = top_changed,
    aes(label   = label),
    size        = 2.8,
    fontface    = "italic",
    color       = "grey20",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey60",
    segment.size  = 0.3,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +

  # Annotate total variant count in upper-left (Gained quadrant)
  annotate(
    "text",
    x     = 0.05,
    y     = 0.92,
    label = count_label,
    hjust = 0,
    vjust = 1,
    size  = 3,
    color = "grey30",
    lineheight = 1.2
  ) +

  scale_color_manual(
    values = COLORS,
    name   = "Clonal trajectory"
  ) +

  # Identical axes: 0–1 with quarter-fraction breaks
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  labs(
    title = "VAF dynamics: Diagnosis \u2192 Relapse",
    x     = "Diagnosis VAF",
    y     = "Relapse VAF"
  ) +

  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "grey40"),
    axis.ticks      = element_line(color = "grey40")
  )

print(p)

# ── 9. Save outputs ───────────────────────────────────────────────────────────
ggsave(
  file.path(BASE_DIR, "Cohorts/Relapse/vaf_scatter_dx_relapse.pdf"),
  plot   = p,
  width  = 5.5,
  height = 5,
  units  = "in"
)

ggsave(
  file.path(BASE_DIR, "Cohorts/Relapse/vaf_scatter_dx_relapse.png"),
  plot   = p,
  width  = 5.5,
  height = 5,
  units  = "in",
  dpi    = 300
)

message("Saved: Cohorts/Relapse/vaf_scatter_dx_relapse.pdf/.png")
