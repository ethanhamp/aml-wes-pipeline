library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# ── Configuration ──────────────────────────────────────────────────────────────
VAF_THRESHOLD <- 0.02

COLORS <- c(
  Persistent = "#2166AC",
  Lost       = "#D6604D",
  Gained     = "#4DAC26"
)

local_env <- new.env()
source(file.path(BASE_DIR, "General/jeff.genes.txt"), local = local_env)
jeff.genes <- local_env$jeff.genes

# ── 1. Load ────────────────────────────────────────────────────────────────────
message("Loading relapse_pass_variants_combined.xlsx ...")
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))
message(sprintf("  %d rows, %d columns", nrow(df), ncol(df)))

# ── 2. Filter and assign family ────────────────────────────────────────────────
# Standardise compound gene names (e.g., "ARHGAP11A;ARHGAP11A-SCG5") to the
# primary gene symbol before family assignment.
df_arhg <- df |>
  filter(
    Timepoint      %in% c("diagnosis", "relapse"),
    !Gene          %in% jeff.genes,
    !Patient_Group %in% c("Patient_23", "Patient_25")
  ) |>
  mutate(
    # Strip compound annotations (take first token before ";" or "-SCG")
    Gene_clean = sub(";.*", "", Gene),
    Gene_clean = sub("-SCG.*", "", Gene_clean)
  ) |>
  filter(startsWith(Gene_clean, "ARHGEF") | startsWith(Gene_clean, "ARHGAP")) |>
  mutate(
    family = if_else(startsWith(Gene_clean, "ARHGEF"), "GEF (ARHGEF)", "GAP (ARHGAP)"),
    family = factor(family, levels = c("GEF (ARHGEF)", "GAP (ARHGAP)"))
  )

message(sprintf("  ARHGEF/ARHGAP variants after filters: %d rows", nrow(df_arhg)))
message(sprintf("  GEF: %d rows, GAP: %d rows",
  sum(df_arhg$family == "GEF (ARHGEF)"),
  sum(df_arhg$family == "GAP (ARHGAP)")
))

# ── 3. Deduplicate ─────────────────────────────────────────────────────────────
variants_dedup <- df_arhg |>
  group_by(Pair_ID, Gene_clean, HGVS, Timepoint, family) |>
  summarise(vaf = max(`T-N VAF`, na.rm = TRUE), .groups = "drop")

# ── 4. Pivot wide ──────────────────────────────────────────────────────────────
variants_wide <- variants_dedup |>
  pivot_wider(
    id_cols     = c(Pair_ID, Gene_clean, HGVS, family),
    names_from  = Timepoint,
    values_from = vaf,
    values_fill = 0
  )

if (!"diagnosis" %in% names(variants_wide)) variants_wide$diagnosis <- 0
if (!"relapse"   %in% names(variants_wide)) variants_wide$relapse   <- 0

variants_wide <- variants_wide |>
  rename(vaf_dx = diagnosis, vaf_rel = relapse) |>
  filter(vaf_dx >= VAF_THRESHOLD | vaf_rel >= VAF_THRESHOLD)

message(sprintf("  Variants after threshold filter: %d", nrow(variants_wide)))

# ── 5. Classify clonal trajectory ─────────────────────────────────────────────
variants_classified <- variants_wide |>
  mutate(
    present_dx  = vaf_dx  >= VAF_THRESHOLD,
    present_rel = vaf_rel >= VAF_THRESHOLD,
    status = case_when(
      present_dx  & present_rel  ~ "Persistent",
      present_dx  & !present_rel ~ "Lost",
      !present_dx & present_rel  ~ "Gained",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(status)) |>
  mutate(status = factor(status, levels = c("Persistent", "Lost", "Gained")))

# ── 6. Summary counts ──────────────────────────────────────────────────────────
summary_tbl <- variants_classified |>
  count(family, status) |>
  group_by(family) |>
  mutate(
    total = sum(n),
    pct   = round(100 * n / total, 1)
  ) |>
  ungroup()

message("\nPersistent/Lost/Gained summary by family:")
print(as.data.frame(summary_tbl))

# ── 7. Fisher's exact test ─────────────────────────────────────────────────────
# Test whether GEF and GAP families differ in status composition (3-level).
# Contingency table: family (rows) x status (cols)
contingency <- variants_classified |>
  count(family, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0) |>
  tibble::column_to_rownames("family")

# Ensure all three columns present
for (col in c("Persistent", "Lost", "Gained")) {
  if (!col %in% names(contingency)) contingency[[col]] <- 0L
}
contingency <- contingency[, c("Persistent", "Lost", "Gained")]

message("\nContingency table (family x status):")
print(contingency)

fisher_result <- fisher.test(contingency, simulate.p.value = TRUE, B = 1e5)
message(sprintf("\nFisher's exact test (GEF vs GAP status composition):\n  p = %.4f", fisher_result$p.value))

# Pairwise: Persistent rate GEF vs GAP
gef_p <- contingency["GEF (ARHGEF)", "Persistent"]; gef_tot <- sum(contingency["GEF (ARHGEF)", ])
gap_p <- contingency["GAP (ARHGAP)", "Persistent"]; gap_tot <- sum(contingency["GAP (ARHGAP)", ])
message(sprintf(
  "  GEF Persistent: %d/%d = %.1f%%  |  GAP Persistent: %d/%d = %.1f%%",
  gef_p, gef_tot, 100*gef_p/gef_tot, gap_p, gap_tot, 100*gap_p/gap_tot
))

fisher_persist <- fisher.test(matrix(
  c(gef_p, gef_tot - gef_p, gap_p, gap_tot - gap_p),
  nrow = 2, dimnames = list(c("GEF", "GAP"), c("Persistent", "Not Persistent"))
))
message(sprintf("  Persistent rate GEF vs GAP Fisher's p = %.4f (OR = %.2f)",
  fisher_persist$p.value, fisher_persist$estimate))

# ── 8. VAF scatter — faceted by family ────────────────────────────────────────
top_changed <- variants_classified |>
  mutate(abs_delta = abs(vaf_rel - vaf_dx)) |>
  group_by(family) |>
  slice_max(abs_delta, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(label = Gene_clean)

# Build per-facet annotation
facet_labels <- variants_classified |>
  count(family, status) |>
  group_by(family) |>
  summarise(
    label = sprintf(
      "n = %d\n  Persistent: %d\n  Lost: %d\n  Gained: %d",
      sum(n),
      n[status == "Persistent"],
      n[status == "Lost"],
      n[status == "Gained"]
    ),
    .groups = "drop"
  )

p_scatter <- ggplot(variants_classified, aes(x = vaf_dx, y = vaf_rel, color = status)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_rug(aes(color = status), alpha = 0.3, linewidth = 0.4, outside = FALSE, show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.6) +
  geom_text_repel(
    data = top_changed,
    aes(label = label),
    size = 2.5, fontface = "italic", color = "grey20",
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "grey60", segment.size = 0.3,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  geom_text(
    data = facet_labels,
    aes(x = 0.04, y = 0.96, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2.8, color = "grey30", lineheight = 1.2
  ) +
  scale_color_manual(values = COLORS, name = "Clonal trajectory") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = expansion(mult = c(0.01, 0.02))) +
  facet_wrap(~ family, ncol = 2) +
  labs(
    title = "VAF dynamics: Diagnosis → Relapse by ARHG family",
    x     = "Diagnosis VAF",
    y     = "Relapse VAF"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    strip.text      = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "grey40"),
    axis.ticks      = element_line(color = "grey40")
  )

# ── 9. Stacked proportion bar ──────────────────────────────────────────────────
p_bar <- ggplot(summary_tbl, aes(x = family, y = pct, fill = status)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
    position = position_stack(vjust = 0.5),
    size = 3.2, color = "white", fontface = "bold", lineheight = 1.1
  ) +
  scale_fill_manual(values = COLORS, name = "Clonal trajectory") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Persistent / Lost / Gained proportions by ARHG family",
    x     = NULL,
    y     = "Percentage of variants (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.line.x     = element_blank(),
    axis.ticks.x    = element_blank()
  ) +
  annotate(
    "text",
    x = 1.5, y = 105,
    label = sprintf("Fisher's p = %.3f", fisher_result$p.value),
    size = 3.2, color = "grey30"
  )

# ── 10. Save outputs ───────────────────────────────────────────────────────────
dir.create(file.path(BASE_DIR, "Cohorts/Relapse/Figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/vaf_scatter_gef_vs_gap.pdf"), plot = p_scatter,
  width = 9, height = 5, units = "in")
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/vaf_scatter_gef_vs_gap.png"), plot = p_scatter,
  width = 9, height = 5, units = "in", dpi = 300)

ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/gainloss_bar_gef_vs_gap.pdf"), plot = p_bar,
  width = 5, height = 4.5, units = "in")
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/gainloss_bar_gef_vs_gap.png"), plot = p_bar,
  width = 5, height = 4.5, units = "in", dpi = 300)

# Save summary table
write.csv(summary_tbl, file.path(BASE_DIR, "Cohorts/Relapse/Figures/gainloss_summary_gef_vs_gap.csv"),
  row.names = FALSE)

message("\nSaved:")
message("  Cohorts/Relapse/Figures/vaf_scatter_gef_vs_gap.pdf/.png")
message("  Cohorts/Relapse/Figures/gainloss_bar_gef_vs_gap.pdf/.png")
message("  Cohorts/Relapse/Figures/gainloss_summary_gef_vs_gap.csv")

message(sprintf("\n[FISHER SUMMARY] GEF vs GAP 3-level status: p = %.4f", fisher_result$p.value))
message(sprintf("[FISHER SUMMARY] Persistent rate GEF vs GAP: p = %.4f  OR = %.2f",
  fisher_persist$p.value, fisher_persist$estimate))
