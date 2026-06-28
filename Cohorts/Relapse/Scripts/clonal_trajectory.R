library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(writexl)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# ── Configuration ──────────────────────────────────────────────────────────────
VAF_THRESHOLD  <- 0.02   # minimum T-N VAF to consider a variant "present"
EXPAND_DELTA   <- 0.10   # VAF increase threshold to call "Expanded"
CONTRACT_DELTA <- 0.10   # VAF decrease threshold to call "Contracted"

COLORS <- c(
  Gained     = "#4DAC26",   # green  — acquired at relapse
  Expanded   = "#B2182B",   # red    — clonal expansion
  Persistent = "#2166AC",   # blue   — stable
  Contracted = "#92C5DE",   # light blue — reduced
  Lost       = "#D6604D"    # coral  — eliminated
)

BENCHMARK_GENES <- c("FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2")

arhg_genes <- c(
  "ARHGEF1","ARHGEF2","ARHGEF3","ARHGEF4","ARHGEF5","ARHGEF6","ARHGEF7",
  "ARHGEF8","ARHGEF9","ARHGEF10","ARHGEF10L","ARHGEF11","ARHGEF12","ARHGEF15",
  "ARHGEF16","ARHGEF17","ARHGEF18","ARHGEF19","ARHGEF25","ARHGEF26","ARHGEF28",
  "ARHGAP1","ARHGAP2","ARHGAP3","ARHGAP4","ARHGAP5","ARHGAP6","ARHGAP7",
  "ARHGAP8","ARHGAP9","ARHGAP10","ARHGAP11A","ARHGAP11B","ARHGAP12","ARHGAP13",
  "ARHGAP14","ARHGAP15","ARHGAP16","ARHGAP17","ARHGAP18","ARHGAP19","ARHGAP20",
  "ARHGAP21","ARHGAP22","ARHGAP23","ARHGAP24","ARHGAP25","ARHGAP26","ARHGAP27",
  "ARHGAP28","ARHGAP29","ARHGAP30","ARHGAP31","ARHGAP32","ARHGAP33","ARHGAP34",
  "ARHGAP35","ARHGAP36","ARHGAP37","ARHGAP38","ARHGAP39","ARHGAP40","ARHGAP41",
  "ARHGAP42","ARHGAP43","ARHGAP44","ARHGAP45"
)

local_env <- new.env()
source(file.path(BASE_DIR, "General/jeff.genes.txt"), local = local_env)
jeff.genes <- local_env$jeff.genes

# ── Output dirs ────────────────────────────────────────────────────────────────
dir.create(file.path(BASE_DIR, "Cohorts/Relapse/Figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BASE_DIR, "Cohorts/Relapse/Output"),  showWarnings = FALSE, recursive = TRUE)

# ── 1. Load and filter ─────────────────────────────────────────────────────────
message("Loading relapse_pass_variants_combined.xlsx ...")
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))
message(sprintf("  %d rows loaded", nrow(df)))

somatic <- df |>
  filter(
    Timepoint      %in% c("diagnosis", "relapse"),
    `T-N VAF`      >= VAF_THRESHOLD,
    `AD Total`     >  20,
    !Gene          %in% jeff.genes,
    !Patient_Group %in% c("Patient_23", "Patient_25")
  )

n_patients <- n_distinct(somatic$Pair_ID)
message(sprintf("  After filters: %d rows, %d patients", nrow(somatic), n_patients))

# ── 2. Dedup and pivot wide ───────────────────────────────────────────────────
variants_dedup <- somatic |>
  group_by(Pair_ID, Patient_Group, Gene, `MANE HGVS`, Effect) |>
  summarise(
    vaf_dx  = { v <- `T-N VAF`[Timepoint == "diagnosis"]; if (length(v) == 0 || all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) },
    vaf_rel = { v <- `T-N VAF`[Timepoint == "relapse"];   if (length(v) == 0 || all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) },
    .groups = "drop"
  ) |>
  mutate(
    vaf_dx  = ifelse(!is.finite(vaf_dx),  NA_real_, vaf_dx),
    vaf_rel = ifelse(!is.finite(vaf_rel), NA_real_, vaf_rel),
    # Treat NA (not observed) as below threshold
    present_dx  = !is.na(vaf_dx)  & vaf_dx  >= VAF_THRESHOLD,
    present_rel = !is.na(vaf_rel) & vaf_rel >= VAF_THRESHOLD
  ) |>
  filter(present_dx | present_rel)   # keep only variants observed at ≥1 timepoint

# ── 3. Classify trajectory ────────────────────────────────────────────────────
# 5-class scheme (GOALS.md specifies 4; "Contracted" added to handle the
# both-present/VAF-decrease case that would otherwise be unclassified):
#   Gained     — Rel only
#   Lost       — Dx only
#   Expanded   — both, delta >= +10 pp
#   Contracted — both, delta <= −10 pp
#   Persistent — both, |delta| < 10 pp
variants_classified <- variants_dedup |>
  mutate(
    delta = coalesce(vaf_rel, 0) - coalesce(vaf_dx, 0),
    trajectory = case_when(
      !present_dx &  present_rel                  ~ "Gained",
       present_dx & !present_rel                  ~ "Lost",
       present_dx &  present_rel & delta >=  EXPAND_DELTA   ~ "Expanded",
       present_dx &  present_rel & delta <= -CONTRACT_DELTA ~ "Contracted",
       present_dx &  present_rel                  ~ "Persistent",
      TRUE ~ NA_character_
    ),
    trajectory = factor(trajectory,
                        levels = c("Gained","Expanded","Persistent","Contracted","Lost")),
    is_arhg = Gene %in% arhg_genes,
    gene_group = case_when(
      Gene %in% arhg_genes      ~ "ARHG",
      Gene %in% BENCHMARK_GENES ~ Gene,
      TRUE                      ~ "Genome-wide"
    ),
    clonal_dx = case_when(
       present_dx & vaf_dx > 0.30 ~ "Clonal",
       present_dx                 ~ "Subclonal",
      !present_dx                 ~ NA_character_
    )
  )

n_total <- nrow(variants_classified)
message("\n=== Trajectory classification (all variants, n=", n_total, ") ===")
print(table(variants_classified$trajectory))

# ── 4. ARHG-specific summary ──────────────────────────────────────────────────
arhg_variants <- variants_classified |> filter(is_arhg)

message("\n=== ARHG variants ===")
message(sprintf("  n = %d variants across %d patients",
                nrow(arhg_variants), n_distinct(arhg_variants$Pair_ID)))
print(table(arhg_variants$trajectory))
print(arhg_variants |> select(Pair_ID, Gene, `MANE HGVS`, vaf_dx, vaf_rel, delta, trajectory))

# ── 5. Save per-variant trajectory CSV ────────────────────────────────────────
write.csv(variants_classified, file.path(BASE_DIR, "Cohorts/Relapse/Output/arhg_trajectory_classified.csv"),
          row.names = FALSE)
message("\nSaved: Cohorts/Relapse/Output/arhg_trajectory_classified.csv")

# ── 6. Trajectory proportion comparison ───────────────────────────────────────
# Groups: ARHG, Genome-wide (non-ARHG, non-benchmark), each benchmark gene
prop_df <- variants_classified |>
  filter(!is.na(trajectory)) |>
  mutate(group = gene_group) |>
  group_by(group, trajectory) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(group) |>
  mutate(
    total = sum(n),
    pct   = n / total * 100
  ) |>
  ungroup() |>
  mutate(
    group = factor(group,
                   levels = c("ARHG", "FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "Genome-wide"))
  )

# Print proportion table
message("\n=== Trajectory proportions by group ===")
print(prop_df |> select(group, trajectory, n, total, pct) |> arrange(group, trajectory))

# Chi-square: ARHG vs genome-wide
arhg_counts    <- variants_classified |> filter(is_arhg, !is.na(trajectory)) |> pull(trajectory)
genome_counts  <- variants_classified |> filter(!is_arhg, !is.na(trajectory)) |> pull(trajectory)
traj_levels    <- levels(variants_classified$trajectory)

mat <- matrix(0, nrow = 2, ncol = length(traj_levels),
              dimnames = list(c("ARHG","Genome-wide"), traj_levels))
for (tl in traj_levels) {
  mat["ARHG", tl]         <- sum(arhg_counts == tl, na.rm = TRUE)
  mat["Genome-wide", tl]  <- sum(genome_counts == tl, na.rm = TRUE)
}
message("\nContingency table (ARHG vs Genome-wide):")
print(mat)
if (all(mat >= 0) && nrow(mat) == 2 && ncol(mat) == length(traj_levels)) {
  chisq_res <- tryCatch(chisq.test(mat), error = function(e) NULL)
  if (!is.null(chisq_res)) {
    message(sprintf("Chi-square: X2=%.3f, df=%d, p=%.4g", chisq_res$statistic,
                    chisq_res$parameter, chisq_res$p.value))
  }
}

p_bar <- ggplot(prop_df |> filter(!is.na(group)),
                aes(x = group, y = pct, fill = trajectory)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = COLORS, name = "Trajectory") +
  scale_y_continuous(labels = scales::label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  geom_text(aes(label = total, y = 102), size = 2.8, vjust = 0, color = "grey40") +
  labs(
    title = "Clonal trajectory proportions by gene group",
    x     = NULL,
    y     = "% of variants",
    caption = "Numbers above bars = total variant count"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/clonal_trajectory_bar_vs_genome.pdf"),
       plot = p_bar, width = 6.5, height = 4.5, units = "in")
message("Saved: Cohorts/Relapse/Figures/clonal_trajectory_bar_vs_genome.pdf")

# ── 7. VAF scatter (ARHG variants, 4-class color) ─────────────────────────────
# Use all-variant background in grey, ARHG on top in color
bg_df <- variants_classified |>
  filter(!is_arhg) |>
  mutate(vaf_dx_plot  = coalesce(vaf_dx,  0),
         vaf_rel_plot = coalesce(vaf_rel, 0))

fg_df <- variants_classified |>
  filter(is_arhg) |>
  mutate(vaf_dx_plot  = coalesce(vaf_dx,  0),
         vaf_rel_plot = coalesce(vaf_rel, 0))

p_scatter <- ggplot() +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_point(data = bg_df,
             aes(x = vaf_dx_plot, y = vaf_rel_plot),
             color = "grey80", size = 1.2, alpha = 0.4) +
  geom_point(data = fg_df,
             aes(x = vaf_dx_plot, y = vaf_rel_plot, color = trajectory),
             size = 3, alpha = 0.85) +
  geom_text_repel(
    data          = fg_df,
    aes(x = vaf_dx_plot, y = vaf_rel_plot, label = Gene, color = trajectory),
    size          = 3,
    fontface      = "italic",
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.3,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +
  scale_color_manual(values = COLORS, name = "Trajectory") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent) +
  annotate("text", x = 0.02, y = 0.98,
           label = sprintf("ARHG: n=%d variants\n%d patients",
                           nrow(fg_df), n_distinct(fg_df$Pair_ID)),
           hjust = 0, vjust = 1, size = 3, color = "grey30") +
  labs(
    title = "Clonal trajectory: Diagnosis -> Relapse",
    x     = "Diagnosis VAF",
    y     = "Relapse VAF",
    subtitle = "ARHG variants highlighted; grey = all other variants"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 9),
    legend.position = "right"
  )

ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/clonal_trajectory_scatter.pdf"),
       plot = p_scatter, width = 5.5, height = 5.5, units = "in")
message("Saved: Cohorts/Relapse/Figures/clonal_trajectory_scatter.pdf")

# ── 8. Swimmer plot (ARHG patients) ───────────────────────────────────────────
swim_df <- arhg_variants |>
  filter(!is.na(trajectory)) |>
  mutate(
    vaf_dx_plot  = coalesce(vaf_dx,  0),
    vaf_rel_plot = coalesce(vaf_rel, 0),
    label        = Gene
  ) |>
  arrange(Pair_ID, Gene)

# One segment per variant: Dx VAF to Rel VAF
# x-axis = VAF, one row per patient, colored by trajectory
p_swim <- ggplot(swim_df) +
  geom_segment(
    aes(x = vaf_dx_plot, xend = vaf_rel_plot,
        y = Pair_ID,      yend = Pair_ID,
        color = trajectory),
    linewidth = 1.2, alpha = 0.8,
    arrow = arrow(length = unit(4, "pt"), type = "closed")
  ) +
  geom_point(aes(x = vaf_dx_plot,  y = Pair_ID, color = trajectory),
             shape = 21, fill = "white", size = 2.5, stroke = 1.2) +
  geom_point(aes(x = vaf_rel_plot, y = Pair_ID, color = trajectory),
             size = 2.5) +
  geom_text_repel(
    aes(x = pmax(vaf_dx_plot, vaf_rel_plot), y = Pair_ID, label = label),
    nudge_x       = 0.03,
    hjust         = 0,
    size          = 2.8,
    fontface      = "italic",
    color         = "grey30",
    segment.color = "grey70",
    segment.size  = 0.3,
    direction     = "y",
    max.overlaps  = Inf
  ) +
  scale_color_manual(values = COLORS, name = "Trajectory") +
  scale_x_continuous(
    limits = c(0, 1.15),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = scales::percent
  ) +
  labs(
    title    = "ARHG clonal evolution: Diagnosis -> Relapse",
    subtitle = "Arrow: Dx VAF (open) -> Relapse VAF (filled). Arrow direction = expansion/contraction.",
    x        = "Variant Allele Frequency",
    y        = "Patient"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 8.5),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    legend.position    = "right"
  )

n_swim_patients <- n_distinct(swim_df$Pair_ID)
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/Figures/clonal_trajectory_swimmer.pdf"),
       plot = p_swim,
       width  = 7,
       height = max(3, n_swim_patients * 0.7 + 1.5),
       units  = "in")
message("Saved: Cohorts/Relapse/Figures/clonal_trajectory_swimmer.pdf")

# ── 9. Clonal dominance at Dx ─────────────────────────────────────────────────
clonal_df <- variants_classified |>
  filter(present_dx, !is.na(clonal_dx))

message("\n=== Clonal dominance at Dx (VAF > 30%) ===")
message("All variants:")
print(table(clonal_df$clonal_dx))
message(sprintf("  Fraction clonal: %.1f%%",
                mean(clonal_df$clonal_dx == "Clonal") * 100))

message("\nARHG variants:")
arhg_clonal <- clonal_df |> filter(is_arhg)
if (nrow(arhg_clonal) > 0) {
  print(table(arhg_clonal$clonal_dx))
  message(sprintf("  Fraction clonal: %.1f%%",
                  mean(arhg_clonal$clonal_dx == "Clonal") * 100))

  # Fisher: ARHG clonal vs genome-wide clonal
  fish_mat <- matrix(c(
    sum(arhg_clonal$clonal_dx == "Clonal"),
    sum(arhg_clonal$clonal_dx == "Subclonal"),
    sum(clonal_df$clonal_dx[!clonal_df$is_arhg] == "Clonal"),
    sum(clonal_df$clonal_dx[!clonal_df$is_arhg] == "Subclonal")
  ), nrow = 2,
  dimnames = list(c("Clonal","Subclonal"), c("ARHG","Genome-wide")))
  message("Fisher's exact (ARHG vs genome-wide clonal status at Dx):")
  print(fish_mat)
  ft <- fisher.test(fish_mat)
  message(sprintf("  OR=%.2f (95%% CI %.2f–%.2f), p=%.4g",
                  ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
} else {
  message("  No ARHG variants present at Dx in this dataset.")
}

# Also compare clonal fraction at Dx vs Rel for ARHG
message("\nARHG clonal at Relapse (VAF > 30%):")
arhg_rel <- variants_classified |> filter(is_arhg, present_rel)
if (nrow(arhg_rel) > 0) {
  rel_clonal <- sum(arhg_rel$vaf_rel > 0.30)
  message(sprintf("  %d/%d (%.1f%%) ARHG variants clonal at Relapse",
                  rel_clonal, nrow(arhg_rel), rel_clonal/nrow(arhg_rel)*100))
}

# ── 10. Flag priority wet lab targets ─────────────────────────────────────────
# ARHG variants with Expanded or Gained trajectory in ≥2 independent patients
message("\n=== Priority wet lab targets (Expanded or Gained in ≥2 patients) ===")
priority_targets <- arhg_variants |>
  filter(trajectory %in% c("Expanded", "Gained")) |>
  group_by(Gene, `MANE HGVS`) |>
  summarise(
    n_patients  = n_distinct(Pair_ID),
    patients    = paste(sort(unique(Pair_ID)), collapse = ", "),
    trajectory  = paste(sort(unique(as.character(trajectory))), collapse = "/"),
    vaf_dx_range  = sprintf("%.2f–%.2f", min(coalesce(vaf_dx, 0)), max(coalesce(vaf_dx, 0))),
    vaf_rel_range = sprintf("%.2f–%.2f", min(coalesce(vaf_rel, 0)), max(coalesce(vaf_rel, 0))),
    .groups = "drop"
  ) |>
  filter(n_patients >= 2) |>
  arrange(desc(n_patients), Gene)

if (nrow(priority_targets) > 0) {
  message(sprintf("  %d variants meet criteria:", nrow(priority_targets)))
  print(priority_targets)
} else {
  message("  No ARHG variants in ≥2 patients with Expanded/Gained trajectory.")
  message("  All Expanded/Gained ARHG variants (any n):")
  print(
    arhg_variants |>
      filter(trajectory %in% c("Expanded", "Gained")) |>
      select(Pair_ID, Gene, `MANE HGVS`, vaf_dx, vaf_rel, delta, trajectory)
  )
}

message("\n=== Done ===")
message("Outputs:")
message("  Cohorts/Relapse/Output/arhg_trajectory_classified.csv")
message("  Cohorts/Relapse/Figures/clonal_trajectory_scatter.pdf")
message("  Cohorts/Relapse/Figures/clonal_trajectory_swimmer.pdf")
message("  Cohorts/Relapse/Figures/clonal_trajectory_bar_vs_genome.pdf")
