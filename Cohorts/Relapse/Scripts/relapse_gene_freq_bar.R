library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# ── Constants ──────────────────────────────────────────────────────────────────
VAF_THRESHOLD <- 0.02
TOP_N_GENES   <- 20

# Artifact / recurrently mutated normal-tissue genes to exclude.
# Identical vector used in gainloss_bar_relapse.R and mutation_burden_relapse.R.
local_env <- new.env()
source(file.path(BASE_DIR, "General/artifact_genes.txt"), local = local_env)
artifact_genes <- local_env$artifact_genes

# ── 1. Load and filter to relapse timepoint ────────────────────────────────────
# We restrict to Timepoint == "relapse" only; germline and diagnosis rows are
# intentionally excluded — this chart describes the relapse mutational landscape.
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))

df_relapse <- df |>
  filter(
    Timepoint == "relapse",
    !Gene %in% artifact_genes,
    !Patient_Group %in% c("Patient_23", "Patient_25")  # hypermutated patients
  ) |>
  # Collapse all ARHG-family genes (ARHGEF*, ARHGAP*) into a single "ARHG" label
  # so we can assess aggregate ARHG pathway involvement per patient
  mutate(Gene = if_else(startsWith(Gene, "ARHG"), "ARHG", Gene))

# ── 2. Deduplicate and apply VAF threshold ─────────────────────────────────────
# Dedup at Pair_ID + Gene + HGVS level before counting patients, so a patient
# with the same variant called in multiple rows only contributes once.
# T-N VAF is the tumor-minus-normal somatic VAF; 0.02 = 2% minimum allele fraction.
variants_relapse <- df_relapse |>
  group_by(Pair_ID, Gene, HGVS) |>
  summarise(vaf = max(`T-N VAF`, na.rm = TRUE), .groups = "drop") |>
  filter(vaf >= VAF_THRESHOLD)

# ── 3. Count unique patients per gene ─────────────────────────────────────────
# One patient is counted once per gene regardless of how many distinct
# variants they carry in that gene (patient-level prevalence, not variant count).
gene_patient_counts <- variants_relapse |>
  distinct(Pair_ID, Gene) |>        # one row per patient-gene pair
  count(Gene, name = "n_patients")

# Total unique patients in the filtered relapse cohort (denominator for %)
n_total_patients <- n_distinct(variants_relapse$Pair_ID)

cat(sprintf("Total relapse patients after filters: %d\n", n_total_patients))

# ── 4. Select top N genes and compute percentage ───────────────────────────────
top_genes_df <- gene_patient_counts |>
  arrange(desc(n_patients)) |>
  slice_head(n = TOP_N_GENES) |>
  mutate(
    pct_patients = n_patients / n_total_patients * 100,
    # Factor with levels ordered most-to-least frequent so bars plot left → right
    Gene = fct_reorder(Gene, n_patients, .desc = TRUE)
  )

cat("Top genes by relapse patient count:\n")
print(top_genes_df)

# ── 5. Plot ────────────────────────────────────────────────────────────────────
# Gain color (#4DAC26) is used as the single fill — consistent with the project
# palette for relapse-gained variants in gainloss_bar_relapse.R.
FILL_COLOR <- "#4DAC26"

# Percentage labels are placed just above each bar
top_genes_df <- top_genes_df |>
  mutate(pct_label = sprintf("%.0f%%", pct_patients))

p <- ggplot(top_genes_df, aes(x = Gene, y = n_patients)) +
  geom_col(fill = FILL_COLOR, width = 0.7) +
  # Percentage annotation above each bar — gives a second scale without a
  # literal dual-axis, which ggplot2 does not support cleanly for a bar chart
  geom_text(
    aes(label = pct_label),
    vjust  = -0.4,
    size   = 2.8,
    colour = "grey30"
  ) +
  # Give enough headroom for the percentage labels
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(
    title   = "Most frequently mutated genes at Relapse",
    subtitle = sprintf("n = %d relapse patients (after exclusions)", n_total_patients),
    x       = NULL,
    y       = "Number of patients",
    caption = sprintf(
      "Filters: Timepoint == relapse | T-N VAF >= %.0f%% | excl. artifact_genes | excl. hypermutators (Patient_23, Patient_25)\nARHG* genes collapsed to 'ARHG'. Percentages shown above bars = %% of %d patients.",
      VAF_THRESHOLD * 100, n_total_patients
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, colour = "grey40", size = 9),
    plot.caption    = element_text(hjust = 0, colour = "grey50", size = 7,
                                   margin = margin(t = 6)),
    axis.text.x     = element_text(face = "italic", angle = 45, hjust = 1,
                                   margin = margin(t = 2)),
    axis.ticks.length = unit(2, "pt"),
    axis.title.y    = element_text(margin = margin(r = 6)),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

print(p)

# ── 6. Save ────────────────────────────────────────────────────────────────────
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/relapse_gene_freq_bar.pdf"),
       plot = p, width = 7, height = 5, units = "in")
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/relapse_gene_freq_bar.png"),
       plot = p, width = 7, height = 5, units = "in", dpi = 300)

message("Saved: Cohorts/Relapse/relapse_gene_freq_bar.pdf/.png")
