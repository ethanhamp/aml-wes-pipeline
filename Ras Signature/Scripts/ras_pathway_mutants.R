################################################################################
# ras_pathway_mutants.R
# ----------------------
# Identifies ARHG-mutant patients who also carry somatic RAS pathway mutations.
# Outputs: ras_pathway_mutants.csv — portable exclusion list for Python/R analyses.
#
# Logic:
#   - Source: arhg_mutant_unified_long.rds (all variants for ARHG-mutant cohort)
#   - Somatic filter: timepoint != "Germline" (keeps Dx, Rel, Tumor)
#   - Quality filters: AD Total > 20, VAF >= 0.02 (per CLAUDE.md project rules)
#   - Gene filter: gene %in% RAS pathway gene list
#
# Output columns:
#   patient_id      — Eisfeld patient ID
#   ras_score_id    — ID used in applied_ras_scores.csv (with _blast suffix if applicable)
#   gene            — RAS pathway gene mutated
#   hgvs_mane       — HGVS notation at MANE transcript
#   vaf             — variant allele frequency
#   ad_total        — total read depth at variant site
#   timepoint       — Dx / Rel / Tumor
#   effect          — predicted protein consequence
#
# RAS pathway genes (same set used in BeatAML LASSO training):
#   NRAS, KRAS, HRAS, PTPN11, NF1, CBL, RRAS, RRAS2, RAF1, BRAF, MAP2K1, MAP2K2
################################################################################

library(dplyr)
library(readr)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
BASE_DIR  <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
RDS_FILE  <- file.path(BASE_DIR, "Cohorts/ARHG/arhg_mutant_unified_long.rds")
SCORE_CSV <- file.path(BASE_DIR, "Ras Signature/applied_ras_scores.csv")
OUT_CSV   <- file.path(BASE_DIR, "Ras Signature/ras_pathway_mutants.csv")

# ---------------------------------------------------------------------------
# RAS PATHWAY GENE LIST
# Must match the gene set used in BeatAML model training.
# ---------------------------------------------------------------------------
RAS_GENES <- c(
  "NRAS", "KRAS", "HRAS",       # canonical RAS oncogenes
  "PTPN11",                      # SHP2 — upstream activator (Noonan/JMML)
  "NF1",                         # RAS GAP — loss-of-function activates RAS
  "CBL",                         # E3 ligase / RAS regulator
  "RRAS", "RRAS2",               # R-RAS family members
  "RAF1", "BRAF",                # RAF kinases (downstream of RAS)
  "MAP2K1", "MAP2K2"             # MEK1/2 (downstream of RAF)
)

cat("============================================================\n")
cat("ras_pathway_mutants.R\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 1. LOAD RDS
# ---------------------------------------------------------------------------
cat("Loading ARHG mutant variant data ...\n")
variants <- readRDS(RDS_FILE)
cat(sprintf("  Total rows: %d  |  Patients: %d\n\n",
            nrow(variants), n_distinct(variants$patient_id)))

# ---------------------------------------------------------------------------
# 2. LOAD APPLIED RAS SCORES (to know which patients have RNA-seq scores)
#    This lets us flag whether a RAS-mutant patient actually appears in the
#    scoring CSV and thus needs to be excluded from compare_cohorts analysis.
# ---------------------------------------------------------------------------
scores <- read_csv(SCORE_CSV, show_col_types = FALSE)
cat(sprintf("Applied RAS scores: %d samples\n\n", nrow(scores)))

# ---------------------------------------------------------------------------
# 3. FILTER TO SOMATIC RAS PATHWAY MUTATIONS
#    Somatic = any timepoint that is NOT Germline.
#    Quality thresholds from CLAUDE.md: AD > 20, VAF >= 0.02.
# ---------------------------------------------------------------------------
ras_somatic <- variants |>
  filter(
    timepoint != "Germline",   # exclude constitutional/germline calls
    ad_total  >  20,            # minimum depth (low-confidence below this)
    vaf       >= 0.02,          # minimum VAF (2% allele fraction floor)
    gene      %in% RAS_GENES    # restrict to RAS pathway genes
  ) |>
  select(patient_id, sample_id, timepoint, gene, hgvs_mane,
         vaf, ad_total, effect) |>
  arrange(patient_id, gene)

cat(sprintf("Somatic RAS pathway variants passing filters: %d\n", nrow(ras_somatic)))
cat(sprintf("Unique patients with somatic RAS mutations:   %d\n\n", n_distinct(ras_somatic$patient_id)))

# ---------------------------------------------------------------------------
# 4. MAP patient_id -> ras_score_id
#    The RNA-seq scoring pipeline uses a different ID scheme:
#      - Numeric IDs (e.g., 4472) become "4472_blast" in expression matrix
#      - C-XX-XXXX IDs appear as-is or with suffixes like "-blasts"
#    Strategy: for each patient_id, find the matching row in applied_ras_scores
#    using a substring search (patient_id is contained within score sample_id).
# ---------------------------------------------------------------------------
score_ids <- scores$sample_id

map_to_score_id <- function(pid) {
  # Exact match first
  exact <- score_ids[score_ids == pid]
  if (length(exact) > 0) return(exact[1])
  # Partial: score ID contains the patient ID (e.g., "4472_blast" contains "4472")
  partial <- score_ids[grepl(pid, score_ids, fixed = TRUE)]
  if (length(partial) > 0) return(partial[1])
  # No match: patient had no RNA-seq data (not in scoring dataset)
  return(NA_character_)
}

ras_somatic <- ras_somatic |>
  mutate(ras_score_id = sapply(patient_id, map_to_score_id))

# ---------------------------------------------------------------------------
# 5. SUMMARISE — one row per patient showing all their RAS mutations
# ---------------------------------------------------------------------------
cat("Per-patient RAS mutation summary:\n")
cat("------------------------------------------------------------\n")

patient_summary <- ras_somatic |>
  group_by(patient_id, ras_score_id) |>
  summarise(
    n_ras_variants = n(),
    genes_mutated  = paste(unique(gene), collapse = ", "),
    variants_detail = paste(
      sprintf("%s %s (VAF=%.3f)", gene, hgvs_mane, vaf),
      collapse = " | "
    ),
    .groups = "drop"
  ) |>
  mutate(
    in_score_csv = !is.na(ras_score_id),
    note = case_when(
      is.na(ras_score_id) ~ "No RNA-seq score (absent from applied_ras_scores.csv)",
      TRUE                ~ "Present in applied_ras_scores.csv — EXCLUDE from RAS-clean analysis"
    )
  ) |>
  arrange(patient_id)

for (i in seq_len(nrow(patient_summary))) {
  r <- patient_summary[i, ]
  cat(sprintf("  %s\n", r$patient_id))
  cat(sprintf("    ras_score_id : %s\n", ifelse(is.na(r$ras_score_id), "<absent>", r$ras_score_id)))
  cat(sprintf("    genes        : %s\n", r$genes_mutated))
  cat(sprintf("    n variants   : %d\n", r$n_ras_variants))
  cat(sprintf("    note         : %s\n\n", r$note))
}

# ---------------------------------------------------------------------------
# 6. SAVE — full variant-level table (one row per variant)
# ---------------------------------------------------------------------------
out_df <- ras_somatic |>
  left_join(
    patient_summary |> select(patient_id, ras_score_id, in_score_csv, note),
    by = c("patient_id", "ras_score_id")
  ) |>
  select(patient_id, ras_score_id, in_score_csv, timepoint, gene,
         hgvs_mane, vaf, ad_total, effect, note)

write_csv(out_df, OUT_CSV)

cat("============================================================\n")
cat(sprintf("Output saved: %s\n", OUT_CSV))
cat(sprintf("  Rows (one per variant)        : %d\n", nrow(out_df)))
cat(sprintf("  Patients with RAS mutations   : %d\n", n_distinct(out_df$patient_id)))
cat(sprintf("  Patients IN score CSV         : %d  <- exclude these from RAS-clean analysis\n",
            sum(patient_summary$in_score_csv)))
cat(sprintf("  Patients NOT in score CSV     : %d  <- no RNA-seq data; already absent from scores\n",
            sum(!patient_summary$in_score_csv)))
cat("============================================================\n")
