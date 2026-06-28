library(dplyr)
library(readxl)
library(writexl)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR    <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
VARIANTS_F  <- file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx")
ARHG_IDS_F  <- file.path(BASE_DIR, "arhg_samples.xlsx")
OUT_XLSX    <- file.path(BASE_DIR, "Exomes/Datasets/arhg_mutant_combined_exomes.xlsx")
OUT_RDS     <- file.path(BASE_DIR, "Exomes/Datasets/arhg_mutant_combined_exomes.rds")
ARTIFACT_GENE_FILE <- file.path(BASE_DIR, "General/artifact_genes.txt")

# ── Artifact gene exclusion list ──────────────────────────────────────────────
artifact_genes <- if (file.exists(ARTIFACT_GENE_FILE)) {
  lines <- readLines(ARTIFACT_GENE_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else { character(0) }

# ── ARHG gene universe (GEFs + GAPs) ─────────────────────────────────────────
arhg_genes <- c(
  "ARHGEF1", "ARHGEF2", "ARHGEF3", "ARHGEF4", "ARHGEF5", "ARHGEF6",
  "ARHGEF7", "ARHGEF8", "ARHGEF9", "ARHGEF10", "ARHGEF11", "ARHGEF12",
  "ARHGEF15", "ARHGEF16", "ARHGEF17", "ARHGEF18", "ARHGEF19", "ARHGEF25",
  "ARHGEF26", "ARHGEF28",
  "ARHGAP1",  "ARHGAP2",  "ARHGAP3",  "ARHGAP4",  "ARHGAP5",  "ARHGAP6",
  "ARHGAP7",  "ARHGAP8",  "ARHGAP9",  "ARHGAP10", "ARHGAP11A","ARHGAP11B",
  "ARHGAP12", "ARHGAP13", "ARHGAP14", "ARHGAP15", "ARHGAP16", "ARHGAP17",
  "ARHGAP18", "ARHGAP19", "ARHGAP20", "ARHGAP21", "ARHGAP22", "ARHGAP23",
  "ARHGAP24", "ARHGAP25", "ARHGAP26", "ARHGAP27", "ARHGAP28", "ARHGAP29",
  "ARHGAP30", "ARHGAP31", "ARHGAP32", "ARHGAP33", "ARHGAP34", "ARHGAP35",
  "ARHGAP36", "ARHGAP37", "ARHGAP38", "ARHGAP39", "ARHGAP40", "ARHGAP41",
  "ARHGAP42", "ARHGAP43", "ARHGAP44", "ARHGAP45"
)

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading combined variant file...")
raw <- read_excel(VARIANTS_F)

# Externally curated list of ARHG-mutant patient sample IDs (Dx timepoint IDs).
# These are mostly from newer sequencing runs not yet in the combined relapse
# file — matched below via Pair_ID (which always stores the Dx sample ID).
arhg_ids <- read_excel(ARHG_IDS_F) %>% pull(sample.id)

# ── QC filters ────────────────────────────────────────────────────────────────
message("Applying QC filters...")
filtered <- raw %>%
  filter(
    `AD Total`        >  20,
    `Alt Percentage`  >= 0.02,     # stored as fraction (0.02 = 2%)
    !Gene %in% artifact_genes
  )

# Hypermutated groups: in the combined relapse file Patient_Group uses
# "Patient_XX" labels. Groups 11 and 4 from the 3801 cohort scripts correspond
# to Patient_11 and Patient_04 in this nomenclature — exclude both.
filtered <- filtered %>%
  filter(!Patient_Group %in% c("Patient_04", "Patient_11"))

# ── Identify ARHG-mutant patients ─────────────────────────────────────────────
# Pair_ID is the Dx-timepoint sample ID and is consistent across all timepoints
# for a given patient, making it the reliable patient key in this file.

# Source 1: patients with an ARHG variant called in this dataset
arhg_by_variant <- filtered %>%
  filter(Gene %in% arhg_genes) %>%
  pull(Pair_ID) %>%
  unique()

# Source 2: patients listed in the curated arhg_samples.xlsx
# (matches on Pair_ID because arhg_samples stores Dx IDs)
arhg_by_list <- intersect(arhg_ids, unique(filtered$Pair_ID))

# Union of both sources — captures patients regardless of which run they came from
arhg_patients <- union(arhg_by_variant, arhg_by_list)

message(sprintf(
  "ARHG-mutant patients identified: %d via variant calls, %d via curated list, %d total union",
  length(arhg_by_variant), length(arhg_by_list), length(arhg_patients)
))

# ── Subset to ARHG-mutant patients (all their variants, not just ARHG variants) ─
out <- filtered %>%
  filter(Pair_ID %in% arhg_patients)

# ── Summary ───────────────────────────────────────────────────────────────────
n_patients  <- length(unique(out$Pair_ID))
n_variants  <- nrow(out)
tp_counts   <- table(out$Timepoint)

message("\n── Output summary ───────────────────────────────────────────────────────────")
message(sprintf("Patients:  %d", n_patients))
message(sprintf("Variants:  %d total", n_variants))
message("Timepoint breakdown:")
for (tp in names(tp_counts)) {
  message(sprintf("  %-12s %d", tp, tp_counts[[tp]]))
}

# ── Write outputs ─────────────────────────────────────────────────────────────
message("\nWriting outputs...")
write_xlsx(out, OUT_XLSX)
saveRDS(out, OUT_RDS)
message(sprintf("XLSX: %s", OUT_XLSX))
message(sprintf("RDS:  %s", OUT_RDS))
