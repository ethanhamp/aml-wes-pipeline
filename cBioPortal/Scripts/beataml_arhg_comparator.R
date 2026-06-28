library(cBioPortalData)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(httr)
library(jsonlite)

BASE_DIR  <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
STUDY_ID  <- "aml_ohsu_2022"   # BeatAML 2.0 (942 samples); swap to "aml_ohsu_2018" for wave 1 only
OUT_DIR   <- file.path(BASE_DIR, "cBioPortal/Datasets")

arhg_genes <- c(
  "ARHGEF1",  "ARHGEF2",  "ARHGEF3",  "ARHGEF4",  "ARHGEF5",
  "ARHGEF6",  "ARHGEF7",  "ARHGEF9",  "ARHGEF10", "ARHGEF11",
  "ARHGEF12", "ARHGEF15", "ARHGEF16", "ARHGEF17", "ARHGEF18",
  "ARHGEF19", "ARHGEF25", "ARHGEF26", "ARHGEF28",
  "ARHGAP1",  "ARHGAP2",  "ARHGAP3",  "ARHGAP4",  "ARHGAP5",  "ARHGAP6",
  "ARHGAP7",  "ARHGAP8",  "ARHGAP9",  "ARHGAP10", "ARHGAP11A","ARHGAP11B",
  "ARHGAP12", "ARHGAP13", "ARHGAP14", "ARHGAP15", "ARHGAP16", "ARHGAP17",
  "ARHGAP18", "ARHGAP19", "ARHGAP20", "ARHGAP21", "ARHGAP22", "ARHGAP23",
  "ARHGAP24", "ARHGAP25", "ARHGAP26", "ARHGAP27", "ARHGAP28", "ARHGAP29",
  "ARHGAP30", "ARHGAP31", "ARHGAP32", "ARHGAP33", "ARHGAP34", "ARHGAP35",
  "ARHGAP36", "ARHGAP37", "ARHGAP38", "ARHGAP39", "ARHGAP40", "ARHGAP41",
  "ARHGAP42", "ARHGAP43", "ARHGAP44", "ARHGAP45"
)

# ---------------------------------------------------------------------------
# 1. Connect and resolve Entrez IDs (needed for API call)
# ---------------------------------------------------------------------------
cbio <- cBioPortal()

entrez_raw <- mapIds(org.Hs.eg.db,
                     keys      = arhg_genes,
                     column    = "ENTREZID",
                     keytype   = "SYMBOL",
                     multiVals = "first")
gene_info <- data.frame(
  hugoGeneSymbol = names(entrez_raw)[!is.na(entrez_raw)],
  entrezGeneId   = as.integer(entrez_raw[!is.na(entrez_raw)]),
  stringsAsFactors = FALSE
)
entrez_ids <- gene_info$entrezGeneId
cat("Resolved", nrow(gene_info), "/", length(arhg_genes), "ARHG genes to Entrez IDs\n")

# ---------------------------------------------------------------------------
# 2. Pull clinical data and identify white samples
# ---------------------------------------------------------------------------
cat("\nPulling clinical data for", STUDY_ID, "...\n")
clin <- clinicalData(cbio, studyId = STUDY_ID)

cat("Clinical columns:", paste(names(clin), collapse = ", "), "\n")
cat("Rows:", nrow(clin), "\n\n")

# Find race column (cBioPortal uses uppercase attribute IDs)
race_col <- intersect(c("RACE", "race", "RACE_CATEGORY", "Race"), names(clin))[1]
if (is.na(race_col)) {
  stop("No race column found. Check clinical columns printed above and update race_col manually.")
}
cat("Race column identified:", race_col, "\n")
cat("Race value counts:\n")
print(table(clin[[race_col]], useNA = "ifany"))

# White samples — exclude Black/African American to keep comparator clean
white_samples <- clin |>
  filter(tolower(.data[[race_col]]) %in% c("white", "white or caucasian")) |>
  pull(sampleId)

black_samples <- clin |>
  filter(grepl("black|african", tolower(.data[[race_col]]))) |>
  pull(sampleId)

cat("\nWhite samples:", length(white_samples), "\n")
cat("Black/African American samples (excluded from comparator):", length(black_samples), "\n")

# ---------------------------------------------------------------------------
# 3. Fetch ARHG mutations via direct cBioPortal REST API call
# ---------------------------------------------------------------------------
cat("\nFetching ARHG mutations for", STUDY_ID, "via API...\n")

PROFILE_ID <- paste0(STUDY_ID, "_mutations")
SAMPLE_LIST <- paste0(STUDY_ID, "_sequenced")
API_URL <- paste0("https://www.cbioportal.org/api/molecular-profiles/",
                  PROFILE_ID, "/mutations/fetch?projection=SUMMARY")

resp <- POST(
  url  = API_URL,
  body = toJSON(list(entrezGeneIds = entrez_ids, sampleListId = SAMPLE_LIST),
                auto_unbox = TRUE),
  content_type("application/json")
)

if (status_code(resp) != 200) {
  stop("API returned ", status_code(resp), ": ", content(resp, "text", encoding = "UTF-8"))
}

all_arhg_muts <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE) |>
  mutate(across(where(is.list), ~ sapply(., function(x) if (length(x) == 1) as.character(x) else NA_character_))) |>
  left_join(gene_info, by = "entrezGeneId")
cat("ARHG mutations retrieved:", nrow(all_arhg_muts), "across",
    n_distinct(all_arhg_muts$sampleId), "samples\n")

# ---------------------------------------------------------------------------
# 4. Filter to white samples and summarize
# ---------------------------------------------------------------------------
white_arhg_muts <- all_arhg_muts |>
  filter(sampleId %in% white_samples)

n_white_total    <- length(unique(white_samples))
n_white_mutant   <- n_distinct(white_arhg_muts$sampleId)
pct_white_mutant <- round(100 * n_white_mutant / n_white_total, 1)

cat("\n========================================\n")
cat("BeatAML", STUDY_ID, "— White Comparator\n")
cat("========================================\n")
cat("Total white patients (sequenced):", n_white_total, "\n")
cat("ARHG-mutant white patients:      ", n_white_mutant, sprintf("(%s%%)\n", pct_white_mutant))
cat("----------------------------------------\n")

# Per-gene breakdown
gene_summary <- white_arhg_muts |>
  group_by(hugoGeneSymbol) |>
  summarise(
    n_mutations = n(),
    n_patients  = n_distinct(sampleId),
    .groups     = "drop"
  ) |>
  arrange(desc(n_patients))

cat("\nPer-gene breakdown (white patients):\n")
print(gene_summary, n = Inf)

# Mutation type breakdown
type_summary <- white_arhg_muts |>
  mutate(mutationType = as.character(mutationType)) |>
  dplyr::count(mutationType, sort = TRUE)
cat("\nMutation types:\n")
print(type_summary)

# ---------------------------------------------------------------------------
# 5. Save outputs
# ---------------------------------------------------------------------------
write.table(
  white_arhg_muts,
  file.path(OUT_DIR, paste0(STUDY_ID, "_white_arhg_mutations.tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  gene_summary,
  file.path(OUT_DIR, paste0(STUDY_ID, "_white_arhg_gene_summary.tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Comparator summary row (for easy merging with your Black cohort stats later)
comparator_summary <- data.frame(
  study       = STUDY_ID,
  race        = "white",
  n_total     = n_white_total,
  n_arhg_mutant = n_white_mutant,
  pct_arhg_mutant = pct_white_mutant
)
write.table(
  comparator_summary,
  file.path(OUT_DIR, paste0(STUDY_ID, "_white_comparator_summary.tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\nOutputs saved to", OUT_DIR, "\n")
cat("  -", paste0(STUDY_ID, "_white_arhg_mutations.tsv"), "\n")
cat("  -", paste0(STUDY_ID, "_white_arhg_gene_summary.tsv"), "\n")
cat("  -", paste0(STUDY_ID, "_white_comparator_summary.tsv"), "\n")
