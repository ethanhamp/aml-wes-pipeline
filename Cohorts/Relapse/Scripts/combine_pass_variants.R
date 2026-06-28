library(readxl)
library(dplyr)
library(readr)
library(openxlsx)   # kept for compatibility; output uses writexl (type-safe)
library(writexl)    # replaces openxlsx::writeData for output — handles numeric/logical correctly

BASE_DIR    <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
relapse_dir <- file.path(BASE_DIR, "Cohorts/Relapse")

local_env <- new.env()
source(file.path(BASE_DIR, "General/artifact_genes.txt"), local = local_env)
artifact_genes <- local_env$artifact_genes

all_files <- list.files(relapse_dir, pattern = "\\.xlsx$", full.names = TRUE)

# Split files into old-style and new-style (GSL-AE-2110 batch + C-05-1267)
is_new_style <- function(f) grepl("^210428_|^C-05-1267_", basename(f))

old_germline_files <- all_files[grepl("_germline",       basename(all_files)) & !is_new_style(all_files)]
old_somatic_files  <- all_files[grepl("_mutect2_somatic", basename(all_files)) & !is_new_style(all_files)]
new_germline_files <- all_files[grepl("_germline",       basename(all_files)) &  is_new_style(all_files)]
new_somatic_files  <- all_files[grepl("_mutect2_somatic", basename(all_files)) &  is_new_style(all_files)]

# ── Pairing table for new-style files ──────────────────────────────────────────
# germline_file_id → pair_id (= diagnosis sample ID) + dx_base for matching T-N VAF cols
# Copy config/relapse_pairs.example.csv → config/relapse_pairs.csv and fill in your IDs.
PAIRS_FILE <- file.path(BASE_DIR, "config/relapse_pairs.csv")
pairing_2110 <- if (file.exists(PAIRS_FILE)) {
  pairs_df <- read.csv(PAIRS_FILE, comment.char = "#", stringsAsFactors = FALSE)
  setNames(
    lapply(seq_len(nrow(pairs_df)), function(i)
      list(pair_id = pairs_df$pair_id[i], dx_base = pairs_df$dx_base[i])),
    pairs_df$relapse_id
  )
} else {
  list()
}

# ── Shared helpers ─────────────────────────────────────────────────────────────
strip_readxl_suffixes <- function(col_names) sub("\\.\\.\\.\\d+$", "", col_names)

drop_duplicate_cols <- function(df) df[, !duplicated(names(df)), drop = FALSE]

filter_min_depth <- function(df, min_depth = 20) {
  ad_cols <- grep("^(BT|Dx|Relapse|Blasts) AD Total$", names(df), value = TRUE)
  if (length(ad_cols) == 0) return(df)
  keep <- rowSums(sapply(ad_cols, function(col) {
    !is.na(df[[col]]) & as.numeric(df[[col]]) < min_depth
  })) == 0
  df[keep, , drop = FALSE]
}

# Strip the primary-sample prefix to produce generic column names
rename_primary_sample_cols <- function(col_names, prefix) {
  col_names <- sub(paste0("^", prefix, " "), "", col_names)
  col_names <- sub(paste0("^", prefix, "_"), "", col_names)
  col_names <- sub("^Alt_Percent$",          "Alt Percentage", col_names)
  col_names <- sub("^T-N VAF \\(Relapse\\)$", "T-N VAF",       col_names)
  col_names <- sub("^T-N VAF \\(Dx\\)$",      "T-N VAF",       col_names)
  col_names
}

# ── Old-style helpers ──────────────────────────────────────────────────────────
extract_sample_id <- function(filepath) {
  sub("_(germline|mutect2_somatic).*", "", basename(filepath))
}

generalize_germline_cols <- function(col_names, sample_id) {
  col_names <- strip_readxl_suffixes(col_names)
  col_names <- gsub(paste0(sample_id, "-BT"), "BT", col_names, fixed = TRUE)
  col_names
}

generalize_somatic_cols <- function(col_names, sample_id) {
  col_names <- strip_readxl_suffixes(col_names)
  col_names <- gsub(paste0(sample_id, "-BT"),     "BT",     col_names, fixed = TRUE)
  col_names <- gsub(paste0(sample_id, "-blasts"), "Blasts", col_names, fixed = TRUE)
  col_names <- gsub("[CU]-\\d{2}-\\d{3,4}(-[A-Za-z0-9]+)?", "Relapse", col_names, perl = TRUE)
  col_names
}

extract_relapse_sample_id <- function(col_names, sample_id) {
  stripped <- strip_readxl_suffixes(col_names)
  matches  <- regmatches(stripped, gregexpr("[CU]-\\d{2}-\\d{3,4}", stripped, perl = TRUE))
  ids <- unique(unlist(matches))
  ids <- ids[!grepl(sample_id, ids, fixed = TRUE)]
  if (length(ids) == 0) return(NA_character_)
  ids[1]
}

# ── New-style helpers ──────────────────────────────────────────────────────────
extract_new_germline_id <- function(filepath) {
  bn <- basename(filepath)
  bn <- sub("^210428_Eisfeld_GSL-AE-2110_", "", bn)
  bn <- sub("_(germline|mutect2_somatic).*", "", bn)
  bn
}

# Rename all sample-specific cols using exact string replacement
rename_new_somatic_cols <- function(col_names, germline_col, dx_col, rel_col) {
  n <- strip_readxl_suffixes(col_names)
  for (pair in list(list(germline_col, "BT"), list(dx_col, "Dx"), list(rel_col, "Relapse"))) {
    old <- pair[[1]]; new <- pair[[2]]
    n <- gsub(paste0(old, " "), paste0(new, " "), n, fixed = TRUE)
    n <- gsub(paste0(old, "_"), paste0(new, "_"), n, fixed = TRUE)
  }
  n <- gsub(paste0("T-N VAF (", dx_col,  ")"), "T-N VAF (Dx)",      n, fixed = TRUE)
  n <- gsub(paste0("T-N VAF (", rel_col, ")"), "T-N VAF (Relapse)", n, fixed = TRUE)
  n
}

rename_new_germline_cols <- function(col_names, germline_col) {
  n <- strip_readxl_suffixes(col_names)
  n <- gsub(paste0(germline_col, " "), "BT ", n, fixed = TRUE)
  n <- gsub(paste0(germline_col, "_"), "BT_", n, fixed = TRUE)
  n
}

# Detect the germline column prefix from an Alt Percentage/Alt_Percent column
detect_germline_col <- function(col_names) {
  stripped <- strip_readxl_suffixes(col_names)
  alt_col  <- stripped[grepl("Alt[_ ]?Percent", stripped)][1]
  if (is.na(alt_col)) return(NULL)
  sub("[ _]Alt[_ ]?Percent.*$", "", alt_col)
}

# ── Read functions ─────────────────────────────────────────────────────────────
read_germline <- function(filepath) {
  df <- tryCatch(read_excel(filepath), error = function(e) { message("ERROR: ", filepath); NULL })
  if (is.null(df)) return(NULL)
  sample_id <- extract_sample_id(filepath)
  names(df) <- generalize_germline_cols(names(df), sample_id)
  df <- drop_duplicate_cols(df)
  n_before <- nrow(df)
  df <- filter_min_depth(df)
  names(df) <- rename_primary_sample_cols(names(df), "BT")
  df <- bind_cols(tibble(Sample_ID = sample_id, Timepoint = "germline", Pair_ID = sample_id), df)
  message(sprintf("  %-50s  %d rows (%d removed <20 reads)", basename(filepath), nrow(df), n_before - nrow(df)))
  df
}

read_somatic_pass <- function(filepath, filter_col = "Filter") {
  df <- tryCatch(read_excel(filepath), error = function(e) { message("ERROR: ", filepath); NULL })
  if (is.null(df)) return(NULL)
  if (!filter_col %in% names(df)) { message("No Filter col: ", basename(filepath)); return(NULL) }

  sample_id  <- extract_sample_id(filepath)
  relapse_id <- extract_relapse_sample_id(names(df), sample_id)

  df <- df[!is.na(df[[filter_col]]) & df[[filter_col]] == "PASS", ]
  if (nrow(df) == 0) { message("No PASS rows: ", basename(filepath)); return(NULL) }

  names(df) <- generalize_somatic_cols(names(df), sample_id)
  df <- drop_duplicate_cols(df)
  n_before <- nrow(df)
  df <- filter_min_depth(df)

  dx_rows       <- df; names(dx_rows)  <- rename_primary_sample_cols(names(df), "Blasts")
  # rename_primary_sample_cols always renames T-N VAF (Relapse) → T-N VAF, so the Dx
  # row's T-N VAF is currently the relapse value. Swap it with T-N VAF (Blasts).
  if ("T-N VAF" %in% names(dx_rows) && "T-N VAF (Blasts)" %in% names(dx_rows)) {
    names(dx_rows)[names(dx_rows) == "T-N VAF"]          <- "T-N VAF (Relapse)"
    names(dx_rows)[names(dx_rows) == "T-N VAF (Blasts)"] <- "T-N VAF"
  }
  dx_rows       <- bind_cols(tibble(Sample_ID = sample_id,  Timepoint = "diagnosis", Pair_ID = sample_id), dx_rows)
  rel_rows      <- df; names(rel_rows) <- rename_primary_sample_cols(names(df), "Relapse")
  rel_rows      <- bind_cols(tibble(Sample_ID = relapse_id, Timepoint = "relapse",   Pair_ID = sample_id), rel_rows)

  out <- bind_rows(dx_rows, rel_rows)
  message(sprintf("  %-50s  %d PASS rows (%d removed <20 reads) → %d rows (dx+rel)",
                  basename(filepath), nrow(df), n_before - nrow(df), nrow(out)))
  out
}

read_new_germline <- function(filepath) {
  df <- tryCatch(read_excel(filepath), error = function(e) { message("ERROR: ", filepath); NULL })
  if (is.null(df)) return(NULL)

  germline_file_id <- extract_new_germline_id(filepath)
  pair_info        <- pairing_2110[[germline_file_id]]
  if (is.null(pair_info)) { message("No pairing for: ", basename(filepath)); return(NULL) }

  germline_col <- detect_germline_col(names(df))
  if (is.null(germline_col)) { message("Cannot detect germline col in: ", basename(filepath)); return(NULL) }

  names(df) <- rename_new_germline_cols(names(df), germline_col)
  df <- drop_duplicate_cols(df)
  n_before <- nrow(df)
  df <- filter_min_depth(df)
  names(df) <- rename_primary_sample_cols(names(df), "BT")
  df <- bind_cols(tibble(Sample_ID = germline_file_id, Timepoint = "germline", Pair_ID = pair_info$pair_id), df)
  message(sprintf("  %-50s  %d rows (%d removed <20 reads)", basename(filepath), nrow(df), n_before - nrow(df)))
  df
}

read_new_somatic_pass <- function(filepath, filter_col = "Filter") {
  df <- tryCatch(read_excel(filepath), error = function(e) { message("ERROR: ", filepath); NULL })
  if (is.null(df)) return(NULL)
  if (!filter_col %in% names(df)) { message("No Filter col: ", basename(filepath)); return(NULL) }

  germline_file_id <- extract_new_germline_id(filepath)
  pair_info        <- pairing_2110[[germline_file_id]]
  if (is.null(pair_info)) { message("No pairing for: ", basename(filepath)); return(NULL) }

  # Detect dx and rel from T-N VAF columns — most reliable approach
  tn_ids <- sub("^T-N VAF \\((.+)\\)$", "\\1",
                names(df)[grepl("^T-N VAF \\(", names(df))])
  tn_ids <- unique(strip_readxl_suffixes(tn_ids))

  if (length(tn_ids) < 2) { message("Cannot detect dx/rel in: ", basename(filepath)); return(NULL) }

  dx_col  <- tn_ids[ grepl(pair_info$dx_base, tn_ids, fixed = TRUE)]
  rel_col <- tn_ids[!grepl(pair_info$dx_base, tn_ids, fixed = TRUE)]

  if (length(dx_col) != 1 || length(rel_col) != 1) {
    message("Ambiguous dx/rel for: ", basename(filepath)); return(NULL)
  }

  # Detect germline col: the Alt Percentage col that is neither dx nor rel
  stripped_names <- strip_readxl_suffixes(names(df))
  alt_cols <- stripped_names[grepl("Alt Percentage$", stripped_names)]
  germline_col_full <- alt_cols[!grepl(paste(c(dx_col, rel_col), collapse = "|"), alt_cols)]
  germline_col <- sub(" Alt Percentage$", "", germline_col_full[1])

  df <- df[!is.na(df[[filter_col]]) & df[[filter_col]] == "PASS", ]
  if (nrow(df) == 0) { message("No PASS rows: ", basename(filepath)); return(NULL) }

  names(df) <- rename_new_somatic_cols(names(df), germline_col, dx_col, rel_col)
  df <- drop_duplicate_cols(df)
  n_before <- nrow(df)
  df <- filter_min_depth(df)

  # New-style files have both "T-N VAF (Dx)" and "T-N VAF (Relapse)" in the same data
  # frame.  rename_primary_sample_cols converts *both* patterns to "T-N VAF", creating
  # a duplicate column pair.  Fix: pre-protect the unwanted column with a placeholder
  # name before calling rename_primary_sample_cols, then drop it afterwards.
  protect_tn_vaf <- function(col_names, protect_pattern, placeholder = "T-N VAF __DROP__") {
    col_names[col_names == protect_pattern] <- placeholder
    col_names
  }

  dx_rows  <- df
  # Protect T-N VAF (Relapse) so rename_primary_sample_cols won't convert it to T-N VAF;
  # only T-N VAF (Dx) → T-N VAF will happen.
  names(dx_rows) <- protect_tn_vaf(names(dx_rows), "T-N VAF (Relapse)")
  names(dx_rows) <- rename_primary_sample_cols(names(dx_rows), "Dx")
  dx_rows <- dx_rows[, names(dx_rows) != "T-N VAF __DROP__", drop = FALSE]
  dx_rows  <- bind_cols(tibble(Sample_ID = pair_info$pair_id, Timepoint = "diagnosis", Pair_ID = pair_info$pair_id), dx_rows)

  rel_rows <- df
  # Protect T-N VAF (Dx) so rename_primary_sample_cols won't convert it to T-N VAF;
  # only T-N VAF (Relapse) → T-N VAF will happen.
  names(rel_rows) <- protect_tn_vaf(names(rel_rows), "T-N VAF (Dx)")
  names(rel_rows) <- rename_primary_sample_cols(names(rel_rows), "Relapse")
  rel_rows <- rel_rows[, names(rel_rows) != "T-N VAF __DROP__", drop = FALSE]
  rel_rows <- bind_cols(tibble(Sample_ID = rel_col, Timepoint = "relapse", Pair_ID = pair_info$pair_id), rel_rows)

  out <- bind_rows(dx_rows, rel_rows)
  message(sprintf("  %-50s  %d PASS rows (%d removed <20 reads) → %d rows (dx+rel)",
                  basename(filepath), nrow(df), n_before - nrow(df), nrow(out)))
  out
}

# ── Run all ────────────────────────────────────────────────────────────────────
message("\n--- Old-style germline files ---")
old_germline_list <- lapply(old_germline_files, read_germline)

message("\n--- Old-style somatic files ---")
old_somatic_list <- lapply(old_somatic_files, read_somatic_pass)

message("\n--- New-style germline files (GSL-AE-2110) ---")
new_germline_list <- lapply(new_germline_files, read_new_germline)

message("\n--- New-style somatic files (GSL-AE-2110) ---")
new_somatic_list <- lapply(new_somatic_files, read_new_somatic_pass)

# ── Filters ────────────────────────────────────────────────────────────────────
EXCLUDED_FILE   <- file.path(BASE_DIR, "config/excluded_samples.txt")
exclude_samples <- if (file.exists(EXCLUDED_FILE)) {
  lines <- readLines(EXCLUDED_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else {
  character(0)
}

# ── Combine ────────────────────────────────────────────────────────────────────
# Somatic lists are placed FIRST so that the T-N VAF columns (present only in
# somatic rows) appear early in the combined data frame.  Both writexl and
# openxlsx infer the Excel column type from the first non-empty cell; if all
# leading rows have NA for T-N VAF (as happens when germline rows come first),
# the column is mistyped as logical/boolean in the output file.
all_lists <- c(old_somatic_list, new_somatic_list, old_germline_list, new_germline_list)
all_variants <- bind_rows(lapply(Filter(Negate(is.null), all_lists),
                                 function(df) mutate(df, across(everything(), as.character))))
all_variants <- suppressMessages(type_convert(all_variants))
all_variants <- all_variants[!all_variants$Sample_ID %in% exclude_samples &
                             !all_variants$Pair_ID   %in% exclude_samples, ]
all_variants <- all_variants[!is.na(all_variants$Gene) & !all_variants$Gene %in% artifact_genes, ]

# ── Patient_Group column ───────────────────────────────────────────────────────
# Rank unique Pair_IDs alphabetically and assign zero-padded sequential labels
# (Patient_01, Patient_02, ...) so all three timepoints for the same patient
# share a single clean label.  Pair_ID already correctly links germline, dx, and
# relapse rows, so deriving Patient_Group from it guarantees consistency.
pair_id_levels <- sort(unique(all_variants$Pair_ID))
n_patients     <- length(pair_id_levels)
pad_width      <- nchar(as.character(n_patients))   # e.g. 2 → "01", 3 → "001"
patient_group_map <- setNames(
  sprintf(paste0("Patient_%0", pad_width, "d"), seq_along(pair_id_levels)),
  pair_id_levels
)
all_variants <- all_variants %>%
  mutate(Patient_Group = patient_group_map[Pair_ID]) %>%
  relocate(Patient_Group, .after = Sample_ID)

message(sprintf("\nTotal: %d rows, %d columns", nrow(all_variants), ncol(all_variants)))
message(sprintf("  Germline:  %d rows, %d samples",
                sum(all_variants$Timepoint == "germline"),
                n_distinct(all_variants$Sample_ID[all_variants$Timepoint == "germline"])))
message(sprintf("  Diagnosis: %d rows, %d samples",
                sum(all_variants$Timepoint == "diagnosis"),
                n_distinct(all_variants$Sample_ID[all_variants$Timepoint == "diagnosis"])))
message(sprintf("  Relapse:   %d rows, %d samples",
                sum(all_variants$Timepoint == "relapse"),
                n_distinct(all_variants$Sample_ID[all_variants$Timepoint == "relapse"])))

out_file <- file.path(relapse_dir, "relapse_pass_variants_combined.xlsx")
# writexl preserves column types faithfully (numeric stays numeric).
# openxlsx ≤4.2.8.x misidentifies [0,1]-range numeric columns as logical (boolean).
write_xlsx(list(All_Variants = all_variants), path = out_file)
message("\nSaved: ", out_file)
