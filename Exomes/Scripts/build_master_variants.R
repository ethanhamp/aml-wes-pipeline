library(tidyverse)
library(readxl)
library(writexl)

# ---- Configuration ----
# Each entry: list(dir, seq_run, out_dir, timepoint_map)
# timepoint_map: named character vector mapping tumor sample IDs → timepoint label.
#   Used for runs with non-standard naming (e.g. trio files where T-N VAF IDs don't
#   carry an explicit Dx/Rel suffix). Keys must match the ID inside "T-N VAF (ID)".
RUN_CONFIG <- list(
  list(
    dir     = "Cohorts/2955",
    seq_run = "GSL-PY-2955",
    out_dir = "Cohorts/2955"
  ),
  list(
    dir     = "Cohorts/2624",
    seq_run = "GSL-PY-2624",
    out_dir = "Cohorts/2624"
  ),
  list(
    dir     = "Cohorts/2434",
    seq_run = "GSL-PY-2434",
    out_dir = "Cohorts/2434"
  ),
  list(
    dir     = "Cohorts/2435",
    seq_run = "GSL-PY-2435",
    out_dir = "Cohorts/2435"
  ),
  list(
    dir     = "Cohorts/3399",
    seq_run = "GSL-EM-3399",
    out_dir = "Cohorts/3399"
  ),
  list(
    dir     = "Cohorts/2110",
    seq_run = "GSL-AE-2110",
    out_dir = "Cohorts/2110",
    # timepoint_map is required for runs where sample IDs do not carry an
    # explicit Dx/Rel suffix. Keys must match the ID inside "T-N VAF (ID)".
    # Replace with your own sample IDs and timepoint labels ("Dx" or "Rel").
    timepoint_map = c(
      "SAMPLE_001_A" = "Dx",  "SAMPLE_002_A" = "Rel"
      # add one entry per tumor sample in this run
    )
  ),
  list(
    dir     = "Cohorts/3801",
    seq_run = "GSL-EM-3801",
    out_dir = "Cohorts/3801"
  ),
  list(
    dir     = "Cohorts/3854",
    seq_run = "GSL-JB-3854",
    out_dir = "Cohorts/3854"
  ),
  list(
    dir     = "Cohorts/4394",
    seq_run = "GSL-EM-4394",
    out_dir = "Cohorts/4394"
  ),
  # UPenn second batch — files are labeled GSL-EM-4403 but belong to run 4401
  list(
    dir     = "Cohorts/4401",
    seq_run = "GSL-EM-4401",
    out_dir = "Cohorts/4401"
  ),
  list(
    dir     = "Cohorts/4403",
    seq_run = "GSL-EM-4403",
    out_dir = "Cohorts/4403"
  ),
  list(
    dir     = "Cohorts/2770",
    seq_run = "GSL-PY-2770",
    out_dir = "Cohorts/2770"
  ),
  list(
    dir     = "Cohorts/3363",
    seq_run = "GSL-PY-3363",
    out_dir = "Cohorts/3363"
  )
)

# Set BASE_DIR to your project root, or set the ARHG_BASE_DIR environment
# variable before launching R (e.g. in .Renviron or via Sys.setenv()).
BASE_DIR    <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
OUTPUT_DIR  <- file.path(BASE_DIR, "Exomes/Datasets/by_run")
MASTER_OUT  <- file.path(BASE_DIR, "Exomes/Datasets/master_variants.xlsx")
ARTIFACT_GENE_FILE   <- file.path(BASE_DIR, "General/artifact_genes.txt")

# Patient IDs to flag as hypermutated. Copy config/hypermutators.example.txt
# to config/hypermutators.txt, fill in your IDs (one per line), and that file
# will be read automatically. If absent, no samples are flagged.
HYPERMUTATOR_FILE <- file.path(BASE_DIR, "config/hypermutators.txt")
HYPERMUTATORS <- if (file.exists(HYPERMUTATOR_FILE)) {
  lines <- readLines(HYPERMUTATOR_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else {
  character(0)
}

# ---- Load artifact_genes exclusion list ----
local_env <- new.env()
source(ARTIFACT_GENE_FILE, local = local_env)
artifact_genes <- local_env$artifact_genes
cat(sprintf("Loaded %d artifact genes from artifact_genes.txt\n", length(artifact_genes)))

# ---- Helpers ----

# Extract canonical patient ID from a filename stem.
# Handles: C-XX-XXXX / U-XX-XXXX / PSXX-XXXX / case_XXXX / GSL-XX-XXXX_CASEID prefixes.
extract_patient_id <- function(stem) {
  m <- regmatches(stem, regexpr("(C|U)-\\d{2}-\\d{3,4}|PS\\d{2}-\\d{3,4}", stem))
  if (length(m) > 0) return(m[1])
  # Strip GSL-XX-XXXX_ run prefix (e.g. "GSL-EM-4403_1894" -> "1894")
  m2 <- sub("^GSL-[A-Z]{2}-\\d+_", "", stem)
  if (m2 != stem) return(m2)
  stem
}

# Rename columns for one tumor sample within a somatic file.
# target_tn_col: the specific "T-N VAF (ID)" column to use as tumor (for trios).
# timepoint_map: named vector mapping tumor IDs → timepoint labels.
# Returns: list(df = renamed df, timepoint = "Dx"/"Rel"/"Tumor"/"blasts")
normalize_somatic_cols <- function(df, patient_id, target_tn_col = NULL,
                                   timepoint_map = NULL) {
  cols <- names(df)

  # --- Format A: PATIENT_ID_TIMEPOINT Alt Percentage ---
  tumor_alt_col <- grep(
    paste0("^", patient_id, "_(Dx|Rel|Tumor|blasts) Alt Percentage$"), cols, value = TRUE
  )

  if (length(tumor_alt_col) > 0) {
    timepoint <- sub(paste0("^", patient_id, "_(.+) Alt Percentage$"), "\\1", tumor_alt_col[1])

    tumor_ad_cols   <- grep(paste0("^", patient_id, "_(Dx|Rel|Tumor|blasts) AD Total"), cols, value = TRUE)
    germ_alt_cols   <- grep(paste0("^", patient_id, "_germline Alt Percentage$"), cols, value = TRUE)
    germ_ad_cols    <- grep(paste0("^", patient_id, "_germline AD Total"), cols, value = TRUE)
    tn_vaf_col      <- grep("^T-N VAF", cols, value = TRUE)
    somatic_flt_col <- grep(paste0("^", patient_id, "_.+ Somatic Filter$"), cols, value = TRUE)

    renames <- c(alt_percentage = tumor_alt_col[1])
    if (length(tumor_ad_cols)   > 0) renames["ad_total"]          <- tumor_ad_cols[1]
    if (length(germ_alt_cols)   > 0) renames["germline_alt_pct"]  <- germ_alt_cols[1]
    if (length(germ_ad_cols)    > 0) renames["germline_ad_total"] <- germ_ad_cols[1]
    if (length(tn_vaf_col)      > 0) renames["tn_vaf"]            <- tn_vaf_col[1]
    if (length(somatic_flt_col) > 0) renames["somatic_filter"]    <- somatic_flt_col[1]

    geno_unhandled <- setdiff(grep(paste0("^", patient_id, "_"), cols, value = TRUE), renames)
    for (gc in geno_unhandled) {
      std <- sub(paste0("^", patient_id, "_germline"), "germline", gc)
      std <- sub(paste0("^", patient_id, "_", timepoint, "(?=[ _]|$)"), "tumor", std, perl = TRUE)
      std <- gsub("[ /()]", "_", std)
      std <- sub("_+$", "", std)
      if (!std %in% names(renames)) renames[std] <- gc
    }

  } else {
    # --- Format B: tumor ID from "T-N VAF (TUMOR_ID)" ---
    all_tn_cols <- grep("^T-N VAF", cols, value = TRUE)
    if (length(all_tn_cols) == 0) stop("No tumor Alt Percentage column found and no T-N VAF column")

    tn_vaf_col <- if (!is.null(target_tn_col)) target_tn_col else all_tn_cols[1]
    tumor_id   <- sub("^T-N VAF \\((.+)\\)$", "\\1", tn_vaf_col)
    tumor_id_re <- gsub("([.+^${}()|\\[\\]\\\\])", "\\\\\\1", tumor_id)

    # Timepoint: check map first, then infer from ID suffix
    if (!is.null(timepoint_map) && tumor_id %in% names(timepoint_map)) {
      timepoint <- timepoint_map[[tumor_id]]
    } else {
      timepoint <- "Dx"
      if (grepl("rel|relapse", tumor_id, ignore.case = TRUE))                          timepoint <- "Rel"
      if (grepl("(?<![a-z])dx|diag|de.?novo", tumor_id, ignore.case = TRUE, perl = TRUE)) timepoint <- "Dx"
      if (grepl("blast", tumor_id, ignore.case = TRUE))                                timepoint <- "blasts"
      if (grepl("tumor",  tumor_id, ignore.case = TRUE))                               timepoint <- "Tumor"
    }

    all_alt_cols   <- grep("Alt Percent", cols, value = TRUE)
    tumor_alt_col  <- grep(paste0("^", tumor_id_re, " Alt Percent"), cols, value = TRUE)
    other_alt_cols <- setdiff(all_alt_cols, tumor_alt_col)
    tumor_ad_cols  <- grep(paste0("^", tumor_id_re, " AD Total"), cols, value = TRUE)
    # germline AD Total: exclude tumor cols and duplicate ...N cols
    other_ad_cols  <- setdiff(
      grep("AD Total", cols, value = TRUE),
      grep(paste0("^", tumor_id_re, " AD Total"), cols, value = TRUE)
    )
    other_ad_cols  <- other_ad_cols[!grepl("\\.\\.\\.\\d+$", other_ad_cols)]
    if (length(other_ad_cols) == 0)
      other_ad_cols <- setdiff(grep("AD Total", cols, value = TRUE),
                               grep(paste0("^", tumor_id_re, " AD Total"), cols, value = TRUE))
    somatic_flt_col <- grep(paste0("^", tumor_id_re, " Somatic Filter"), cols, value = TRUE)

    if (length(tumor_alt_col) == 0) stop(sprintf("No Alt Percentage col for tumor '%s'", tumor_id))

    renames <- c(alt_percentage = tumor_alt_col[1], tn_vaf = tn_vaf_col)
    if (length(tumor_ad_cols)   > 0) renames["ad_total"]          <- tumor_ad_cols[1]
    if (length(other_alt_cols)  > 0) renames["germline_alt_pct"]  <- other_alt_cols[1]
    if (length(other_ad_cols)   > 0) renames["germline_ad_total"] <- other_ad_cols[1]
    if (length(somatic_flt_col) > 0) renames["somatic_filter"]    <- somatic_flt_col[1]
  }

  df <- rename(df, !!!renames)
  dup_cols <- grep("\\.\\.\\.\\d+$", names(df), value = TRUE)
  df <- select(df, -any_of(dup_cols))

  list(df = df, timepoint = timepoint)
}

normalize_germline_cols <- function(df, patient_id) {
  cols <- names(df)

  # Format A: PATIENT_ID_germline_Alt_Percent
  germ_alt_col <- grep(paste0("^", patient_id, "_germline_?Alt_?Percent"), cols, value = TRUE)

  # Format B fallback: any col with _Alt_Percent
  if (length(germ_alt_col) == 0) germ_alt_col <- grep("_Alt_Percent", cols, value = TRUE)
  if (length(germ_alt_col) == 0) stop("No germline Alt% column found")

  germ_ad_cols <- grep(paste0("^", patient_id, "_germline AD Total"), cols, value = TRUE)
  if (length(germ_ad_cols) == 0) {
    germ_ad_cols <- grep("AD Total", cols, value = TRUE)
    germ_ad_cols <- germ_ad_cols[!grepl("\\.\\.\\.\\d+$", germ_ad_cols)]
    if (length(germ_ad_cols) == 0) germ_ad_cols <- grep("AD Total", cols, value = TRUE)
  }

  renames <- c(alt_percentage = germ_alt_col[1])
  if (length(germ_ad_cols) > 0) renames["ad_total"] <- germ_ad_cols[1]

  geno_unhandled <- setdiff(grep(paste0("^", patient_id, "_"), cols, value = TRUE), renames)
  for (gc in geno_unhandled) {
    std <- sub(paste0("^", patient_id, "_germline"), "germline", gc)
    std <- sub(paste0("^", patient_id, "_"), "", std)
    std <- gsub("[ /()]", "_", std)
    std <- sub("_+$", "", std)
    if (!std %in% names(renames)) renames[std] <- gc
  }

  df <- rename(df, !!!renames)
  dup_cols <- grep("\\.\\.\\.\\d+$", names(df), value = TRUE)
  df <- select(df, -any_of(dup_cols))
  df
}

# ---- File readers ----

read_somatic_file <- function(filepath, seq_run, artifact_genes, timepoint_map = NULL) {
  stem       <- sub("(_mutect2_somatic(_re)?|_retry_somatic)\\.xlsx$", "", basename(filepath))
  patient_id <- extract_patient_id(stem)
  raw        <- suppressMessages(read_xlsx(filepath))
  n_raw      <- nrow(raw)

  # Determine how many tumors are in this file (1 = pair, 2 = trio)
  tn_cols <- grep("^T-N VAF", names(raw), value = TRUE)

  # Format A files have 0 T-N VAF cols; Format B/trio have 1 or 2
  targets <- if (length(tn_cols) >= 2) tn_cols else list(NULL)

  all_dfs <- compact(lapply(targets, function(target_tn) {
    result <- tryCatch(
      normalize_somatic_cols(raw, patient_id, target_tn_col = target_tn,
                             timepoint_map = timepoint_map),
      error = function(e) stop(sprintf("[%s somatic] %s", patient_id, e$message))
    )

    result$df %>%
      mutate(
        patient_id        = patient_id,
        sample_type       = result$timepoint,
        call_type         = "Somatic",
        seq_run           = seq_run,
        hypermutator_flag = patient_id %in% HYPERMUTATORS,
        Filter            = as.character(Filter),
        .before = 1
      ) %>%
      filter(
        Filter == "PASS",
        !is.na(ad_total),
        as.numeric(ad_total) > 20,
        as.numeric(alt_percentage) >= 0.02,
        !Gene %in% artifact_genes
      )
  }))

  df <- bind_rows(all_dfs)
  attr(df, "n_raw") <- n_raw
  df
}

read_germline_file <- function(filepath, seq_run, artifact_genes) {
  stem       <- sub("_germline(_re)?\\.xlsx$", "", basename(filepath))
  patient_id <- extract_patient_id(stem)
  raw        <- suppressMessages(read_xlsx(filepath))
  n_raw      <- nrow(raw)

  df <- tryCatch(
    normalize_germline_cols(raw, patient_id),
    error = function(e) stop(sprintf("[%s germline] %s", patient_id, e$message))
  )

  df <- df %>%
    mutate(
      patient_id        = patient_id,
      sample_type       = "Germline",
      call_type         = "Germline",
      seq_run           = seq_run,
      hypermutator_flag = patient_id %in% HYPERMUTATORS,
      .before = 1
    ) %>%
    filter(
      !is.na(ad_total),
      as.numeric(ad_total) > 20,
      as.numeric(alt_percentage) >= 0.02,
      !Gene %in% artifact_genes
    )

  attr(df, "n_raw") <- n_raw
  df
}

# ---- Per-run builder ----

build_run_combined <- function(run_dir, seq_run, artifact_genes, out_dir = NULL,
                               timepoint_map = NULL) {
  run_dir <- file.path(BASE_DIR, run_dir)
  somatic_files  <- list.files(run_dir, pattern = "_mutect2_somatic(_re)?\\.xlsx$|_retry_somatic\\.xlsx$", full.names = TRUE, recursive = TRUE)
  germline_files <- list.files(run_dir, pattern = "_germline(_re)?\\.xlsx$",        full.names = TRUE, recursive = TRUE)
  # Exclude any previously generated combined files from being picked up
  somatic_files  <- somatic_files[!grepl("_combined\\.xlsx$", somatic_files)]
  germline_files <- germline_files[!grepl("_combined\\.xlsx$", germline_files)]

  cat(sprintf("\n=== %s ===\n", seq_run))
  cat(sprintf("Input dir: %s\n", run_dir))
  cat(sprintf("Files: %d somatic, %d germline\n", length(somatic_files), length(germline_files)))

  read_somatic_logged <- function(f) {
    tryCatch({
      df <- read_somatic_file(f, seq_run, artifact_genes, timepoint_map)
      cat(sprintf("  %-55s  raw=%d  pass=%d\n", basename(f), attr(df, "n_raw"), nrow(df)))
      df
    }, error = function(e) {
      cat(sprintf("  ERROR %-55s  %s\n", basename(f), e$message))
      NULL
    })
  }

  read_germline_logged <- function(f) {
    tryCatch({
      df <- read_germline_file(f, seq_run, artifact_genes)
      cat(sprintf("  %-55s  raw=%d  pass=%d\n", basename(f), attr(df, "n_raw"), nrow(df)))
      df
    }, error = function(e) {
      cat(sprintf("  ERROR %-55s  %s\n", basename(f), e$message))
      NULL
    })
  }

  somatic_list  <- compact(map(somatic_files,  read_somatic_logged))
  germline_list <- compact(map(germline_files, read_germline_logged))
  all_list <- c(somatic_list, germline_list)

  if (length(all_list) == 0) {
    warning("No files processed for run: ", seq_run)
    return(NULL)
  }

  col_sets   <- map(all_list, names)
  all_cols   <- Reduce(union, col_sets)
  extra_cols <- setdiff(all_cols, Reduce(intersect, col_sets))
  if (length(extra_cols) > 0) {
    cat(sprintf("  Schema note: %d columns vary across files (will be NA where absent):\n", length(extra_cols)))
    cat(sprintf("    %s\n", paste(head(extra_cols, 10), collapse = ", ")))
  }

  combined <- bind_rows(map(all_list, \(df) mutate(df, across(everything(), as.character))))
  combined <- combined %>% mutate(
    alt_percentage    = as.numeric(alt_percentage),
    ad_total          = as.numeric(ad_total),
    hypermutator_flag = as.logical(hypermutator_flag)
  )

  cat(sprintf("\nSummary for %s:\n", seq_run))
  cat(sprintf("  Total PASS rows:  %d\n", nrow(combined)))
  cat(sprintf("  Somatic PASS:     %d\n", sum(combined$call_type == "Somatic")))
  cat(sprintf("  Germline PASS:    %d\n", sum(combined$call_type == "Germline")))
  cat(sprintf("  Patients:         %d\n", n_distinct(combined$patient_id)))
  cat(sprintf("  Timepoints:       %s\n", paste(sort(unique(combined$sample_type)), collapse = ", ")))
  cat(sprintf("  Hypermutator flag: %d rows flagged\n", sum(combined$hypermutator_flag)))

  save_dir <- if (!is.null(out_dir)) file.path(BASE_DIR, out_dir) else OUTPUT_DIR
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(save_dir, paste0(seq_run, "_combined.xlsx"))
  write_xlsx(combined, out_path)
  cat(sprintf("  Saved: %s\n", out_path))

  combined
}

# ---- Main ----

all_runs <- map(RUN_CONFIG, function(cfg) {
  build_run_combined(cfg$dir, cfg$seq_run, artifact_genes, cfg$out_dir, cfg$timepoint_map)
})
all_runs <- compact(all_runs)

if (length(all_runs) > 1) {
  master <- bind_rows(all_runs)
  write_xlsx(master, MASTER_OUT)
  cat(sprintf("\nMaster variants saved: %s\n", MASTER_OUT))
  cat(sprintf("  Total rows: %d across %d runs\n", nrow(master), length(all_runs)))
} else if (length(all_runs) == 1) {
  cat(sprintf("\nOnly one run processed — master_variants.xlsx not written yet.\n"))
  cat(sprintf("Add more runs to RUN_CONFIG and re-run to build master.\n"))
}
