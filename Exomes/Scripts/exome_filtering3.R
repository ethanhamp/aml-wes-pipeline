library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggrepel)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

local_env <- new.env()
source(file.path(BASE_DIR, "General/jeff.genes.txt"), local = local_env)
jeff.genes <- local_env$jeff.genes

EXCLUDED_FILE <- file.path(BASE_DIR, "config/excluded_samples.txt")
EXCLUDED_SAMPLES <- if (file.exists(EXCLUDED_FILE)) {
  lines <- readLines(EXCLUDED_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else {
  character(0)
}

all_samples_clean <- read_excel(file.path(BASE_DIR, "Exomes/Datasets/3801_combined_pass_withBT.xlsx"))
combined_samples_clean <- read_excel(file.path(BASE_DIR, "Exomes/Datasets/combined_exomes_relapse.xlsx"))
combined <- merge(all_samples_clean, combined_samples_clean, all = TRUE)

combined <- combined[combined$`AD Total` > 20, ]

combined <- combined[!combined$`Sample ID` %in% EXCLUDED_SAMPLES, ]
combined <- combined[!(combined$Group %in% c("11","4")), ]
combined <- combined[!combined$Gene %in% jeff.genes, ]

# all_samples_clean <- all_samples_clean[all_samples_clean$`Alt Percentage` >0.02, ]

variant_key <- "MANE HGVS"

clean <- combined %>%
  mutate(
    Timepoint = tolower(trimws(as.character(Timepoint))),     # "diagnosis"/"relapse"
    Group     = as.character(Group),
    HGVS      = str_squish(as.character(`MANE HGVS`)),
    # optional: normalize transcript version so NM_XXXX.6:c.xxx == NM_XXXX:c.xxx
    HGVS      = str_replace(HGVS, "(NM_\\d+)\\.\\d+(:.*)", "\\1\\2"),
    VAF_raw   = as.character(`Alt Percentage`),
    VAF       = parse_number(VAF_raw),
    VAF       = ifelse(str_detect(VAF_raw, "%"), VAF/100, VAF)
  ) %>%
  filter(Timepoint %in% c("germline","diagnosis","relapse"),
         !is.na(HGVS) & HGVS != "")

# Collapse duplicates per (Group × HGVS × Timepoint)
collapsed <- clean %>%
  group_by(Group, Gene, HGVS, Timepoint, `AD Total`) %>%
  summarise(VAF = max(VAF, na.rm = TRUE), .groups = "drop")

write.csv(collapsed, file.path(BASE_DIR, "Exomes/Datasets/filtered_collapsed.csv"), row.names = FALSE)

# Pivot to wide and compute delta
# pairs_hgvs <- collapsed %>%
#   pivot_wider(
#     id_cols    = c(Group, Gene, HGVS),
#     names_from = Timepoint,                 # -> diagnosis, relapse
#     values_from = c(VAF, `AD Total`),
#     values_fill = NA_real_
#   ) %>%
#   filter(!is.na(diagnosis) & !is.na(relapse)) %>%
#   mutate(delta = relapse - diagnosis)
# 
# pairs_hgvs <- pairs_hgvs %>%
#   mutate(Group = factor(Group, levels = as.character(sort(as.numeric(unique(Group))))))

library(dplyr)
library(tidyr)
library(stringr)

pairs_hgvs <- collapsed %>%
  # Reconstruct counts from VAF and total depth (DP)
  mutate(
    DP  = `AD Total`,
    ALT = round(VAF * DP),
    REF = pmax(DP - ALT, 0)
  ) %>%
  # If there can be duplicate rows per (Group, Gene, HGVS, Timepoint),
  # aggregate them first
  group_by(Group, Gene, HGVS, Timepoint) %>%
  summarise(
    ALT = sum(ALT, na.rm = TRUE),
    REF = sum(REF, na.rm = TRUE),
    DP  = sum(DP,  na.rm = TRUE),
    VAF = ifelse(ALT + REF > 0, ALT / (ALT + REF), NA_real_),
    .groups = "drop"
  ) %>%
  # Wide format: create per-timepoint columns
  pivot_wider(
    id_cols    = c(Group, Gene, HGVS),
    names_from = Timepoint,                      # -> diagnosis, relapse
    values_from = c(VAF, ALT, REF, DP),
    names_sep  = "_",
    values_fill = NA_real_
  ) %>%
  # Keep only variants seen at both timepoints (optional)
  filter(!is.na(VAF_diagnosis), !is.na(VAF_relapse)) %>%
  # Useful delta
  mutate(delta = VAF_relapse - VAF_diagnosis) %>%
  # Make Group an ordered factor without coercion issues
  mutate(Group = factor(Group, levels = str_sort(unique(Group), numeric = TRUE)))

pairs_hgvs <- pairs_hgvs %>%
  rowwise() %>%
  mutate(
    p_fisher = fisher.test(
      matrix(c(ALT_diagnosis, REF_diagnosis, ALT_relapse, REF_relapse),
             nrow = 2, byrow = TRUE)
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_fisher, method = "BH"))

write.csv(pairs_hgvs, file.path(BASE_DIR, "Exomes/Datasets/combined_exomes_withpvalues.csv"), row.names = FALSE)

sig_variants <- pairs_hgvs %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  mutate(
    delta = if (!"delta" %in% names(.)) VAF_relapse - VAF_diagnosis else delta,
    direction = case_when(
      delta >  0 ~ "increase_at_relapse",
      delta <  0 ~ "decrease_at_relapse",
      TRUE        ~ "no_change"
    ),
    abs_delta = abs(delta),
    # Optional: log2 fold-change with a tiny pseudocount to avoid /0
    log2FC = log2( pmax(VAF_relapse, 1e-6) / pmax(VAF_diagnosis, 1e-6) )
  ) %>%
  arrange(FDR, desc(abs_delta))

pi_version <- pairs_hgvs %>%
  # STEP 1: keep rows where either Dx OR Relapse VAF > 2%
  filter(
    (!is.na(VAF_diagnosis) & VAF_diagnosis > 0.02) |
      (!is.na(VAF_relapse)   & VAF_relapse   > 0.02)
  ) %>%
  # STEP 2: high VAF label
  mutate(
    high.VAF.variant = case_when(
      (
        ((!is.na(VAF_diagnosis) & VAF_diagnosis > 0.30) |
           (!is.na(VAF_relapse)   & VAF_relapse   > 0.30)) &
          (is.na(VAF_germline) | VAF_germline < 0.05)
      ) ~ "y",
      TRUE ~ "n"
    ),
    # STEP 3: low VAF label
    low.VAF.variant = case_when(
      (
        (
          (!is.na(VAF_diagnosis) & VAF_diagnosis >= 0.02 & VAF_diagnosis < 0.30) |
            (!is.na(VAF_relapse)   & VAF_relapse   >= 0.02 & VAF_relapse   < 0.30)
        ) &
          (is.na(VAF_germline) | VAF_germline < 0.01)
      ) ~ "y",
      TRUE ~ "n"
    )
  )

write.csv(pi_version, file.path(BASE_DIR, "Exomes/Datasets/combined_exomes_filtered.csv"), row.names = FALSE)
