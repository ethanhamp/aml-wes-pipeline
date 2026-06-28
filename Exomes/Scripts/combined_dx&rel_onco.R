library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggrepel)
library(ComplexHeatmap)
library(grid)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
ARTIFACT_GENE_FILE <- file.path(BASE_DIR, "General/artifact_genes.txt")
EXCLUDED_FILE <- file.path(BASE_DIR, "config/excluded_samples.txt")
HYPERMUTATOR_FILE <- file.path(BASE_DIR, "config/hypermutators.txt")

artifact_genes <- if (file.exists(ARTIFACT_GENE_FILE)) {
  lines <- readLines(ARTIFACT_GENE_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else { character(0) }

EXCLUDED_SAMPLES <- if (file.exists(EXCLUDED_FILE)) {
  lines <- readLines(EXCLUDED_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else { character(0) }

HYPERMUTATORS <- if (file.exists(HYPERMUTATOR_FILE)) {
  lines <- readLines(HYPERMUTATOR_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else { character(0) }

all_samples_clean <- read_excel(file.path(BASE_DIR, "Exomes/Datasets/3801_combined_pass_withBT.xlsx"))
combined_samples_clean <- read_excel(file.path(BASE_DIR, "Exomes/Datasets/combined_exomes_relapse.xlsx"))
combined <- merge(all_samples_clean, combined_samples_clean, all = TRUE)

combined <- combined[combined$`AD Total` >20, ]

combined <- combined[!combined$`Sample ID` %in% EXCLUDED_SAMPLES, ]
combined <- combined[!combined$Group %in% HYPERMUTATORS, ]

combined <- combined[!combined$Gene %in% artifact_genes, ]

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
  group_by(Group, Gene, HGVS, Timepoint) %>%
  summarise(VAF = max(VAF, na.rm = TRUE), .groups = "drop")

genes <- c("NPM1","FLT3","NRAS","DNMT3A", "ARHG", "PHIP")

counts_long <- collapsed %>%
  filter(Gene %in% genes,
         Timepoint %in% c("diagnosis","relapse"),
         VAF >= 0.02) %>%
  distinct(Gene, Group) %>%
  count(Gene, name = "patients_any")

counts_long


# Pivot to wide and compute delta
pairs_hgvs <- collapsed %>%
  pivot_wider(
    id_cols    = c(Group, Gene, HGVS),
    names_from = Timepoint,                 # -> diagnosis, relapse
    values_from = VAF,
    values_fill = NA_real_
  ) %>%
  filter(!is.na(diagnosis) & !is.na(relapse)) %>%
  mutate(delta = relapse - diagnosis)

pairs_hgvs <- pairs_hgvs %>%
  mutate(Group = factor(Group, levels = as.character(sort(as.numeric(unique(Group))))))

pi_version <- pairs_hgvs %>%
  # STEP 1: keep rows where either Dx OR Relapse VAF > 2%
  filter(
    (!is.na(diagnosis) & diagnosis > 0.02) |
      (!is.na(relapse)   & relapse   > 0.02)
  ) %>%
  # STEP 2: high VAF label
  mutate(
    high.VAF.variant = case_when(
      (
        ((!is.na(diagnosis) & diagnosis > 0.30) |
           (!is.na(relapse)   & relapse   > 0.30)) &
          (is.na(germline) | germline < 0.05)
      ) ~ "y",
      TRUE ~ "n"
    ),
    # STEP 3: low VAF label
    low.VAF.variant = case_when(
      (
        (
          (!is.na(diagnosis) & diagnosis >= 0.02 & diagnosis < 0.30) |
            (!is.na(relapse)   & relapse   >= 0.02 & relapse   < 0.30)
        ) &
          (is.na(germline) | germline < 0.01)
      ) ~ "y",
      TRUE ~ "n"
    )
  )

# 1) Keep genes/samples that pass your Step 1 (>2% in Dx or Rel)
mat_base <- pi_version %>%
  filter(pmax(coalesce(diagnosis, 0), coalesce(relapse, 0)) > 0.02)

# 2) Build per-timepoint alteration labels using PI thresholds
label_call <- function(tp, vaf, germline) {
  if (is.na(vaf)) return(NA_character_)
  if (vaf > 0.30 && (is.na(germline) || germline < 0.05)) {
    return(paste0(tp, "_high"))      # e.g., "Dx_high" / "Rel_high"
  } else if (vaf >= 0.02 && vaf < 0.30 && (is.na(germline) || germline < 0.01)) {
    return(paste0(tp, "_low"))       # e.g., "Dx_low" / "Rel_low"
  } else if (vaf >= 0.02)
  return(NA_character_)
}

mat_calls <- mat_base %>%
  mutate(
    # Dx labels
    Dx_call = case_when(
      is.na(diagnosis) ~ NA_character_,
      diagnosis > 0.30 & (is.na(germline) | germline < 0.05) ~ "Dx_high",
      diagnosis >= 0.02 & diagnosis < 0.30 & (is.na(germline) | germline < 0.01) ~ "Dx_low",
      TRUE ~ NA_character_
    ),
    # Rel labels
    Rel_call = case_when(
      is.na(relapse) ~ NA_character_,
      relapse > 0.30 & (is.na(germline) | germline < 0.05) ~ "Rel_high",
      relapse >= 0.02 & relapse < 0.30 & (is.na(germline) | germline < 0.01) ~ "Rel_low",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(c(Dx_call, Rel_call), names_to = "tp", values_to = "call") %>%
  mutate(
    TimepointCol = if_else(tp == "Dx_call", "Dx", "Rel"),
    ColID = paste0(Group, "_", TimepointCol)
  ) %>%
  filter(!is.na(call)) %>%
  group_by(Gene, ColID) %>%
  summarise(call = paste(sort(unique(call)), collapse = ";"), .groups = "drop")

# 3) Make a complete Gene x ColID grid (missing = "")
all_genes <- sort(unique(mat_calls$Gene))
all_cols  <- mat_base %>%
  transmute(ColID_Dx  = paste0(Group, "_Dx"),
            ColID_Rel = paste0(Group, "_Rel")) %>%
  pivot_longer(everything(), values_to = "ColID") %>%
  distinct(ColID) %>%
  pull(ColID) %>%
  sort()

mat_full <- tidyr::complete(
  mat_calls, Gene = all_genes, ColID = all_cols, fill = list(call = "")
) %>%
  arrange(Gene, ColID)

# 4) Convert to matrix (rows=genes, columns=sample_timepoint)
oncomat <- mat_full %>%
  pivot_wider(names_from = ColID, values_from = call) %>%
  arrange(Gene) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

# (Optional) Limit to top N genes by prevalence
keep_top <- 30
present_counts <- rowSums(oncomat != "")
oncomat <- oncomat[order(present_counts, decreasing = TRUE)[1:keep_top], , drop = FALSE]

# 5) Define colors and how to draw each alteration
col_fun <- c(
  "Dx_high"  = "#1f77b4",
  "Dx_low"   = "#aec7e8",
  "Rel_high" = "#d62728",
  "Rel_low"  = "#ff9896"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#f7f7f7", col = NA))
  },
  Dx_high = function(x, y, w, h)  grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col_fun["Dx_high"],  col = NA)),
  Dx_low  = function(x, y, w, h)  grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col_fun["Dx_low"],   col = NA)),
  Rel_high= function(x, y, w, h)  grid.circle(x, y, r = min(w,h)*0.45,   gp = gpar(fill = col_fun["Rel_high"], col = NA)),
  Rel_low = function(x, y, w, h)  grid.circle(x, y, r = min(w,h)*0.45,   gp = gpar(fill = col_fun["Rel_low"],  col = NA))
)

# 6) Order columns: group samples together (Dx next to Rel)
sample_order <- unique(str_remove(colnames(oncomat), "_(Dx|Rel)$"))
col_order <- as.vector(rbind(paste0(sample_order, "_Dx"), paste0(sample_order, "_Rel")))
col_order <- col_order[col_order %in% colnames(oncomat)]

# 7) Draw
oncoPrint(
  oncomat,
  alter_fun = alter_fun,
  col = col_fun,
  row_names_side = "left",
  pct_side = "right",
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  show_column_names = FALSE,
  column_order = col_order,
  top_annotation = HeatmapAnnotation(
    Timepoint = factor(ifelse(str_ends(colnames(oncomat), "Dx"), "Dx", "Rel"),
                       levels = c("Dx", "Rel")),
    annotation_name_side = "left"
  ),
  column_split = ifelse(str_ends(colnames(oncomat), "Dx"), "Dx", "Rel")
)

# overlaying dx and relapse
# Build exactly one Dx and one Rel label per Gene×Group
calls_one_per_tp <- mat_base %>%
  mutate(
    Dx_call  = dplyr::case_when(
      is.na(diagnosis) ~ NA_character_,
      diagnosis > 0.30 & (is.na(germline) | germline < 0.05) ~ "Dx_high",
      diagnosis >= 0.02 & diagnosis < 0.30 & (is.na(germline) | germline < 0.01) ~ "Dx_low",
      TRUE ~ NA_character_
    ),
    Rel_call = dplyr::case_when(
      is.na(relapse) ~ NA_character_,
      relapse > 0.30 & (is.na(germline) | germline < 0.05) ~ "Rel_high",
      relapse >= 0.02 & relapse < 0.30 & (is.na(germline) | germline < 0.01) ~ "Rel_low",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Gene, Group) %>%
  summarise(
    # pick ONE Dx and ONE Rel per cell, prioritizing "high" over "low"
    Dx = dplyr::case_when(
      any(Dx_call == "Dx_high") ~ "Dx_high",
      any(Dx_call == "Dx_low")  ~ "Dx_low",
      TRUE ~ NA_character_
    ),
    Rel = dplyr::case_when(
      any(Rel_call == "Rel_high") ~ "Rel_high",
      any(Rel_call == "Rel_low")  ~ "Rel_low",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  transmute(Gene, Group, call = paste(na.omit(c(Dx, Rel)), collapse = ";"))

# Complete grid and pivot: ONE column per Group
all_genes  <- sort(unique(calls_one_per_tp$Gene))
all_groups <- sort(unique(calls_one_per_tp$Group))

mat_full_onecol <- tidyr::complete(
  calls_one_per_tp, Gene = all_genes, Group = all_groups, fill = list(call = "")
) %>%
  arrange(Gene, Group)

oncomat_onecol <- mat_full_onecol %>%
  tidyr::pivot_wider(names_from = Group, values_from = call) %>%
  arrange(Gene) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

# Optional: keep top genes
keep_top <- 30
present_counts <- rowSums(oncomat_onecol != "")
oncomat_onecol <- oncomat_onecol[order(present_counts, decreasing = TRUE)[1:keep_top], , drop = FALSE]
# sanity check: no duplicate patient columns
stopifnot(anyDuplicated(colnames(oncomat_onecol)) == 0)

gap_mm <- unit(0.5, "mm")  # or 0 mm if you want no gap

ht <- oncoPrint(
  oncomat_onecol,
  name = "onc",
  alter_fun = alter_fun,    # your left/right half drawing funcs
  col = col_fun,
  column_gap = gap_mm,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  row_names_side = "left",
  pct_side = "right",
  show_column_names = TRUE
)

ht <- draw(ht)

decorate_heatmap_body("onc", {
  nc <- ncol(oncomat_onecol)
  if (nc > 1) {
    # draw at exact column boundaries: x = j + 0.5 in "native" coords
    for (j in seq_len(nc - 1)) {
      grid.segments(
        x0 = unit(j + 0.5, "native"), y0 = unit(0, "npc"),
        x1 = unit(j + 0.5, "native"), y1 = unit(1, "npc"),
        gp = gpar(col = "grey70", lwd = 0.6)
      )
    }
  }
})

