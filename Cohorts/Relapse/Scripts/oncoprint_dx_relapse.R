library(readxl)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(grid)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# ── Configuration ──────────────────────────────────────────────────────────────
VAF_THRESHOLD <- 0.02   # 2% — variant is "present" at a timepoint if VAF >= this
TOP_N_GENES   <- 20

# Color scheme: persistent = blue, lost = red/coral, gained = green
COLORS <- c(
  Persistent = "#2166AC",
  Lost       = "#D6604D",
  Gained     = "#4DAC26"
)

# Artifact genes excluded from all analyses (sourced from combine_pass_variants.R)
local_env <- new.env()
source(file.path(BASE_DIR, "General/artifact_genes.txt"), local = local_env)
artifact_genes <- local_env$artifact_genes

# ── 1. Load combined variants ──────────────────────────────────────────────────
message("Loading relapse_pass_variants_combined.xlsx ...")
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))
message(sprintf("  %d rows, %d columns", nrow(df), ncol(df)))

df |>
  filter(Timepoint %in% c("diagnosis", "relapse")) |>                                                                                                                                                     
  distinct(Pair_ID, Timepoint, Sample_ID) |>                                                                                                                                                              
  arrange(Pair_ID) 

# ── 2. Filter to somatic timepoints only; exclude artifact genes
#       (combined file may predate artifact_genes exclusion in combine script)
df_somatic <- df |>
  filter(Timepoint %in% c("diagnosis", "relapse"),
         !Gene %in% artifact_genes)

# ── 2b. Collapse all ARHG* genes into a single "ARHG" label ──────────────────
df_somatic <- df_somatic |>
  mutate(Gene = if_else(startsWith(Gene, "ARHG"), "ARHG", Gene))

df_somatic <- df_somatic %>% 
  filter(!Patient_Group %in% c('Patient_23', 'Patient_25'))

# ── 3. Deduplicate: one row per Pair_ID × Gene × HGVS × Timepoint
#       (take max VAF when duplicate entries exist)
variants_dedup <- df_somatic |>
  group_by(Pair_ID, Gene, HGVS, Timepoint) |>
  summarise(vaf = max(`T-N VAF`, na.rm = TRUE), .groups = "drop")

# ── 4. Pivot to wide: one row per Pair_ID × Gene × HGVS with Dx and Rel VAFs
variants_wide <- variants_dedup |>
  pivot_wider(
    names_from  = Timepoint,
    values_from = vaf,
    values_fill = 0       # 0 VAF if variant not seen at that timepoint
  )

# Ensure both columns exist even if all variants are one-sided
if (!"diagnosis" %in% names(variants_wide)) variants_wide$diagnosis <- 0
if (!"relapse"   %in% names(variants_wide)) variants_wide$relapse   <- 0

# ── 5. Classify each variant as Persistent / Lost / Gained ────────────────────
variants_classified <- variants_wide |>
  mutate(
    present_dx  = (diagnosis >= VAF_THRESHOLD),
    present_rel = (relapse   >= VAF_THRESHOLD),
    status = case_when(
      present_dx  & present_rel  ~ "Persistent",
      present_dx  & !present_rel ~ "Lost",
      !present_dx & present_rel  ~ "Gained",
      TRUE ~ NA_character_        # not present at either timepoint — exclude
    )
  ) |>
  filter(!is.na(status))

message(sprintf("  Classified variants — Persistent: %d | Lost: %d | Gained: %d",
                sum(variants_classified$status == "Persistent"),
                sum(variants_classified$status == "Lost"),
                sum(variants_classified$status == "Gained")))

# ── 6. Top N genes by unique patients affected (across all statuses) ───────────
top_genes <- variants_classified |>
  group_by(Gene) |>
  summarise(n_pts = n_distinct(Pair_ID), .groups = "drop") |>
  arrange(desc(n_pts)) |>
  slice_head(n = TOP_N_GENES) |>
  pull(Gene)

message("Top ", TOP_N_GENES, " genes: ", paste(top_genes, collapse = ", "))

# ── 7. Collapse to one status per Pair_ID × Gene
#       Priority: Persistent > Gained > Lost
#       (e.g., if a patient has one persistent + one lost variant in a gene,
#        the gene is shown as Persistent — most biologically important)
vs_top <- variants_classified |>
  filter(Gene %in% top_genes) |>
  group_by(Pair_ID, Gene) |>
  summarise(
    gene_status = case_when(
      any(status == "Persistent") ~ "Persistent",
      any(status == "Gained")     ~ "Gained",
      any(status == "Lost")       ~ "Lost"
    ),
    .groups = "drop"
  )

# ── 8. Build full gene × patient grid (fill empty = "None") ───────────────────
all_patients <- sort(unique(variants_classified$Pair_ID))

full_grid <- expand.grid(
  Gene    = top_genes,
  Pair_ID = all_patients,
  stringsAsFactors = FALSE
) |>
  left_join(vs_top, by = c("Gene", "Pair_ID")) |>
  replace_na(list(gene_status = "None"))

# ── 9. Build character matrices for Dx and Relapse panels ─────────────────────
# Dx panel shows: Persistent (present at both) and Lost (present only at Dx)
# Relapse panel shows: Persistent and Gained (present only at Relapse)

make_mat <- function(full_grid, patients, genes, keep_statuses) {
  mat <- matrix(
    "",
    nrow = length(genes),
    ncol = length(patients),
    dimnames = list(genes, patients)
  )
  sub <- full_grid[full_grid$gene_status %in% keep_statuses, ]
  for (i in seq_len(nrow(sub))) {
    g <- sub$Gene[i]
    p <- sub$Pair_ID[i]
    if (g %in% genes && p %in% patients) {
      mat[g, p] <- sub$gene_status[i]
    }
  }
  mat
}

mat_dx  <- make_mat(full_grid, all_patients, top_genes, c("Persistent", "Lost"))
mat_rel <- make_mat(full_grid, all_patients, top_genes, c("Persistent", "Gained"))

# ── 10. Determine display order ────────────────────────────────────────────────
# Genes: by total affected patients (most on top) — computed from full_grid
gene_freq <- full_grid |>
  filter(gene_status != "None") |>
  count(Gene) |>
  arrange(desc(n))
gene_order <- match(gene_freq$Gene, top_genes)  # index into top_genes

# Patients: sort by number of Dx mutations (most mutated patients on left)
dx_burden   <- colSums(mat_dx != "")
patient_order <- order(dx_burden, decreasing = TRUE)

# Apply ordering
mat_dx  <- mat_dx[ gene_order, patient_order]
mat_rel <- mat_rel[gene_order, patient_order]

n_pts  <- ncol(mat_dx)
n_genes <- nrow(mat_dx)
message(sprintf("Matrix: %d genes × %d patients", n_genes, n_pts))

# ── 11. Define alter_fun and colors for ComplexHeatmap oncoPrint ───────────────
# alter_fun lists determine how each alteration type is drawn inside each cell.
# Both panels share the "Persistent" blue; only their second type differs.

make_alter_fun <- function(types, colors) {
  bg_fun <- list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#E0E0E0", col = NA))
    }
  )
  type_funs <- lapply(types, function(tp) {
    col <- colors[tp]
    function(x, y, w, h) {
      grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col, col = NA))
    }
  })
  names(type_funs) <- types
  c(bg_fun, type_funs)
}

alter_fun_dx  <- make_alter_fun(c("Persistent", "Lost"),    COLORS)
alter_fun_rel <- make_alter_fun(c("Persistent", "Gained"),  COLORS)

col_dx  <- COLORS[c("Persistent", "Lost")]
col_rel <- COLORS[c("Persistent", "Gained")]

# ── 12. Build oncoPrint objects ────────────────────────────────────────────────
# column_order and row_order are both fixed (seq 1:n) so both panels use the
# identical patient and gene ordering we set above.

op_dx <- oncoPrint(
  mat_dx,
  alter_fun           = alter_fun_dx,
  col                 = col_dx,
  show_pct            = TRUE,
  pct_side            = "right",
  pct_gp              = gpar(fontsize = 8),
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 10, fontface = "italic"),
  column_title        = "Diagnosis",
  column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
  show_column_names   = FALSE,
  column_order        = seq_len(n_pts),     # lock patient order
  row_order           = seq_len(n_genes),   # lock gene order
  remove_empty_columns = FALSE,             # keep all patients aligned between panels
  remove_empty_rows    = FALSE,
  show_heatmap_legend  = FALSE              # suppress auto-legend; we draw a unified one below
)

op_rel <- oncoPrint(
  mat_rel,
  alter_fun           = alter_fun_rel,
  col                 = col_rel,
  show_pct            = TRUE,
  pct_side            = "right",
  pct_gp              = gpar(fontsize = 8),
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 10, fontface = "italic"),
  column_title        = "Relapse",
  column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
  show_column_names   = FALSE,
  column_order        = seq_len(n_pts),
  row_order           = seq_len(n_genes),
  remove_empty_columns = FALSE,
  remove_empty_rows    = FALSE,
  show_heatmap_legend  = FALSE
)

# ── 13. Unified legend ─────────────────────────────────────────────────────────
# Build one combined legend for all three statuses so both panels share it.
unified_lgd <- Legend(
  labels     = c("Persistent", "Lost", "Gained"),
  legend_gp  = gpar(fill = COLORS[c("Persistent", "Lost", "Gained")]),
  title      = "Mutation Status",
  title_gp   = gpar(fontsize = 10, fontface = "bold"),
  labels_gp  = gpar(fontsize = 9),
  grid_height = unit(5, "mm"),
  grid_width  = unit(5, "mm")
)

# ── 14. Draw and save ──────────────────────────────────────────────────────────
out_dir  <- "Cohorts/Relapse"
out_file <- file.path(out_dir, "oncoprint_dx_vs_relapse_top20.pdf")

pdf(out_file, width = 18, height = 7)
draw(
  op_dx + op_rel,
  ht_gap                 = unit(10, "mm"),
  padding                = unit(c(5, 5, 10, 5), "mm"),
  annotation_legend_list = list(unified_lgd)
)
dev.off()

message("\nSaved: ", out_file)
