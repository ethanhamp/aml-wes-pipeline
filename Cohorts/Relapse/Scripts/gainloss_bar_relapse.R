library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

VAF_THRESHOLD <- 0.02
TOP_N_GENES   <- 20

local_env <- new.env()
source(file.path(BASE_DIR, "General/jeff.genes.txt"), local = local_env)
jeff.genes <- local_env$jeff.genes

# ── 1. Load and filter ─────────────────────────────────────────────────────────
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))

df_somatic <- df |>
  filter(
    Timepoint %in% c("diagnosis", "relapse"),
    !Gene %in% jeff.genes,
    !Patient_Group %in% c("Patient_23", "Patient_25")
  ) |>
  mutate(Gene = if_else(startsWith(Gene, "ARHG"), "ARHG", Gene))

# ── 2. Deduplicate and pivot to wide VAFs ──────────────────────────────────────
variants_dedup <- df_somatic |>
  group_by(Pair_ID, Gene, HGVS, Timepoint) |>
  summarise(vaf = max(`T-N VAF`, na.rm = TRUE), .groups = "drop")

variants_wide <- variants_dedup |>
  pivot_wider(names_from = Timepoint, values_from = vaf, values_fill = 0)

if (!"diagnosis" %in% names(variants_wide)) variants_wide$diagnosis <- 0
if (!"relapse"   %in% names(variants_wide)) variants_wide$relapse   <- 0

# ── 3. Classify each variant ───────────────────────────────────────────────────
variants_classified <- variants_wide |>
  mutate(
    present_dx  = diagnosis >= VAF_THRESHOLD,
    present_rel = relapse   >= VAF_THRESHOLD,
    status = case_when(
      present_dx  & present_rel  ~ "persistent",
      present_dx  & !present_rel ~ "loss",
      !present_dx & present_rel  ~ "gain",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(status))

# ── 4. Collapse to one status per Pair_ID × Gene (Priority: persistent > gain > loss)
vs_gene <- variants_classified |>
  group_by(Pair_ID, Gene) |>
  summarise(
    status = case_when(
      any(status == "persistent") ~ "persistent",
      any(status == "gain")       ~ "gain",
      any(status == "loss")       ~ "loss"
    ),
    .groups = "drop"
  )

# ── 5. Select top N genes by number of affected patients ──────────────────────
top_genes <- vs_gene |>
  count(Gene, name = "n_pts") |>
  arrange(desc(n_pts)) |>
  slice_head(n = TOP_N_GENES) |>
  pull(Gene)

# ── 6. Compute % of affected patients per gene per status ─────────────────────
# Denominator = total unique patients with any variant in that gene
gene_counts <- vs_gene |>
  filter(Gene %in% top_genes) |>
  group_by(Gene) |>
  mutate(total_pts = n_distinct(Pair_ID)) |>
  ungroup() |>
  count(Gene, status, total_pts) |>
  mutate(pct = n / total_pts * 100)

# ── 7. Gene order: by total affected patients (most at top of plot) ────────────
gene_order <- vs_gene |>
  filter(Gene %in% top_genes) |>
  count(Gene, name = "n_pts") |>
  arrange(n_pts) |>   # ascending so fct_inorder puts most-frequent at top
  pull(Gene)

gene_counts <- gene_counts |>
  mutate(
    Gene   = factor(Gene, levels = gene_order),
    status = factor(status, levels = c("persistent", "loss", "gain"))
  )

# ── 8. Plot ────────────────────────────────────────────────────────────────────
COLORS <- c(
  persistent = "#2166AC",   # blue  — canonical palette (vaf_scatter, oncoprint)
  loss       = "#D6604D",   # red
  gain       = "#4DAC26"    # green
)

p <- ggplot(gene_counts, aes(x = pct, y = Gene, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = COLORS,
    labels = c("persistent", "loss", "gain"),
    name   = "Dx-Relapse category"
  ) +
  scale_x_reverse(
    limits = c(100, 0),
    breaks = seq(100, 0, by = -25),
    labels = seq(100, 0, by = -25),
    expand = expansion(mult = c(0.01, 0))
  ) +
  labs(
    title = "Gain vs Loss vs Persistent (Dx \u2192 Relapse)",
    x     = "% of variants",
    y     = "Gene"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.y     = element_text(face = "italic"),
    legend.position = "right",
    panel.grid      = element_blank()
  )

print(p)

ggsave(file.path(BASE_DIR, "Cohorts/Relapse/gainloss_bar_relapse.pdf"),
       plot = p, width = 7, height = 6, units = "in")
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/gainloss_bar_relapse.png"),
       plot = p, width = 7, height = 6, units = "in", dpi = 300)

message("Saved: Cohorts/Relapse/gainloss_bar_relapse.pdf/.png")
