library(tidyverse)
library(readxl)
library(pheatmap)
library(ggrepel)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
fig_dir  <- file.path(BASE_DIR, "RNAseq/Figures")

# ── ARHG gene list ────────────────────────────────────────────────────────────
arhg_gefs <- c(paste0("ARHGEF", c(1:19, 25, 26, 28, 37, 38, 39, 40)), "ARHGEF10L")
arhg_gaps <- paste0("ARHGAP", c(1:10, "11A", "11B", 12:45))
arhg_all  <- c(arhg_gefs, arhg_gaps)

# ── Load NOISeq expression matrix (RDS cache for fast reloads) ───────────────
noiseq_rds <- file.path(BASE_DIR, "RNAseq/NOISeqAdjusted1429_Sept2025.rds")
if (file.exists(noiseq_rds)) {
  cat("Loading NOISeq matrix from RDS cache...\n")
  expr_raw <- readRDS(noiseq_rds)
} else {
  cat("Loading NOISeq matrix from text (first time — will cache as RDS)...\n")
  expr_raw <- read.delim(file.path(BASE_DIR, "RNAseq/NOISeqAdjusted1429_Sept2025.txt"),
                         row.names = 1, check.names = FALSE)
  gene_sym           <- sapply(strsplit(rownames(expr_raw), "\\|"), `[`, 3)
  rownames(expr_raw) <- make.unique(gene_sym)
  saveRDS(expr_raw, noiseq_rds)
  cat("Cached to", noiseq_rds, "\n")
}
cat("Matrix:", nrow(expr_raw), "genes x", ncol(expr_raw), "samples\n")

arhg_present <- intersect(arhg_all, rownames(expr_raw))
cat("ARHG genes detected:", length(arhg_present),
    "(", sum(startsWith(arhg_present, "ARHGEF")), "GEFs,",
    sum(startsWith(arhg_present, "ARHGAP")), "GAPs)\n")

expr_arhg <- expr_raw[arhg_present, , drop = FALSE]
expr_log  <- log2(as.matrix(expr_arhg) + 0.01)

# ── Build carrier map from curated ARHG variant catalog ───────────────────────
# Source: Variants/arhg_variants_clean.xlsx — the authoritative list of confirmed
# ARHG somatic mutations across all patients, including those whose per-case exome
# files are not yet in Cohorts/ARHG/Datasets/. Using this as carrier source for
# expression figures so mutation status (not exome pipeline filters) determines
# who gets labeled.
variants_catalog <- read_excel(file.path(BASE_DIR, "Variants/arhg_variants_clean.xlsx")) %>%
  mutate(patient_id = sub("_blast$|_blasts$|-BT$", "", Sample_ID)) %>%
  select(patient_id, gene = Hugo_Symbol) %>%
  distinct()

cat("ARHG carrier-gene pairs (catalog):", nrow(variants_catalog),
    "| Unique patients:", n_distinct(variants_catalog$patient_id), "\n")

# Map patient IDs to NOISeq UIDs via RNA-seq manifest
rna_manifest <- read_excel(file.path(BASE_DIR, "RNAseq/RNAseq samples 06122025.xlsx")) %>%
  filter(!is.na(access)) %>%
  select(uid, access)

carrier_uid <- inner_join(variants_catalog, rna_manifest,
                          by = c("patient_id" = "access"),
                          relationship = "many-to-many") %>%
  filter(uid %in% colnames(expr_log)) %>%
  select(uid, patient_id, gene) %>%
  distinct() %>%
  # One UID per patient per gene — avoids duplicate labels when manifest has
  # multiple UIDs for the same access ID
  group_by(patient_id, gene) %>%
  slice(1) %>%
  ungroup()

cat("Carriers with NOISeq data:", n_distinct(carrier_uid$patient_id),
    "| UID-gene pairs:", nrow(carrier_uid), "\n")

# ── Long expression table ─────────────────────────────────────────────────────
expr_long <- as.data.frame(t(expr_log)) %>%
  rownames_to_column("uid") %>%
  pivot_longer(-uid, names_to = "gene", values_to = "log2expr") %>%
  mutate(family = ifelse(startsWith(gene, "ARHGEF"), "GEF", "GAP")) %>%
  left_join(carrier_uid %>% mutate(is_carrier = TRUE),
            by = c("uid", "gene")) %>%
  mutate(
    is_carrier = replace_na(is_carrier, FALSE),
    patient_id = replace_na(patient_id, "")
  )

fam_colors <- c("GEF" = "#4575b4", "GAP" = "#d73027")

# Gene order: ascending by median expression (matches UPenn figure orientation)
gene_order <- expr_long %>%
  group_by(gene) %>%
  summarise(med = median(log2expr), .groups = "drop") %>%
  arrange(med) %>%
  pull(gene)

# ── Figure 1: Carrier strip plot — all 1429 samples, carriers highlighted ─────
# Replicates the UPenn carrier figure (upenn_arhg_expression.R Figure 5) using
# the Eisfeld NOISeq 1429 cohort.  Each column = one ARHG gene, sorted left→right
# by ascending median log2 expression.  Small faint points = all 1429 samples.
# Large outlined circles = patients carrying a somatic mutation in that gene.

plot1 <- expr_long %>%
  mutate(gene = factor(gene, levels = gene_order))

carriers1 <- filter(plot1, is_carrier)

p1 <- ggplot(plot1, aes(x = gene, y = log2expr, color = family)) +
  geom_jitter(data = filter(plot1, !is_carrier),
              width = 0.2, size = 0.6, alpha = 0.2, shape = 16) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, linewidth = 0.45,
               position = position_dodge(0)) +
  geom_point(data = carriers1,
             aes(fill = family), shape = 21, size = 3.5,
             stroke = 1.1, color = "black") +
  geom_text_repel(data = carriers1,
                  aes(label = patient_id), size = 2.4,
                  fontface = "bold", color = "black",
                  box.padding = 0.45, segment.size = 0.3,
                  show.legend = FALSE) +
  scale_color_manual(values = fam_colors,
                     labels = c("GEF" = "ARHGEF", "GAP" = "ARHGAP")) +
  scale_fill_manual(values  = fam_colors,
                    labels  = c("GEF" = "ARHGEF", "GAP" = "ARHGAP")) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = NULL, y = "log2(NOISeq + 0.01)",
       color = "Family", fill = "Family",
       title = "ARHG expression — Eisfeld 1429 cohort (carriers highlighted)",
       subtitle = sprintf(
         "Outlined points: patient carries somatic mutation in that gene  |  %d carriers with RNA-seq (of 25 ARHG-mutant patients)",
         n_distinct(carrier_uid$uid)
       )) +
  theme_classic(base_size = 9) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "top") +
  guides(fill = "none")

pdf(file.path(fig_dir, "arhg_expression_carriers_1429.pdf"), width = 18, height = 6)
print(p1)
dev.off()
cat("Figure 1 saved: arhg_expression_carriers_1429.pdf\n")

# ── Figure 2: Heatmap — carriers + 200 random WT ─────────────────────────────
set.seed(42)
wt_uids   <- setdiff(colnames(expr_log), unique(carrier_uid$uid))
show_uids <- c(unique(carrier_uid$uid), sample(wt_uids, min(200, length(wt_uids))))
heat_mat  <- expr_log[arhg_present, show_uids]

col_ann <- data.frame(
  ARHG_mutant = ifelse(colnames(heat_mat) %in% unique(carrier_uid$uid), "Mutant", "WT"),
  row.names   = colnames(heat_mat)
)
row_ann <- data.frame(
  Family    = ifelse(startsWith(rownames(heat_mat), "ARHGEF"), "GEF", "GAP"),
  row.names = rownames(heat_mat)
)
ann_colors <- list(
  ARHG_mutant = c(Mutant = "#d73027", WT = "grey85"),
  Family      = c(GEF = "#4575b4", GAP = "#d73027")
)

pdf(file.path(fig_dir, "arhg_expression_heatmap_1429.pdf"), width = 14, height = 10)
pheatmap(
  heat_mat,
  annotation_col    = col_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  color             = colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(100),
  fontsize_row      = 7,
  main              = sprintf(
    "ARHG expression — %d carriers + 200 random WT (log2 NOISeq)",
    n_distinct(carrier_uid$uid)
  )
)
dev.off()
cat("Figure 2 saved: arhg_expression_heatmap_1429.pdf\n")

# ── Figure 3: Gene-specific carrier stripplot ─────────────────────────────────
# Like Figure 1 but only shows genes where at least one carrier exists,
# so the carrier points are easier to see.
carrier_genes <- intersect(unique(carrier_uid$gene), arhg_present)
strip_data    <- filter(expr_long, gene %in% carrier_genes) %>%
  mutate(gene = factor(gene, levels = intersect(gene_order, carrier_genes)))

p3 <- ggplot(strip_data, aes(x = gene, y = log2expr, color = family)) +
  geom_jitter(data = filter(strip_data, !is_carrier),
              width = 0.25, size = 0.5, alpha = 0.2, shape = 16) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, linewidth = 0.5) +
  geom_point(data = filter(strip_data, is_carrier),
             aes(fill = family), shape = 21, size = 4.5,
             stroke = 1.1, color = "black") +
  geom_text_repel(
    data        = filter(strip_data, is_carrier),
    aes(label   = patient_id),
    size        = 2.6, fontface = "bold", color = "black",
    box.padding = 0.5, segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = fam_colors,
                     labels = c("GEF" = "ARHGEF", "GAP" = "ARHGAP")) +
  scale_fill_manual(values  = fam_colors,
                    labels  = c("GEF" = "ARHGEF", "GAP" = "ARHGAP")) +
  labs(x = NULL, y = "log2(NOISeq + 0.01)",
       color = "Family", fill = "Family",
       title = "ARHG carrier expression — mutated gene only",
       subtitle = "Filled circle: patient carries somatic mutation in that gene  |  background: all 1429 patients") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  guides(fill = "none")

fig3_w <- max(6, length(carrier_genes) * 1.0)
pdf(file.path(fig_dir, "arhg_carrier_stripplot_1429.pdf"), width = fig3_w, height = 6)
print(p3)
dev.off()
cat("Figure 3 saved: arhg_carrier_stripplot_1429.pdf\n")

cat("\nDone. Figures in", fig_dir, "\n")
