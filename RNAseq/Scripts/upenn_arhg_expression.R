library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)

BASE_DIR   <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
fig_dir    <- file.path(BASE_DIR, "RNAseq/UPenn/Figures")
out_dir    <- file.path(BASE_DIR, "RNAseq/UPenn/Output")
exome_file <- file.path(BASE_DIR, "Cohorts/UPenn/Output/upenn_all_combined.xlsx")

arhg_gefs  <- c(paste0("ARHGEF", c(1:19, 25, 26, 28)))
arhg_gaps  <- paste0("ARHGAP", c(1:10, "11A", "11B", 12:45))
arhg_genes <- c(arhg_gefs, arhg_gaps)

# ── Load VST matrix and metadata ──────────────────────────────────────────────
vst_mat <- as.matrix(read.csv(file.path(out_dir, "upenn_vst_normalized.csv"),
                              row.names = 1, check.names = FALSE))

manifest    <- read.delim(file.path(BASE_DIR, "Cohorts/Manifests/master_sample_manifest.tsv"),
                          stringsAsFactors = FALSE)
upenn_man   <- manifest[!is.na(manifest$cohort_UPenn) & manifest$cohort_UPenn == TRUE, ]
upenn_man$sample_id <- as.character(upenn_man$sample_id)
upenn_man$disease_status <- upenn_man$timepoint
upenn_man$disease_status[upenn_man$disease_status == "Dx"]       <- "de novo"
upenn_man$disease_status[upenn_man$disease_status == "blast"]    <- "de novo"
upenn_man$disease_status[upenn_man$disease_status == "Relapse"]  <- "relapse"
upenn_man$disease_status[upenn_man$disease_status == "Unknown"]  <- "unknown"
status_df <- unique(upenn_man[, c("sample_id", "disease_status")])

# ── Map ARHG symbols to Ensembl ───────────────────────────────────────────────
arhg_ens <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = arhg_genes,
         column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
)
detected_ens <- intersect(arhg_ens[!is.na(arhg_ens)], rownames(vst_mat))
ens_to_sym   <- setNames(names(arhg_ens[!is.na(arhg_ens)]), arhg_ens[!is.na(arhg_ens)])

arhg_sub        <- vst_mat[detected_ens, , drop = FALSE]
rownames(arhg_sub) <- ens_to_sym[rownames(arhg_sub)]

cat("Detected ARHG genes:", nrow(arhg_sub),
    "(", sum(startsWith(rownames(arhg_sub), "ARHGEF")), "GEFs,",
    sum(startsWith(rownames(arhg_sub), "ARHGAP")), "GAPs )\n")

# ── Build long table ───────────────────────────────────────────────────────────
arhg_long <- as.data.frame(t(arhg_sub)) |>
  tibble::rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "gene", values_to = "vst") |>
  mutate(family = ifelse(startsWith(gene, "ARHGEF"), "GEF", "GAP")) |>
  left_join(status_df, by = "sample_id")

# Per-gene summary (median + IQR across all samples)
gene_summary <- arhg_long |>
  group_by(gene, family) |>
  summarise(
    median_vst = median(vst),
    q25        = quantile(vst, 0.25),
    q75        = quantile(vst, 0.75),
    .groups    = "drop"
  ) |>
  arrange(desc(median_vst))

# ── Figure 1: GAP vs GEF — per-gene median distribution ──────────────────────
# Each point = one gene; shows expression landscape of each family

wt  <- wilcox.test(
  gene_summary$median_vst[gene_summary$family == "GAP"],
  gene_summary$median_vst[gene_summary$family == "GEF"]
)
p_str <- if (wt$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", wt$p.value)

# For the y annotation, place above the max value
y_ann <- max(gene_summary$median_vst) + 0.5

p1 <- ggplot(gene_summary, aes(x = family, y = median_vst, fill = family)) +
  geom_violin(alpha = 0.35, trim = FALSE, width = 0.7, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7, color = "grey30") +
  geom_jitter(aes(color = family), width = 0.12, size = 2.5, alpha = 0.85,
              shape = 21, stroke = 0.3) +
  geom_text_repel(aes(label = gene, color = family), size = 2.4,
                  max.overlaps = 20, show.legend = FALSE,
                  segment.size = 0.3, segment.alpha = 0.5) +
  annotate("segment", x = 1, xend = 2,
           y = y_ann - 0.1, yend = y_ann - 0.1, color = "grey40") +
  annotate("text", x = 1.5, y = y_ann + 0.1,
           label = paste("Wilcoxon", p_str), size = 3.5, color = "grey30") +
  scale_fill_manual(values  = c("GAP" = "#d73027", "GEF" = "#4575b4")) +
  scale_color_manual(values = c("GAP" = "#d73027", "GEF" = "#4575b4")) +
  scale_x_discrete(labels = c(
    "GAP" = paste0("ARHGAP\n(n=", sum(gene_summary$family == "GAP"), ")"),
    "GEF" = paste0("ARHGEF\n(n=", sum(gene_summary$family == "GEF"), ")")
  )) +
  labs(x = NULL, y = "Median VST expression (across 43 samples)",
       title = "ARHGAP vs ARHGEF expression") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

pdf(file.path(fig_dir, "arhg_gap_vs_gef.pdf"), width = 6, height = 6)
print(p1)
dev.off()
cat("Figure 1 written: arhg_gap_vs_gef.pdf (Wilcoxon", p_str, ")\n")

# ── Figure 2: Top 5 expressed ARHG genes ─────────────────────────────────────
top5_genes <- head(gene_summary$gene, 5)
cat("Top 5 expressed:", paste(top5_genes, collapse = ", "), "\n")

top5_long <- arhg_long |>
  filter(gene %in% top5_genes) |>
  mutate(
    gene   = factor(gene, levels = top5_genes),
    family = ifelse(startsWith(as.character(gene), "ARHGEF"), "GEF", "GAP")
  )

# Median labels for geom_text
top5_medians <- top5_long |>
  group_by(gene) |>
  summarise(med = median(vst), .groups = "drop")

p2 <- ggplot(top5_long, aes(x = gene, y = vst, fill = family)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = "grey30") +
  geom_jitter(aes(color = family), width = 0.15, size = 1.8,
              alpha = 0.6, shape = 21, stroke = 0.3) +
  geom_text(data = top5_medians,
            aes(x = gene, y = med, label = round(med, 1)),
            inherit.aes = FALSE, vjust = -0.6, size = 3.2, fontface = "bold") +
  scale_fill_manual(values  = c("GAP" = "#d73027", "GEF" = "#4575b4"),
                    labels  = c("GAP" = "ARHGAP", "GEF" = "ARHGEF")) +
  scale_color_manual(values = c("GAP" = "#d73027", "GEF" = "#4575b4"),
                     labels = c("GAP" = "ARHGAP", "GEF" = "ARHGEF")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "VST expression", fill = "Family", color = "Family",
       title = "Top 5 expressed ARHG genes",
       subtitle = "Ranked by median VST across 43 UPenn RNA-seq samples") +
  theme_classic(base_size = 13) +
  theme(legend.title = element_text(size = 11))

pdf(file.path(fig_dir, "arhg_top5_expressed.pdf"), width = 6, height = 5)
print(p2)
dev.off()
cat("Figure 2 written: arhg_top5_expressed.pdf\n")

# ── Figure 3: Full ARHG ranked expression — lollipop ─────────────────────────
# All 43 detected genes ranked by median VST, colored by family
gene_summary_ranked <- gene_summary |>
  mutate(gene = factor(gene, levels = gene[order(median_vst)]))

p3 <- ggplot(gene_summary_ranked,
             aes(x = gene, y = median_vst, color = family)) +
  geom_segment(aes(xend = gene, y = 0, yend = median_vst),
               linewidth = 0.7, alpha = 0.6) +
  geom_point(size = 3, alpha = 0.9) +
  geom_errorbar(aes(ymin = q25, ymax = q75),
                width = 0.3, linewidth = 0.5, alpha = 0.5) +
  coord_flip() +
  scale_color_manual(values = c("GAP" = "#d73027", "GEF" = "#4575b4"),
                     labels = c("GAP" = "ARHGAP", "GEF" = "ARHGEF")) +
  labs(x = NULL, y = "Median VST expression (IQR bars)",
       color = "Family",
       title = "ARHG gene expression - UPenn RNA-seq",
       subtitle = "All 43 detected genes ranked by median VST") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 8))

pdf(file.path(fig_dir, "arhg_all_ranked_expression.pdf"), width = 6, height = 9)
print(p3)
dev.off()
cat("Figure 3 written: arhg_all_ranked_expression.pdf\n")

# ── Figure 4: GAP vs GEF per-sample median (paired within sample) ─────────────
# Each sample has one GAP median and one GEF median — shows systematic difference
per_sample_family <- arhg_long |>
  filter(!is.na(disease_status), disease_status != "unknown") |>
  group_by(sample_id, family, disease_status) |>
  summarise(median_vst = median(vst), .groups = "drop") |>
  mutate(disease_status = factor(disease_status,
                                  levels = c("de novo", "refractory", "relapse")))

wt2   <- wilcox.test(
  per_sample_family$median_vst[per_sample_family$family == "GAP"],
  per_sample_family$median_vst[per_sample_family$family == "GEF"],
  paired = FALSE
)
p_str2 <- if (wt2$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", wt2$p.value)

p4 <- ggplot(per_sample_family,
             aes(x = family, y = median_vst, fill = family)) +
  geom_violin(alpha = 0.35, trim = FALSE, color = NA, width = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7, color = "grey30") +
  geom_jitter(aes(color = disease_status), width = 0.1, size = 2.5,
              alpha = 0.85, shape = 21, stroke = 0.4) +
  scale_fill_manual(values  = c("GAP" = "#d73027", "GEF" = "#4575b4")) +
  scale_color_manual(values = c("de novo" = "#4575b4", "refractory" = "#d73027",
                                "relapse" = "#f46d43"),
                     name = "Disease status") +
  scale_x_discrete(labels = c("GAP" = "ARHGAP", "GEF" = "ARHGEF")) +
  annotate("text", x = 1.5,
           y = max(per_sample_family$median_vst) + 0.15,
           label = paste("Wilcoxon", p_str2), size = 3.5, color = "grey30") +
  labs(x = NULL, y = "Per-sample median VST across family members",
       title = "GAP vs GEF expression per sample",
       subtitle = "Each point = one sample's median across all detected family members") +
  theme_classic(base_size = 13) +
  guides(fill = "none")

pdf(file.path(fig_dir, "arhg_gap_vs_gef_per_sample.pdf"), width = 6, height = 5.5)
print(p4)
dev.off()
cat("Figure 4 written: arhg_gap_vs_gef_per_sample.pdf\n")

cat("\n== Summary table: all detected ARHG genes ==\n")
print(gene_summary, n = Inf)

# ── Figure 5: ARHG expression with gene-specific carrier highlighting ──────────
# For each gene on the x-axis, the sample(s) carrying a mutation IN THAT GENE
# are highlighted — so ARHGAP20 column highlights sample 4219, ARHGEF17 highlights
# 8476, etc. This lets you ask: does the carrier have outlier expression of their
# mutated gene?

exome     <- suppressWarnings(read_excel(exome_file))
exome_som <- exome[exome$Call_Type == "somatic", ]

# One row per sample-gene mutation pair (a sample may have multiple ARHG mutations)
gene_carrier_map <- exome_som[exome_som$Gene %in% arhg_genes,
                               c("Sample", "Gene")] |>
  mutate(sample_id = as.character(Sample), gene = Gene) |>
  dplyr::select(sample_id, gene) |>
  distinct() |>
  mutate(is_carrier = TRUE)

# Join carrier flag onto arhg_long (gene-specific: only flags the matching column)
arhg_long_carriers <- arhg_long |>
  left_join(gene_carrier_map, by = c("sample_id", "gene")) |>
  mutate(is_carrier = replace(is_carrier, is.na(is_carrier), FALSE))

carriers_only <- filter(arhg_long_carriers, is_carrier)

status_colors <- c("de novo"     = "#4575b4",
                   "refractory"  = "#d73027",
                   "relapse"     = "#f46d43",
                   "unknown"     = "grey60")

# Gene order: same ranked order as Figure 3
gene_order <- levels(gene_summary_ranked$gene)
present_order <- gene_order[gene_order %in% unique(arhg_long_carriers$gene)]

arhg_long_carriers <- arhg_long_carriers |>
  filter(!is.na(disease_status), disease_status != "relapse") |>
  mutate(gene = factor(gene, levels = present_order))

carriers_only <- filter(arhg_long_carriers, is_carrier)

p5 <- ggplot(arhg_long_carriers,
             aes(x = gene, y = vst, color = disease_status)) +
  # Background: all samples, small and faint
  geom_jitter(data = filter(arhg_long_carriers, !is_carrier),
              width = 0.2, size = 1.2, alpha = 0.4, shape = 16) +
  # Crossbar at median per gene × status
  stat_summary(fun = median, geom = "crossbar",
               width = 0.45, linewidth = 0.45,
               position = position_dodge(0.4)) +
  # Carriers: filled circle outlined in black, larger
  geom_point(data = carriers_only,
             aes(fill = disease_status),
             shape = 21, size = 3.5, stroke = 1.1,
             color = "black") +
  # Label each carrier with its sample ID
  geom_text_repel(data = carriers_only,
                  aes(label = sample_id),
                  size = 2.6, fontface = "bold",
                  color = "black",
                  box.padding = 0.4,
                  segment.size = 0.3,
                  show.legend = FALSE) +
  scale_color_manual(values = status_colors) +
  scale_fill_manual(values  = status_colors) +
  scale_x_discrete(drop = TRUE) +
  labs(x = NULL, y = "VST expression",
       color = "Status", fill = "Status",
       title = "ARHG expression — UPenn cohort (carriers highlighted)",
       subtitle = "Outlined points: sample carries mutation in that specific gene") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  guides(fill = "none")

pdf(file.path(fig_dir, "arhg_expression_upenn_carriers.pdf"), width = 16, height = 6)
print(p5)
dev.off()
cat("Figure 5 written: arhg_expression_upenn_carriers.pdf\n")
cat("  Carrier annotations:\n")
print(gene_carrier_map)
