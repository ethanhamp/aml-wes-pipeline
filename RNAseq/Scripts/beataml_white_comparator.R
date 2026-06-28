library(limma)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR    <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
beataml_dir <- Sys.getenv("BEATAML_DIR", unset = file.path(BASE_DIR, "../BeatAML"))
expr_file   <- file.path(beataml_dir, "beataml_waves1to4_norm_exp_dbgap.txt")
clin_file   <- file.path(beataml_dir, "BEATAML1.0-COHORT.clinical.tsv.gz.tsv")
out_dir     <- file.path(BASE_DIR, "RNAseq/BeatAML/Output")
fig_dir     <- file.path(BASE_DIR, "RNAseq/BeatAML/Figures")

# ── Step 1: Load expression matrix ────────────────────────────────────────────
cat("== Step 1: Loading BeatAML expression ==\n")
expr_raw <- as.data.frame(data.table::fread(expr_file, sep = "\t",
                                                    check.names = FALSE,
                                                    quote = ""))
cat("  Dimensions:", nrow(expr_raw), "genes x", ncol(expr_raw) - 4, "samples\n")

# gene metadata in first 4 columns; expression in the rest
gene_meta <- expr_raw[, 1:4]
expr_mat  <- as.matrix(expr_raw[, 5:ncol(expr_raw)])
rownames(expr_mat) <- gene_meta$stable_id  # Ensembl IDs (no version)

# strip trailing 'R' from sample IDs to match clinical BA-IDs
colnames(expr_mat) <- sub("R$", "", colnames(expr_mat))
cat("  Sample ID example (after strip):", colnames(expr_mat)[1], "\n")

# ── Step 2: Load clinical metadata ────────────────────────────────────────────
cat("== Step 2: Loading clinical metadata ==\n")
clin <- read.table(gzfile(clin_file), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "")
cat("  Clinical rows:", nrow(clin), "| cols:", ncol(clin), "\n")

# keep one row per sample (some patients have multiple rows)
clin <- clin[!duplicated(clin$sample), ]
rownames(clin) <- clin$sample

# ── Step 3: Filter to white patients, Primary vs Recurrence ───────────────────
cat("== Step 3: Filtering to white patients ==\n")
shared_ids <- intersect(colnames(expr_mat), rownames(clin))
cat("  Samples with both expression + clinical:", length(shared_ids), "\n")

clin_shared <- clin[shared_ids, ]
white_ids   <- shared_ids[clin_shared$race.demographic == "white" &
                          clin_shared$tumor_descriptor.samples %in%
                            c("Primary", "Recurrence")]
clin_white  <- clin_shared[white_ids, ]
expr_white  <- expr_mat[, white_ids]

cat("  White Primary:", sum(clin_white$tumor_descriptor.samples == "Primary"), "\n")
cat("  White Recurrence:", sum(clin_white$tumor_descriptor.samples == "Recurrence"), "\n")

# ── Step 4: Limma DE — Primary vs Recurrence ──────────────────────────────────
cat("== Step 4: Limma DE ==\n")
tumor_desc <- factor(clin_white$tumor_descriptor.samples,
                     levels = c("Primary", "Recurrence"))
design <- model.matrix(~ tumor_desc)
colnames(design) <- c("Intercept", "Recurrence_vs_Primary")

fit  <- lmFit(expr_white, design)
fit  <- eBayes(fit)
res  <- topTable(fit, coef = "Recurrence_vs_Primary", number = Inf, sort.by = "none")
res$ensembl_id  <- rownames(res)
res$gene_symbol <- gene_meta$display_label[match(rownames(res), gene_meta$stable_id)]
res <- res[order(res$adj.P.Val), ]
write.csv(res, file.path(out_dir, "de_primary_vs_recurrence_white.csv"), row.names = FALSE)

sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05, ]
cat("  Significant DE genes (FDR < 0.05):", nrow(sig), "\n")
cat("  Up in Recurrence:", sum(sig$logFC > 0), "\n")
cat("  Down in Recurrence:", sum(sig$logFC < 0), "\n")

# volcano
res$color_group <- case_when(
  !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$logFC > 1  ~ "Up in Recurrence",
  !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$logFC < -1 ~ "Down in Recurrence",
  TRUE ~ "NS")
res$sig_label <- ifelse(res$color_group != "NS", res$gene_symbol, NA)

pdf(file.path(fig_dir, "volcano_primary_vs_recurrence_white.pdf"), width = 7, height = 6)
ggplot(res, aes(logFC, -log10(adj.P.Val), color = color_group)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_text(aes(label = sig_label), size = 2.5, hjust = -0.1, na.rm = TRUE,
            check_overlap = TRUE) +
  scale_color_manual(values = c("Up in Recurrence" = "#d73027",
                                "Down in Recurrence" = "#4575b4", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(x = "log2 fold change (Recurrence vs Primary)", y = "-log10 adjusted p-value",
       color = NULL,
       title = "BeatAML white: Recurrence vs Primary",
       subtitle = paste0("n=", sum(tumor_desc == "Primary"), " Primary vs n=",
                         sum(tumor_desc == "Recurrence"), " Recurrence")) +
  theme_classic(base_size = 12)
dev.off()

# ── Step 5: GSEA — Hallmarks ranked by limma t-statistic ──────────────────────
cat("== Step 5: GSEA Hallmarks ==\n")
hallmarks <- msigdbr(species = "Homo sapiens", collection = "H")
pathways  <- split(hallmarks$ensembl_gene, hallmarks$gs_name)

rank_vec <- setNames(res$t, res$ensembl_id)
rank_vec <- sort(rank_vec[!is.na(rank_vec)], decreasing = TRUE)
cat("  Rank vector range: [", round(min(rank_vec), 2), ",",
    round(max(rank_vec), 2), "]\n")

set.seed(42)
gsea_res <- fgsea(pathways = pathways, stats = rank_vec,
                  minSize = 15, maxSize = 500, nPermSimple = 10000)
gsea_res <- gsea_res[order(gsea_res$padj), ]
write.csv(as.data.frame(gsea_res[, c("pathway", "pval", "padj", "NES", "size")]),
          file.path(out_dir, "gsea_hallmarks_primary_vs_recurrence_white.csv"),
          row.names = FALSE)

sig_gsea <- gsea_res[!is.na(gsea_res$padj) & gsea_res$padj < 0.25, ]
cat("  Significant pathways (FDR < 0.25):", nrow(sig_gsea), "\n")
print(as.data.frame(sig_gsea[, c("pathway", "NES", "padj")]), row.names = FALSE)

# ── Step 6: Side-by-side comparison with Black cohort (UPenn pilot) ────────────
cat("== Step 6: Black vs white GSEA comparison ==\n")
black_gsea_file <- "RNAseq/UPenn/Output/gsea_hallmarks_denovo_vs_refractory.csv"
if (file.exists(black_gsea_file)) {
  black_gsea <- read.csv(black_gsea_file)
  colnames(black_gsea)[colnames(black_gsea) == "NES"]  <- "NES_black"
  colnames(black_gsea)[colnames(black_gsea) == "padj"] <- "padj_black"
  white_gsea <- as.data.frame(gsea_res[, c("pathway", "NES", "padj")])
  colnames(white_gsea)[colnames(white_gsea) == "NES"]  <- "NES_white"
  colnames(white_gsea)[colnames(white_gsea) == "padj"] <- "padj_white"

  comparison <- merge(black_gsea[, c("pathway", "NES_black", "padj_black")],
                      white_gsea[, c("pathway", "NES_white", "padj_white")],
                      by = "pathway", all = TRUE)
  comparison$sig_black <- !is.na(comparison$padj_black) & comparison$padj_black < 0.25
  comparison$sig_white <- !is.na(comparison$padj_white) & comparison$padj_white < 0.25
  comparison$category  <- case_when(
    comparison$sig_black & comparison$sig_white  ~ "Shared",
    comparison$sig_black & !comparison$sig_white ~ "Black-specific",
    !comparison$sig_black & comparison$sig_white ~ "White-specific",
    TRUE ~ "NS in both"
  )
  write.csv(comparison,
            file.path(out_dir, "gsea_hallmarks_black_vs_white_comparison.csv"),
            row.names = FALSE)

  cat("  Shared (sig in both):", sum(comparison$category == "Shared"), "\n")
  cat("  Black-specific:", sum(comparison$category == "Black-specific"), "\n")
  cat("  White-specific:", sum(comparison$category == "White-specific"), "\n")

  # side-by-side bar chart for pathways significant in either cohort
  plot_data <- comparison[comparison$category != "NS in both", ] %>%
    mutate(pathway_clean = gsub("HALLMARK_", "", pathway)) %>%
    pivot_longer(cols = c("NES_black", "NES_white"),
                 names_to = "cohort", values_to = "NES") %>%
    mutate(cohort = ifelse(cohort == "NES_black", "Black (UPenn pilot)", "White (BeatAML)"))

  pdf(file.path(fig_dir, "gsea_hallmarks_black_vs_white_comparison.pdf"),
      width = 10, height = 7)
  print(
    ggplot(plot_data,
           aes(x = reorder(pathway_clean, NES), y = NES, fill = cohort)) +
      geom_col(position = position_dodge(0.7), width = 0.65) +
      coord_flip() +
      scale_fill_manual(values = c("Black (UPenn pilot)" = "#2166ac",
                                   "White (BeatAML)" = "#d6604d")) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
      labs(x = NULL, y = "Normalized Enrichment Score",
           fill = NULL,
           title = "GSEA Hallmarks: Recurrence vs Primary",
           subtitle = "Pathways significant (FDR < 0.25) in at least one cohort") +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom")
  )
  dev.off()
  cat("  Comparison figure written.\n")
} else {
  cat("  Black cohort GSEA not found — skipping comparison\n")
}

cat("\n== Done ==\n")
cat("Outputs in:", out_dir, "and", fig_dir, "\n")
