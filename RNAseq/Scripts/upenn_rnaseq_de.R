library(DESeq2)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(ggrepel)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR    <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
counts_dir  <- file.path(BASE_DIR, "Cohorts/UPenn/RNASeq/Counts")
exome_file  <- file.path(BASE_DIR, "Cohorts/UPenn/Output/upenn_all_combined.xlsx")
seq_summary <- file.path(BASE_DIR, "Cohorts/UPenn/upenn_seq_summary.xlsx")
out_dir     <- file.path(BASE_DIR, "RNAseq/UPenn/Output")
fig_dir     <- file.path(BASE_DIR, "RNAseq/UPenn/Figures")

arhg_gefs <- c(paste0("ARHGEF", c(1:19, 25, 26, 28)))
arhg_gaps <- paste0("ARHGAP", c(1:10, "11A", "11B", 12:45))
arhg_genes <- c(arhg_gefs, arhg_gaps)

# ── Step 1: Merge per-sample count files ──────────────────────────────────────
cat("== Step 1: Merging count files ==\n")
files    <- list.files(counts_dir, pattern = "\\.txt$", full.names = TRUE)
first    <- read.table(files[1], header = TRUE, sep = "\t", check.names = FALSE)
gene_ids <- first[[1]]

counts_list <- lapply(files, function(f) {
  d <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  stopifnot(identical(d[[1]], gene_ids))
  d[[2]]
})
counts_mat <- do.call(cbind, counts_list)
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- basename(files)
storage.mode(counts_mat) <- "integer"

rownames(counts_mat) <- sub("\\.[0-9]+$", "", rownames(counts_mat))
colnames(counts_mat) <- gsub("-RNA.*\\.txt$", "", colnames(counts_mat))
cat("  Ensembl version stripping confirmed (e.g.,", head(rownames(counts_mat), 1), ")\n")
cat("  Genes:", nrow(counts_mat), "| Samples:", ncol(counts_mat), "\n")

# ── Step 2: Load metadata ──────────────────────────────────────────────────────
cat("== Step 2: Loading metadata ==\n")

# Primary source: master manifest (more complete than seq summary)
manifest <- read.delim(file.path(BASE_DIR, "Cohorts/Manifests/master_sample_manifest.tsv"),
                       stringsAsFactors = FALSE)
upenn_manifest <- manifest[!is.na(manifest$cohort_UPenn) & manifest$cohort_UPenn == TRUE, ]
upenn_manifest$sample_id <- as.character(upenn_manifest$sample_id)

# Normalize timepoint labels to match seq summary conventions
upenn_manifest$disease_status <- upenn_manifest$timepoint
upenn_manifest$disease_status[upenn_manifest$disease_status == "Dx"]       <- "de novo"
upenn_manifest$disease_status[upenn_manifest$disease_status == "blast"]    <- "de novo"
upenn_manifest$disease_status[upenn_manifest$disease_status == "Relapse"]  <- "relapse"
status_df <- upenn_manifest[, c("sample_id", "disease_status")]
status_df <- status_df[!duplicated(status_df$sample_id), ]

# Fallback: seq summary for the 3 samples not in manifest (127, 8740, 9044)
ss <- suppressMessages(read_excel(seq_summary))
ss_df <- rbind(
  data.frame(sample_id = as.character(as.integer(ss[["sample.id...2"]])),
             disease_status = ss[["status...3"]], stringsAsFactors = FALSE),
  data.frame(sample_id = as.character(as.integer(ss[["sample.id...11"]])),
             disease_status = ss[["status...12"]], stringsAsFactors = FALSE)
)
ss_df <- ss_df[!is.na(ss_df$sample_id) & !is.na(ss_df$disease_status), ]
ss_df <- ss_df[!duplicated(ss_df$sample_id), ]
missing_from_manifest <- setdiff(colnames(counts_mat), status_df$sample_id)
status_df <- rbind(status_df,
                   ss_df[ss_df$sample_id %in% missing_from_manifest, ])
cat("  Samples resolved from manifest:",
    sum(colnames(counts_mat) %in% upenn_manifest$sample_id), "\n")
cat("  Samples resolved from seq summary fallback:",
    sum(colnames(counts_mat) %in% ss_df$sample_id &
        !colnames(counts_mat) %in% upenn_manifest$sample_id), "\n")

exome     <- suppressWarnings(read_excel(exome_file))
exome_som <- exome[exome$Call_Type == "somatic", ]
arhg_mut_ids    <- unique(exome_som$Sample[exome_som$Gene %in% arhg_genes])
all_exome_ids   <- unique(exome$Sample)

meta <- data.frame(sample_id = colnames(counts_mat), stringsAsFactors = FALSE)
meta <- merge(meta, status_df, by = "sample_id", all.x = TRUE)
meta$has_exome   <- meta$sample_id %in% as.character(all_exome_ids)
meta$arhg_mutant <- meta$sample_id %in% as.character(arhg_mut_ids)
meta$library_size <- colSums(counts_mat)[meta$sample_id]
rownames(meta) <- meta$sample_id
meta <- meta[colnames(counts_mat), ]

# Samples explicitly excluded from disease-status analyses
unknown_ids  <- meta$sample_id[is.na(meta$disease_status)]   # 127
relapse_id   <- meta$sample_id[!is.na(meta$disease_status) & meta$disease_status == "relapse"]  # 6714
# Duplicate-patient exclusions: one sample per patient retained
# R2928: keep 6286 (ARHG-mutant), drop 5838
# R2127: keep 4313, drop 4393
duplicate_drops <- c("5838", "4393")

cat("  Unknown status (excluded from DE):", paste(unknown_ids, collapse = ", "), "\n")
cat("  Relapse sample (excluded from primary DE):", relapse_id, "\n")
cat("  Duplicate-patient drops:", paste(duplicate_drops, collapse = ", "), "\n")
cat("  Disease status (all 43 samples):\n")
print(table(meta$disease_status, useNA = "ifany"))
cat("  ARHG-mutant (has exome):", sum(meta$arhg_mutant), "\n")

# ── Step 3: DESeq2 — all 43 samples, VST for PCA/normalization ────────────────
cat("== Step 3: DESeq2 VST normalization (all samples) ==\n")

# LOW-COUNT FILTER: ≥10 counts in ≥10 samples (calibrated to smallest group n=13)
keep <- rowSums(counts_mat >= 10) >= 10
cat("  Genes after low-count filter (≥10 counts in ≥10 samples):", sum(keep), "\n")
counts_filt <- counts_mat[keep, ]

dds_all <- DESeqDataSetFromMatrix(countData = counts_filt, colData = meta, design = ~ 1)
dds_all <- estimateSizeFactors(dds_all)
vst_all  <- vst(dds_all, blind = TRUE)
vst_mat  <- assay(vst_all)
write.csv(as.data.frame(vst_mat), file.path(out_dir, "upenn_vst_normalized.csv"))

# ── Step 4: QC ────────────────────────────────────────────────────────────────
cat("== Step 4: QC figures ==\n")

# library sizes
lib_df <- data.frame(
  sample_id      = colnames(counts_mat),
  library_size   = colSums(counts_mat),
  disease_status = meta$disease_status,
  arhg_mutant    = meta$arhg_mutant
)
pdf(file.path(fig_dir, "library_sizes.pdf"), width = 12, height = 5)
ggplot(lib_df, aes(x = reorder(sample_id, library_size), y = library_size / 1e6,
                   fill = disease_status)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Library size (M reads)", fill = "Status",
       title = "UPenn RNA-seq library sizes") +
  theme_classic(base_size = 11)
dev.off()

# PCA — check for technical + biological structure
pca_data <- plotPCA(vst_all, intgroup = c("disease_status", "arhg_mutant"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
pca_data$label        <- pca_data$name
pca_data$library_size <- meta[pca_data$name, "library_size"]

cat("  PCA variance explained: PC1 =", pct_var[1], "%, PC2 =", pct_var[2], "%\n")
cat("  Spearman r (PC1 vs library size):",
    round(cor(pca_data$PC1, pca_data$library_size, method = "spearman"), 3), "\n")
cat("  Spearman r (PC2 vs library size):",
    round(cor(pca_data$PC2, pca_data$library_size, method = "spearman"), 3), "\n")

pdf(file.path(fig_dir, "pca_disease_status.pdf"), width = 7.5, height = 5.5)
ggplot(pca_data, aes(PC1, PC2, color = disease_status, shape = arhg_mutant)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "WT", "TRUE" = "ARHG-mutant")) +
  labs(x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)"),
       color = "Status", shape = "ARHG mutation",
       title = "UPenn RNA-seq PCA") +
  theme_classic(base_size = 12)
dev.off()

pdf(file.path(fig_dir, "pca_library_size.pdf"), width = 7, height = 5.5)
ggplot(pca_data, aes(PC1, PC2, color = library_size)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
  scale_color_viridis_c(labels = scales::comma) +
  labs(x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)"),
       color = "Library size", title = "UPenn RNA-seq PCA colored by library size") +
  theme_classic(base_size = 12)
dev.off()

# ── Step 5: DE — de novo vs refractory (relapse + unknowns excluded) ──────────
cat("== Step 5: DE — de novo vs refractory ==\n")
cat("  Excluded:", paste(c(unknown_ids, relapse_id), collapse = ", "),
    "(unknown status or relapse)\n")

excl       <- c(unknown_ids, relapse_id, duplicate_drops)
meta_de1   <- meta[!meta$sample_id %in% excl &
                   meta$disease_status %in% c("de novo", "refractory"), ]
counts_de1 <- counts_filt[, rownames(meta_de1)]
meta_de1$disease_status <- factor(meta_de1$disease_status,
                                   levels = c("de novo", "refractory"))
cat("  n de novo:", sum(meta_de1$disease_status == "de novo"),
    "| n refractory:", sum(meta_de1$disease_status == "refractory"), "\n")

dds_de1 <- DESeqDataSetFromMatrix(countData = counts_de1, colData = meta_de1,
                                   design = ~ disease_status)
dds_de1 <- DESeq(dds_de1)
res_de1  <- lfcShrink(dds_de1, coef = "disease_status_refractory_vs_de.novo",
                       type = "apeglm")
# keep unshrunken stat for GSEA ranking
res_de1_raw <- results(dds_de1, name = "disease_status_refractory_vs_de.novo")

sym_map <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = rownames(res_de1),
         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
)
res_de1_df <- as.data.frame(res_de1)
res_de1_df$ensembl_id  <- rownames(res_de1_df)
res_de1_df$gene_symbol <- sym_map[rownames(res_de1_df)]
res_de1_df$stat        <- res_de1_raw$stat[rownames(res_de1_df)]
res_de1_df <- res_de1_df[order(res_de1_df$padj, na.last = TRUE), ]
write.csv(res_de1_df, file.path(out_dir, "de_denovo_vs_refractory.csv"), row.names = FALSE)

sig_de1 <- res_de1_df[!is.na(res_de1_df$padj) & res_de1_df$padj < 0.05, ]
cat("  Significant DE genes (FDR < 0.05):", nrow(sig_de1), "\n")
cat("  Up in refractory:", sum(sig_de1$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Down in refractory:", sum(sig_de1$log2FoldChange < 0, na.rm = TRUE), "\n")

# volcano
res_de1_df$sig_label <- ifelse(
  !is.na(res_de1_df$padj) & res_de1_df$padj < 0.05 & abs(res_de1_df$log2FoldChange) > 1,
  res_de1_df$gene_symbol, NA)
res_de1_df$color_group <- case_when(
  !is.na(res_de1_df$padj) & res_de1_df$padj < 0.05 & res_de1_df$log2FoldChange > 1  ~ "Up in refractory",
  !is.na(res_de1_df$padj) & res_de1_df$padj < 0.05 & res_de1_df$log2FoldChange < -1 ~ "Down in refractory",
  TRUE ~ "NS")
pdf(file.path(fig_dir, "volcano_denovo_vs_refractory.pdf"), width = 7, height = 6)
ggplot(res_de1_df, aes(log2FoldChange, -log10(padj), color = color_group)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_text_repel(aes(label = sig_label), size = 2.8, max.overlaps = 20,
                  show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = c("Up in refractory" = "#d73027",
                                "Down in refractory" = "#4575b4", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(x = "log2 fold change (refractory vs de novo)",
       y = "-log10 adjusted p-value", color = NULL,
       title = "DE: refractory vs de novo",
       subtitle = paste0("n=28 de novo vs n=12 refractory; relapse + unknowns excluded")) +
  theme_classic(base_size = 12)
dev.off()

# ── Step 5b: Outlier check — top 10 DE hits ───────────────────────────────────
cat("== Step 5b: Outlier boxplots for top DE hits ==\n")
top_hits <- head(res_de1_df[!is.na(res_de1_df$padj) & res_de1_df$padj < 0.05, ], 10)
top_ids  <- top_hits$ensembl_id
top_syms <- ifelse(!is.na(top_hits$gene_symbol), top_hits$gene_symbol, top_hits$ensembl_id)

norm_counts <- counts(dds_de1, normalized = TRUE)
outlier_df <- as.data.frame(t(norm_counts[top_ids, , drop = FALSE]))
colnames(outlier_df) <- top_syms
outlier_df$sample_id      <- rownames(outlier_df)
outlier_df$disease_status <- meta_de1[rownames(outlier_df), "disease_status"]
outlier_long <- pivot_longer(outlier_df, cols = all_of(top_syms),
                             names_to = "gene", values_to = "norm_count")

pdf(file.path(fig_dir, "top10_hits_outlier_check.pdf"), width = 14, height = 5)
ggplot(outlier_long, aes(x = disease_status, y = norm_count + 1, color = disease_status)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  labs(x = NULL, y = "Normalized counts + 1 (log10)", color = NULL,
       title = "Top 10 DE hits: outlier check") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
cat("  Outlier boxplots written.\n")

# ── Step 5c: Sensitivity — rerun including relapse in refractory group ─────────
cat("== Step 5c: Sensitivity — relapse merged into refractory ==\n")
meta_sens <- meta[!meta$sample_id %in% c(unknown_ids, duplicate_drops) &
                  meta$disease_status %in% c("de novo", "refractory", "relapse"), ]
meta_sens$disease_status_sens <- ifelse(
  meta_sens$disease_status %in% c("refractory", "relapse"), "refractory", "de novo")
meta_sens$disease_status_sens <- factor(meta_sens$disease_status_sens,
                                         levels = c("de novo", "refractory"))
counts_sens <- counts_filt[, rownames(meta_sens)]
dds_sens <- DESeqDataSetFromMatrix(countData = counts_sens, colData = meta_sens,
                                    design = ~ disease_status_sens)
dds_sens    <- DESeq(dds_sens)
res_sens    <- results(dds_sens, name = "disease_status_sens_refractory_vs_de.novo")
sig_sens    <- rownames(res_sens)[!is.na(res_sens$padj) & res_sens$padj < 0.05]
sig_primary <- sig_de1$ensembl_id
overlap_pct <- round(100 * length(intersect(sig_primary, sig_sens)) / length(sig_primary), 1)
cat("  Primary (relapse excluded): n =", nrow(sig_de1), "sig genes\n")
cat("  Sensitivity (relapse in refractory): n =", length(sig_sens), "sig genes\n")
cat("  Overlap:", length(intersect(sig_primary, sig_sens)),
    "genes (", overlap_pct, "% of primary)\n")

# ── Step 6: GSEA — Hallmarks using signed Wald statistic ──────────────────────
cat("== Step 6: GSEA Hallmarks (ranked by DESeq2 Wald stat) ==\n")

hallmarks <- msigdbr(species = "Homo sapiens", collection = "H")
pathways  <- split(hallmarks$ensembl_gene, hallmarks$gs_name)

# pull Wald stats directly from unshrunken DESeq2 results — avoids row-order issues
rank_raw <- as.data.frame(res_de1_raw)
rank_raw <- rank_raw[!is.na(rank_raw$stat), ]
rank_vec <- setNames(rank_raw$stat, rownames(rank_raw))
cat("  GSEA rank vector range: [", round(min(rank_vec), 2), ",",
    round(max(rank_vec), 2), "] — should span negative to positive\n")

set.seed(42)
gsea_res <- fgsea(pathways = pathways, stats = rank_vec,
                  minSize = 15, maxSize = 500, nPermSimple = 10000)
gsea_res <- gsea_res[order(gsea_res$padj), ]
write.csv(as.data.frame(gsea_res[, c("pathway", "pval", "padj", "NES", "size")]),
          file.path(out_dir, "gsea_hallmarks_denovo_vs_refractory.csv"), row.names = FALSE)

sig_gsea <- gsea_res[!is.na(gsea_res$padj) & gsea_res$padj < 0.25, ]
cat("  Significant pathways (FDR < 0.25):", nrow(sig_gsea), "\n")
cat("  Note: FDR < 0.25 => ~", round(0.25 * nrow(sig_gsea)),
    "expected false discoveries among significant pathways\n")
if (nrow(sig_gsea) > 0) {
  print(sig_gsea[, c("pathway", "NES", "padj")], row.names = FALSE)
}

if (nrow(sig_gsea) > 0) {
  top_paths <- rbind(
    head(sig_gsea[sig_gsea$NES > 0, ], 10),
    head(sig_gsea[sig_gsea$NES < 0, ], 10)
  )
  top_paths$pathway_clean <- gsub("HALLMARK_", "", top_paths$pathway)
  top_paths$direction <- ifelse(top_paths$NES > 0, "Up in refractory", "Down in refractory")
  pdf(file.path(fig_dir, "gsea_hallmarks_bar.pdf"), width = 8, height = 6)
  print(
    ggplot(top_paths, aes(x = reorder(pathway_clean, NES), y = NES, fill = direction)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c("Up in refractory" = "#d73027",
                                   "Down in refractory" = "#4575b4")) +
      labs(x = NULL, y = "Normalized Enrichment Score", fill = NULL,
           title = "GSEA Hallmarks: refractory vs de novo",
           subtitle = "FDR < 0.25; ranked by DESeq2 Wald statistic") +
      theme_classic(base_size = 11)
  )
  dev.off()
}

# ── Step 7: ARHG-mutant vs WT ─────────────────────────────────────────────────
cat("== Step 7: ARHG-mutant vs WT ==\n")
cat("  CONFOUND WARNING: All 5 ARHG-mutant samples are de novo; WT group is mixed.\n")
cat("  Disease-status covariate does not resolve this collinearity.\n")
cat("  Results are retained as a hypothesis catalog ONLY.\n")
cat("  No gene-level or pathway-level causal claims are valid from this comparison.\n")

meta_de2 <- meta[meta$has_exome &
                 !meta$sample_id %in% c(unknown_ids, relapse_id, duplicate_drops) &
                 meta$disease_status %in% c("de novo", "refractory"), ]
counts_de2 <- counts_filt[, rownames(meta_de2)]
meta_de2$arhg_group <- factor(ifelse(meta_de2$arhg_mutant, "ARHG_mutant", "WT"),
                               levels = c("WT", "ARHG_mutant"))
cat("  n ARHG-mutant:", sum(meta_de2$arhg_group == "ARHG_mutant"),
    "| n WT:", sum(meta_de2$arhg_group == "WT"), "\n")

dds_de2 <- DESeqDataSetFromMatrix(countData = counts_de2, colData = meta_de2,
                                   design = ~ disease_status + arhg_group)
dds_de2 <- DESeq(dds_de2)
res_de2  <- lfcShrink(dds_de2, coef = "arhg_group_ARHG_mutant_vs_WT", type = "apeglm")
res_de2_df <- as.data.frame(res_de2)
res_de2_df$ensembl_id  <- rownames(res_de2_df)
res_de2_df$gene_symbol <- sym_map[rownames(res_de2_df)]
res_de2_df <- res_de2_df[order(res_de2_df$padj, na.last = TRUE), ]
write.csv(res_de2_df, file.path(out_dir, "de_arhg_mutant_vs_wt.csv"), row.names = FALSE)

sig_de2 <- res_de2_df[!is.na(res_de2_df$padj) & res_de2_df$padj < 0.05, ]
cat("  Nominal sig genes (FDR < 0.05):", nrow(sig_de2),
    "-- DO NOT INTERPRET (confounded)\n")

res_de2_df$color_group <- case_when(
  !is.na(res_de2_df$padj) & res_de2_df$padj < 0.05 & res_de2_df$log2FoldChange > 1  ~ "Up",
  !is.na(res_de2_df$padj) & res_de2_df$padj < 0.05 & res_de2_df$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS")
res_de2_df$sig_label <- ifelse(res_de2_df$color_group != "NS", res_de2_df$gene_symbol, NA)
pdf(file.path(fig_dir, "volcano_arhg_mutant_vs_wt.pdf"), width = 7, height = 6)
ggplot(res_de2_df, aes(log2FoldChange, -log10(padj), color = color_group)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_text_repel(aes(label = sig_label), size = 2.8, max.overlaps = 20,
                  show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(x = "log2 fold change (ARHG-mutant vs WT)", y = "-log10 adjusted p-value",
       color = NULL,
       title = "ARHG-mutant vs WT (CONFOUNDED: all mutants are de novo)",
       subtitle = "Hypothesis catalog only — no causal inference valid") +
  theme_classic(base_size = 12)
dev.off()

# ── Step 8: ARHG expression ───────────────────────────────────────────────────
cat("== Step 8: ARHG expression ==\n")

arhg_ensembl <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = arhg_genes,
         column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
)
no_mapping    <- arhg_genes[is.na(arhg_ensembl)]
has_mapping   <- arhg_ensembl[!is.na(arhg_ensembl)]
in_matrix     <- intersect(has_mapping, rownames(vst_mat))
not_in_matrix <- setdiff(has_mapping, rownames(vst_mat))
ens_to_sym    <- setNames(names(has_mapping), has_mapping)

cat("  No Ensembl mapping (aliases/pseudogenes):",
    paste(no_mapping, collapse = ", "), "\n")
cat("  Has Ensembl ID but below expression threshold:",
    paste(ens_to_sym[not_in_matrix], collapse = ", "), "\n")
cat("  Detected in matrix:", length(in_matrix), "of", length(has_mapping), "\n")

arhg_expr <- as.data.frame(t(vst_mat[in_matrix, , drop = FALSE]))
colnames(arhg_expr) <- ens_to_sym[colnames(arhg_expr)]
arhg_expr$sample_id      <- rownames(arhg_expr)
arhg_expr$disease_status  <- meta[rownames(arhg_expr), "disease_status"]
arhg_expr$arhg_mutant     <- meta[rownames(arhg_expr), "arhg_mutant"]

present_syms <- intersect(names(has_mapping), colnames(arhg_expr))
arhg_long <- arhg_expr %>%
  filter(!is.na(disease_status), disease_status != "relapse") %>%
  pivot_longer(cols = all_of(present_syms), names_to = "gene", values_to = "vst")

pdf(file.path(fig_dir, "arhg_expression_upenn.pdf"), width = 16, height = 6)
print(
  ggplot(arhg_long, aes(x = gene, y = vst, color = disease_status)) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.5,
                 position = position_dodge(0.4)) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "VST expression", color = "Status",
         title = "ARHG family expression — UPenn cohort")
)
dev.off()

# ── Step 9: Save counts matrix ────────────────────────────────────────────────
write.csv(as.data.frame(counts_mat),
          file.path(out_dir, "upenn_counts_matrix.csv"))

cat("\n== Done ==\n")
cat("Outputs in:", out_dir, "and", fig_dir, "\n")
