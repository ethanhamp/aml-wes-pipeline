library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(forcats)

BASE_DIR    <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
fig_dir     <- file.path(BASE_DIR, "RNAseq/UPenn/Figures")
seq_summary <- file.path(BASE_DIR, "Cohorts/UPenn/upenn_seq_summary.xlsx")
exome_file  <- file.path(BASE_DIR, "Cohorts/UPenn/Output/upenn_all_combined.xlsx")

arhg_gefs  <- c(paste0("ARHGEF", c(1:19, 25, 26, 28)))
arhg_gaps  <- paste0("ARHGAP", c(1:10, "11A", "11B", 12:45))
arhg_genes <- c(arhg_gefs, arhg_gaps)

# ── Load metadata ──────────────────────────────────────────────────────────────
manifest    <- read.delim(file.path(BASE_DIR, "Cohorts/Manifests/master_sample_manifest.tsv"),
                          stringsAsFactors = FALSE)
upenn_man   <- manifest[!is.na(manifest$cohort_UPenn) & manifest$cohort_UPenn == TRUE, ]
upenn_man$sample_id <- as.character(upenn_man$sample_id)

# Normalize timepoint → disease_status labels consistent with DE script
upenn_man$disease_status <- upenn_man$timepoint
upenn_man$disease_status[upenn_man$disease_status == "Dx"]       <- "de novo"
upenn_man$disease_status[upenn_man$disease_status == "blast"]    <- "de novo"
upenn_man$disease_status[upenn_man$disease_status == "Relapse"]  <- "relapse"
upenn_man$disease_status[upenn_man$disease_status == "Unknown"]  <- "unknown"

# RNAseq sample IDs from the counts matrix columns
counts_mat      <- read.csv(file.path(BASE_DIR, "RNAseq/UPenn/Output/upenn_counts_matrix.csv"),
                            row.names = 1, check.names = FALSE)
rnaseq_ids      <- colnames(counts_mat)

# Exome-derived ARHG mutation status
exome           <- suppressWarnings(read_excel(exome_file))
exome_som       <- exome[exome$Call_Type == "somatic", ]
arhg_mut_ids    <- unique(as.character(exome_som$Sample[exome_som$Gene %in% arhg_genes]))
all_exome_ids   <- unique(as.character(exome$Sample))

# Build per-sample table for all 52 UPenn manifest entries
upenn_man$has_rnaseq   <- as.logical(upenn_man$has_rnaseq)
upenn_man$has_wes      <- as.logical(upenn_man$has_wes)
upenn_man$arhg_mutant  <- upenn_man$sample_id %in% arhg_mut_ids

# Fallback status for RNA-seq samples not in manifest (127)
ss <- suppressMessages(read_excel(seq_summary))
ss_df <- rbind(
  data.frame(sample_id = as.character(as.integer(ss[["sample.id...2"]])),
             disease_status = ss[["status...3"]], stringsAsFactors = FALSE),
  data.frame(sample_id = as.character(as.integer(ss[["sample.id...11"]])),
             disease_status = ss[["status...12"]], stringsAsFactors = FALSE)
)
ss_df <- ss_df[!is.na(ss_df$sample_id) & !is.na(ss_df$disease_status), ]
ss_df <- ss_df[!duplicated(ss_df$sample_id), ]

# Build RNA-seq–only metadata table (43 samples)
meta_rna <- data.frame(sample_id = rnaseq_ids, stringsAsFactors = FALSE)
man_status <- upenn_man[, c("sample_id", "disease_status", "arhg_mutant")]
meta_rna   <- merge(meta_rna, man_status, by = "sample_id", all.x = TRUE)
missing    <- is.na(meta_rna$disease_status)
if (any(missing)) {
  fill     <- ss_df[ss_df$sample_id %in% meta_rna$sample_id[missing], ]
  meta_rna$disease_status[missing] <- fill$disease_status[match(
    meta_rna$sample_id[missing], fill$sample_id)]
  meta_rna$arhg_mutant[missing]    <- meta_rna$sample_id[missing] %in% arhg_mut_ids
}
meta_rna$disease_status[is.na(meta_rna$disease_status)] <- "unknown"
meta_rna$arhg_mutant[is.na(meta_rna$arhg_mutant)]       <- FALSE

status_colors <- c(
  "de novo"    = "#4575b4",
  "refractory" = "#d73027",
  "relapse"    = "#f46d43",
  "unknown"    = "grey60"
)

# ── Figure 1: Samples by disease status ───────────────────────────────────────
status_counts <- meta_rna %>%
  count(disease_status) %>%
  mutate(disease_status = factor(disease_status,
                                  levels = c("de novo", "refractory", "relapse", "unknown")))

p1 <- ggplot(status_counts, aes(x = disease_status, y = n, fill = disease_status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = status_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of samples",
       title = "UPenn RNA-seq cohort",
       subtitle = paste0("n = ", nrow(meta_rna), " samples")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

pdf(file.path(fig_dir, "summary_disease_status.pdf"), width = 5, height = 4.5)
print(p1)
dev.off()
cat("Figure 1 written: summary_disease_status.pdf\n")

# ── Figure 2: Data availability — WES vs RNA-seq ──────────────────────────────
# Use all 52 manifest entries for this
avail_df <- upenn_man %>%
  mutate(
    data_type = case_when(
       has_rnaseq &  has_wes ~ "RNA-seq + WES",
       has_rnaseq & !has_wes ~ "RNA-seq only",
      !has_rnaseq &  has_wes ~ "WES only",
      TRUE                   ~ "Neither"
    ),
    data_type = factor(data_type,
                        levels = c("RNA-seq + WES", "RNA-seq only", "WES only", "Neither"))
  ) %>%
  count(data_type, disease_status) %>%
  mutate(disease_status = factor(disease_status,
                                  levels = c("de novo", "refractory", "relapse", "unknown")))

p2 <- ggplot(avail_df, aes(x = data_type, y = n, fill = disease_status)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = status_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Number of samples", fill = "Disease status",
       title = "Data availability - UPenn cohort",
       subtitle = paste0("n = ", nrow(upenn_man), " total samples in manifest")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 11))

pdf(file.path(fig_dir, "summary_data_availability.pdf"), width = 5.5, height = 4.5)
print(p2)
dev.off()
cat("Figure 2 written: summary_data_availability.pdf\n")

# ── Figure 3: ARHG mutation status in RNA-seq cohort ─────────────────────────
arhg_df <- meta_rna %>%
  filter(disease_status != "unknown") %>%
  mutate(
    arhg_label    = ifelse(arhg_mutant, "ARHG-mutant", "ARHG-WT"),
    disease_status = factor(disease_status,
                             levels = c("de novo", "refractory", "relapse"))
  ) %>%
  count(disease_status, arhg_label)

p3 <- ggplot(arhg_df, aes(x = disease_status, y = n, fill = arhg_label)) +
  geom_col(width = 0.6, position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("ARHG-mutant" = "#1b7837", "ARHG-WT" = "#bababa")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Number of samples", fill = NULL,
       title = "ARHG mutation status by disease status",
       subtitle = "RNA-seq cohort (n = 42, unknown excluded)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 12))

pdf(file.path(fig_dir, "summary_arhg_mutation_status.pdf"), width = 5, height = 4.5)
print(p3)
dev.off()
cat("Figure 3 written: summary_arhg_mutation_status.pdf\n")

# ── Figure 4: Combined overview — dot/bubble summary ─────────────────────────
# Lollipop-style cohort overview combining all three metrics in one panel
overview <- meta_rna %>%
  mutate(disease_status = factor(disease_status,
                                  levels = c("de novo", "refractory", "relapse", "unknown"))) %>%
  group_by(disease_status) %>%
  summarise(
    n_total    = n(),
    n_arhg_mut = sum(arhg_mutant),
    .groups = "drop"
  ) %>%
  mutate(n_arhg_wt = n_total - n_arhg_mut) %>%
  pivot_longer(cols = c(n_arhg_wt, n_arhg_mut),
               names_to = "group", values_to = "count") %>%
  mutate(group = recode(group, n_arhg_wt = "ARHG-WT", n_arhg_mut = "ARHG-mutant"),
         group = factor(group, levels = c("ARHG-WT", "ARHG-mutant")))

p4 <- ggplot(overview, aes(x = disease_status, y = count, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(data = overview %>% group_by(disease_status) %>%
              summarise(total = sum(count), .groups = "drop"),
            aes(x = disease_status, y = total, label = paste0("n=", total)),
            inherit.aes = FALSE, vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("ARHG-WT" = "#bababa", "ARHG-mutant" = "#1b7837")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = "Number of samples", fill = NULL,
       title = "UPenn RNA-seq cohort overview",
       subtitle = "Stacked by ARHG mutation status") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 12))

pdf(file.path(fig_dir, "summary_cohort_overview.pdf"), width = 5.5, height = 5)
print(p4)
dev.off()
cat("Figure 4 written: summary_cohort_overview.pdf\n")

cat("\n== Summary ==\n")
cat("RNA-seq samples total:", nrow(meta_rna), "\n")
print(table(meta_rna$disease_status))
cat("ARHG-mutant:", sum(meta_rna$arhg_mutant), "\n")
