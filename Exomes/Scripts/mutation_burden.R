library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(readxl)
library(readr)
library(stringr)
library(ggrepel)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

local_env <- new.env()
source(file.path(BASE_DIR, "General/artifact_genes.txt"), local = local_env)
artifact_genes <- local_env$artifact_genes

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
  group_by(Group, Gene, HGVS, Timepoint, `AD Total`) %>%
  summarise(VAF = max(VAF, na.rm = TRUE), .groups = "drop")

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

# ---- Settings ----
PRESENCE_VAF <- 0      # detection threshold
METRIC <- "count"         # "count" or "sum"
# -------------------

# 1) Build per-patient mutation burden at Dx and Rel
burden <- pi_version %>%
  select(Group, Gene, VAF_diagnosis, VAF_relapse) %>%
  pivot_longer(cols = starts_with("VAF_"),
               names_to = "timepoint", values_to = "VAF") %>%
  mutate(timepoint = recode(timepoint,
                            VAF_diagnosis = "Diagnosis",
                            VAF_relapse   = "Relapse")) %>%
  filter(!is.na(VAF)) %>%
  group_by(Group, timepoint) %>%
  summarise(
    burden_count = sum(VAF > PRESENCE_VAF),
    burden_sum   = sum(VAF),
    .groups = "drop"
  ) %>%
  mutate(timepoint = factor(timepoint, levels = c("Diagnosis","Relapse")))

# Choose which metric to plot
plot_df <- burden %>%
  transmute(
    Group, timepoint,
    value = if (METRIC == "sum") burden_sum else burden_count
  )

# 2) Optional: paired test for overall shift
dx  <- plot_df %>% filter(timepoint == "Diagnosis") %>% arrange(Group) %>% pull(value)
rel <- plot_df %>% filter(timepoint == "Relapse")   %>% arrange(Group) %>% pull(value)
wil <- wilcox.test(rel, dx, paired = TRUE)  # nonparametric paired test
cat(sprintf("Wilcoxon paired test (Relapse vs Diagnosis): W=%s, p=%.3g\n",
            wil$statistic, wil$p.value))

# 3) Plot: per-patient lines + cohort trend
p_val <- wil$p.value
y_max <- max(plot_df$value, na.rm = TRUE)

ggplot(plot_df, aes(x = timepoint, y = value, group = Group)) +
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.4, size = 1.6) +
  stat_summary(aes(group = 1), fun = median, geom = "line", size = 1.4) +
  stat_summary(aes(group = 1), fun = median, geom = "point", size = 3) +
  stat_summary(aes(group = 1), fun.data = mean_cl_boot, geom = "errorbar", width = 0.08) +
  annotate(
    "text", 
    x = 1.5,                          # midway between Diagnosis and Relapse
    y = max(plot_df$value) * 1.03,    # slightly above your data range
    label = sprintf("Wilcoxon p = %.3g", p_val),
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, y_max * 1.06),                   # 0 anchored at bottom
    expand = expansion(mult = c(0, 0.02))   # no extra white space below; tiny at top
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(0.3, 0.3))      # remove left/right padding around the two x points
  ) +
  labs(
    x = NULL,
    y = if (METRIC == "sum") "Sum of VAFs per patient"
    else sprintf("Mutations per patient", PRESENCE_VAF),
    title = "Mutation burden from Diagnosis to Relapse"
  ) +
  theme_classic(base_size = 10) +
  theme (plot.margin      = margin(4, 6, 4, 6),  # squeeze outer whitespace (t, r, b, l) in pts
         axis.ticks.length= unit(2, "pt"),
         axis.title.x     = element_text(margin = margin(t = 4)),
         axis.title.y     = element_text(margin = margin(r = 6)),
         axis.text.x      = element_text(margin = margin(t = 2)),
         panel.grid       = element_blank(),
         legend.position  = "none")

avg_by_time <- burden %>%
  summarise(
    mean_mutations = mean(burden_count, na.rm = TRUE),
    sd             = sd(burden_count,   na.rm = TRUE),
    median         = median(burden_count, na.rm = TRUE),
    IQR            = IQR(burden_count,    na.rm = TRUE),
    n_patients     = dplyr::n(),
    .by = timepoint
  )

avg_by_time

ggplot(pi_version, aes(x = VAF_relapse)) +
  geom_histogram(bins = 50)

ggplot(pi_version, aes(x = VAF_diagnosis)) +
  geom_density()

ks.test(pi_version$VAF_diagnosis, pi_version$VAF_relapse) # p = 7.011e-09, sig diff

d1 <- density(pi_version$VAF_diagnosis, na.rm = TRUE)
d2 <- density(pi_version$VAF_relapse, na.rm = TRUE)

overlap_density <- function(x1, x2, n = 2048) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)  # need data
  
  # common range/grid so no NA from extrapolation
  rng <- range(c(x1, x2))
  d1  <- density(x1, from = rng[1], to = rng[2], n = n, na.rm = TRUE)
  d2  <- density(x2, from = rng[1], to = rng[2], n = n, na.rm = TRUE)
  
  x  <- d1$x
  y1 <- d1$y
  y2 <- d2$y
  dx <- x[2] - x[1]
  
  # normalize to integrate to 1 (use dx, not raw sums)
  A1 <- sum(y1 * dx); if (!is.finite(A1) || A1 == 0) return(NA_real_)
  A2 <- sum(y2 * dx); if (!is.finite(A2) || A2 == 0) return(NA_real_)
  y1 <- y1 / A1
  y2 <- y2 / A2
  
  # overlap area in [0,1]
  sum(pmin(y1, y2) * dx)
}

# Example:
ov <- overlap_density(pi_version$VAF_diagnosis, pi_version$VAF_relapse)
ov # 0.8655188 very similar shape; lots of overlap


ymax <- 1.05 * max(d1$y, d2$y)   # add a little headroom

plot(d1, col="blue", lwd=2, main="Diagnosis vs Relapse VAF Distributions", xlab="VAF", 
     ylab="Density", ylim=c(0, ymax))
lines(d2, col="red", lwd=2)
legend("topright", legend=c("Diagnosis","Relapse"), col=c("blue","red"), lwd=2)

dev.off()

long <- pi_version %>%
  pivot_longer(c(VAF_diagnosis, VAF_relapse, VAF_germline), names_to="Timepoint", values_to="VAF") %>%
  filter(is.finite(VAF))

dx  <- long %>% filter(Timepoint=="VAF_diagnosis") %>% arrange(VAF) %>% mutate(idx = row_number())
rel <- long %>% filter(Timepoint=="VAF_relapse")   %>% arrange(VAF) %>% mutate(idx = row_number())

ggplot() +
  geom_segment(data = dx,  aes(x = idx, y = 0, xend = idx, yend = VAF, color = "Diagnosis"),
               linewidth = 0.3, alpha = 0.6) +
  geom_point(  data = dx,  aes(x = idx, y = VAF, color = "Diagnosis"), size = 1.5) +
  geom_segment(data = rel, aes(x = idx, y = 0, xend = idx, yend = VAF, color = "Relapse"),
               linewidth = 0.3, alpha = 0.6) +
  geom_point(  data = rel, aes(x = idx, y = VAF, color = "Relapse"), size = 1.5) +
  scale_color_manual(
    name   = "Timepoint",
    values = c("Diagnosis" = "#1f77b4", "Relapse" = "#d62728")
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 3))) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Variant index", y = "VAF", title = "Diagnosis vs Relapse VAFs") +
  theme_minimal(base_size = 11)
