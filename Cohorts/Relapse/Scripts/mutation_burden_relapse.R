library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

VAF_THRESHOLD <- 0.02

local_env <- new.env()
source(file.path(BASE_DIR, "General/artifact_genes.txt"), local = local_env)
artifact_genes <- local_env$artifact_genes

# ── 1. Load and filter ─────────────────────────────────────────────────────────
df <- read_excel(file.path(BASE_DIR, "Cohorts/Relapse/relapse_pass_variants_combined.xlsx"))

df_somatic <- df |>
  filter(
    Timepoint %in% c("diagnosis", "relapse"),
    !Gene %in% artifact_genes,
    !Patient_Group %in% c("Patient_23", "Patient_25")   # hypermutators
  )

# ── 2. Deduplicate and apply VAF threshold ─────────────────────────────────────
# One row per Pair_ID × Gene × HGVS × Timepoint; keep only variants present
# (T-N VAF >= 2%) at that timepoint
variants <- df_somatic |>
  group_by(Pair_ID, Patient_Group, Gene, HGVS, Timepoint) |>
  summarise(vaf = max(`T-N VAF`, na.rm = TRUE), .groups = "drop") |>
  filter(vaf >= VAF_THRESHOLD)

# ── 3. Count mutations per patient per timepoint ───────────────────────────────
burden <- variants |>
  group_by(Pair_ID, Patient_Group, Timepoint) |>
  summarise(n_mutations = n_distinct(paste(Gene, HGVS)), .groups = "drop") |>
  mutate(Timepoint = factor(
    case_when(Timepoint == "diagnosis" ~ "Diagnosis",
              Timepoint == "relapse"   ~ "Relapse"),
    levels = c("Diagnosis", "Relapse")
  ))

# ── 4. Paired Wilcoxon test ────────────────────────────────────────────────────
# Pivot wide so both timepoints are on the same row; patients with 0 mutations
# at a timepoint won't appear in burden — fill those with 0
burden_wide <- burden |>
  pivot_wider(names_from = Timepoint, values_from = n_mutations, values_fill = 0)

dx   <- burden_wide |> arrange(Pair_ID) |> pull(Diagnosis)
rel  <- burden_wide |> arrange(Pair_ID) |> pull(Relapse)
diffs <- rel - dx

# ── 4a. Assumption diagnostics ────────────────────────────────────────────────
cat(sprintf("Pairs increased: %d | decreased: %d | tied: %d\n",
            sum(diffs > 0), sum(diffs < 0), sum(diffs == 0)))

# Pearson skewness of differences (|value| > 1 warrants using sign test instead)
skew_val <- (mean(diffs) - median(diffs)) / sd(diffs)
cat(sprintf("Skewness of differences (Pearson): %.3f\n", skew_val))

# Sign test as sensitivity check — no distributional assumptions
sign_p <- binom.test(sum(diffs > 0), sum(diffs != 0), p = 0.5)$p.value
cat(sprintf("Sign test (sensitivity):            p = %.4f\n", sign_p))

# ── 4b. Wilcoxon with confidence interval ─────────────────────────────────────
wil <- wilcox.test(rel, dx, paired = TRUE, conf.int = TRUE)
cat(sprintf(
  "\nWilcoxon signed-rank: W = %.0f, p = %.4f\n  Hodges-Lehmann estimate: %.1f [95%% CI: %.1f, %.1f]\n",
  wil$statistic, wil$p.value,
  wil$estimate, wil$conf.int[1], wil$conf.int[2]
))

# ── 4c. Effect size ────────────────────────────────────────────────────────────
# r = |Z| / sqrt(N) where Z is derived from the W statistic directly.
# Thresholds: < 0.1 negligible | 0.1–0.3 small | 0.3–0.5 medium | > 0.5 large
n_pairs <- sum(diffs != 0)   # non-tied pairs only
w       <- as.numeric(wil$statistic)
e_w     <- n_pairs * (n_pairs + 1) / 4
var_w   <- n_pairs * (n_pairs + 1) * (2 * n_pairs + 1) / 24
z_approx <- (w - e_w) / sqrt(var_w)
r_effect  <- abs(z_approx) / sqrt(n_pairs)
magnitude <- cut(r_effect,
                 breaks = c(0, 0.1, 0.3, 0.5, Inf),
                 labels = c("negligible", "small", "medium", "large"),
                 right  = FALSE)
cat(sprintf("Effect size r = %.3f (%s)\n", r_effect, magnitude))

# Rebuild long burden from wide so zero-count patients are included in the plot
burden <- burden_wide |>
  pivot_longer(c(Diagnosis, Relapse), names_to = "Timepoint", values_to = "n_mutations") |>
  mutate(Timepoint = factor(Timepoint, levels = c("Diagnosis", "Relapse")))

# ── 5. Summary stats ───────────────────────────────────────────────────────────
burden |>
  summarise(
    mean   = mean(n_mutations),
    sd     = sd(n_mutations),
    median = median(n_mutations),
    IQR    = IQR(n_mutations),
    n      = n(),
    .by    = Timepoint
  )

# ── 6. Plot ────────────────────────────────────────────────────────────────────
p_label <- sprintf("Wilcoxon p = %.4f", wil$p.value)
y_max   <- max(burden$n_mutations, na.rm = TRUE)

p <- ggplot(burden, aes(x = Timepoint, y = n_mutations, group = Pair_ID)) +
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.4, size = 1.6) +
  stat_summary(aes(group = 1), fun = median, geom = "line",     linewidth = 1.4) +
  stat_summary(aes(group = 1), fun = median, geom = "point",    size = 3) +
  stat_summary(aes(group = 1), fun.data = mean_cl_boot, geom = "errorbar", width = 0.08) +
  annotate("text",
           x     = 1.5,
           y     = y_max * 1.05,
           label = p_label,
           size  = 4) +
  scale_y_continuous(
    limits = c(0, y_max * 1.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.3))) +
  labs(
    x     = NULL,
    y     = "Mutations per patient",
    title = "Mutation burden from Diagnosis to Relapse"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.margin       = margin(4, 6, 4, 6),
    axis.ticks.length = unit(2, "pt"),
    axis.title.y      = element_text(margin = margin(r = 6)),
    axis.text.x       = element_text(margin = margin(t = 2)),
    panel.grid        = element_blank(),
    legend.position   = "none"
  )

print(p)

ggsave(file.path(BASE_DIR, "Cohorts/Relapse/mutation_burden_relapse.pdf"),
       plot = p, width = 4, height = 5, units = "in")
ggsave(file.path(BASE_DIR, "Cohorts/Relapse/mutation_burden_relapse.png"),
       plot = p, width = 4, height = 5, units = "in", dpi = 300)

message("Saved: Cohorts/Relapse/mutation_burden_relapse.pdf/.png")
