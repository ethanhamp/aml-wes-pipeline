library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# ── Configuration ──────────────────────────────────────────────────────────────
# Primary file: most current with domain annotation (hg38 + hg19 coords)
VARIANT_FILE <- file.path(BASE_DIR, "Variants/arhg_variants_clean_hg38.hg19.xlsx")

# Domain colors — canonical
DOMAIN_COLORS <- c(
  "RhoGEF"     = "#E69F00",
  "RhoGAP"     = "#D55E00",
  "PH"         = "#0072B2",
  "CH"         = "#009E73",
  "FF"         = "#CC79A7",
  "Rho_GDI"    = "#56B4E9",
  "PF17838"    = "#F0E442",  # ARHGEF12-specific Pfam domain
  "PF19056"    = "#999999",  # ARHGEF17-specific Pfam domain
  "non-domain" = "grey85"
)

# ── 1. Load ────────────────────────────────────────────────────────────────────
message(sprintf("Loading %s ...", VARIANT_FILE))
df_raw <- read_excel(VARIANT_FILE)
message(sprintf("  %d rows, %d columns", nrow(df_raw), ncol(df_raw)))

df <- df_raw |>
  rename(
    gene           = Hugo_Symbol,
    protein_change = Protein_Change,
    sample_id      = Sample_ID,
    mutation_type  = Mutation_Type,
    domain_hg38    = Domain.hg38,
    domain_hg19    = Domain.hg19
  ) |>
  mutate(
    gene = str_trim(gene),
    # Family assignment
    family = case_when(
      str_detect(gene, "^ARHGEF") ~ "GEF",
      str_detect(gene, "^ARHGAP") ~ "GAP",
      str_detect(gene, "^ARHGDI") ~ "GDI",
      TRUE ~ "Other"
    ),
    # Use hg38 domain as primary; fill NA from hg19
    domain = coalesce(
      if_else(is.na(domain_hg38) | domain_hg38 == "", NA_character_, domain_hg38),
      if_else(is.na(domain_hg19) | domain_hg19 == "", NA_character_, domain_hg19),
      "non-domain"
    ),
    domain = factor(domain, levels = c("RhoGEF", "RhoGAP", "PH", "CH", "FF",
                                        "Rho_GDI", "PF17838", "PF19056", "non-domain"))
  )

message(sprintf("  Families: GEF=%d, GAP=%d, GDI=%d",
  sum(df$family == "GEF"), sum(df$family == "GAP"), sum(df$family == "GDI")))

# ── 2. Domain summary ─────────────────────────────────────────────────────────
domain_summary <- df |>
  count(family, domain) |>
  group_by(family) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup()

message("\nDomain distribution by family:")
print(as.data.frame(domain_summary))

n_in_domain <- df |> filter(domain != "non-domain")
message(sprintf("\n%d of %d variants in a named domain (%.0f%%)",
  nrow(n_in_domain), nrow(df), 100 * nrow(n_in_domain) / nrow(df)))

# ── 3. Stacked bar: domain frequency by family ────────────────────────────────
p_domain_bar <- domain_summary |>
  filter(family %in% c("GEF", "GAP")) |>  # exclude GDI (n=1) for main figure
  ggplot(aes(x = family, y = n, fill = domain)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = if_else(n > 0, sprintf("%d\n(%.0f%%)", n, pct), "")),
    position = position_stack(vjust = 0.5),
    size = 3.2, lineheight = 1.0,
    color = ifelse(domain_summary |> filter(family %in% c("GEF","GAP")) |> pull(domain) == "non-domain",
                   "grey40", "white")
  ) +
  scale_fill_manual(values = DOMAIN_COLORS, name = "Domain", drop = FALSE) +
  scale_x_discrete(labels = c("GEF" = "GEF (ARHGEF)", "GAP" = "GAP (ARHGAP)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title  = "ARHG variant domain distribution",
    subtitle = sprintf("n=%d variants total; hg38 primary annotation", nrow(df)),
    x = NULL,
    y = "Number of variants"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 9),
    legend.position = "right",
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# ── 4. Per-gene domain bar (horizontal) ───────────────────────────────────────
p_gene_domain <- df |>
  mutate(gene = factor(gene, levels = rev(sort(unique(gene))))) |>
  ggplot(aes(x = gene, fill = domain)) +
  geom_bar(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = DOMAIN_COLORS, name = "Domain", drop = FALSE) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  facet_grid(family ~ ., scales = "free_y", space = "free_y",
    labeller = labeller(family = c(GEF = "GEF", GAP = "GAP", GDI = "GDI"))) +
  labs(
    title = "ARHG variants per gene with domain annotation",
    x     = NULL,
    y     = "Number of variants"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text.y    = element_text(face = "bold", angle = 0),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.position = "right",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

# ── 5. Domain-in vs out comparison ────────────────────────────────────────────
domain_binary <- df |>
  filter(family %in% c("GEF", "GAP")) |>
  mutate(in_domain = if_else(domain != "non-domain", "In domain", "Non-domain")) |>
  count(family, in_domain) |>
  group_by(family) |>
  mutate(pct = round(100 * n / sum(n), 1))

# Fisher's exact: is domain-hit rate different between GEF and GAP?
ct_wide <- domain_binary |>
  select(family, in_domain, n) |>
  pivot_wider(names_from = in_domain, values_from = n, values_fill = 0)
ct <- matrix(
  c(ct_wide[["In domain"]], ct_wide[["Non-domain"]]),
  nrow = 2,
  dimnames = list(ct_wide$family, c("In domain", "Non-domain"))
)

if (all(c("In domain", "Non-domain") %in% colnames(ct))) {
  ft <- fisher.test(ct[, c("In domain", "Non-domain")])
  message(sprintf("\nFisher's exact (domain-hit rate GEF vs GAP): p=%.3f, OR=%.2f",
    ft$p.value, ft$estimate))

  p_binary <- domain_binary |>
    ggplot(aes(x = family, y = pct, fill = in_domain)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
      position = position_stack(vjust = 0.5),
      size = 3.5, lineheight = 1.0,
      color = c("white","grey40","white","grey40")
    ) +
    scale_fill_manual(
      values = c("In domain" = "#2166AC", "Non-domain" = "grey85"),
      name = NULL
    ) +
    scale_x_discrete(labels = c("GEF" = "GEF (ARHGEF)", "GAP" = "GAP (ARHGAP)")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    annotate("text", x = 1.5, y = 105,
      label = sprintf("Fisher's p = %.3f", ft$p.value),
      size = 3.2, color = "grey30") +
    labs(
      title = "Domain-targeting variants: GEF vs GAP",
      x = NULL, y = "Percentage of variants (%)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      axis.line.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# ── 6. Save ───────────────────────────────────────────────────────────────────
dir.create(file.path(BASE_DIR, "Variants/Figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(BASE_DIR, "Variants/Figures/domain_bar_gef_vs_gap.pdf"), plot = p_domain_bar,
  width = 5, height = 4.5, units = "in")
ggsave(file.path(BASE_DIR, "Variants/Figures/domain_bar_gef_vs_gap.png"), plot = p_domain_bar,
  width = 5, height = 4.5, units = "in", dpi = 300)

ggsave(file.path(BASE_DIR, "Variants/Figures/domain_bar_per_gene.pdf"), plot = p_gene_domain,
  width = 6, height = 7, units = "in")
ggsave(file.path(BASE_DIR, "Variants/Figures/domain_bar_per_gene.png"), plot = p_gene_domain,
  width = 6, height = 7, units = "in", dpi = 300)

if (exists("p_binary")) {
  ggsave(file.path(BASE_DIR, "Variants/Figures/domain_binary_gef_vs_gap.pdf"), plot = p_binary,
    width = 4.5, height = 4, units = "in")
  ggsave(file.path(BASE_DIR, "Variants/Figures/domain_binary_gef_vs_gap.png"), plot = p_binary,
    width = 4.5, height = 4, units = "in", dpi = 300)
}

message("\nSaved:")
message("  Variants/Figures/domain_bar_gef_vs_gap.pdf/.png")
message("  Variants/Figures/domain_bar_per_gene.pdf/.png")
message("  Variants/Figures/domain_binary_gef_vs_gap.pdf/.png")
