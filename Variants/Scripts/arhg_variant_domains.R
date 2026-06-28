library(readxl)
library(dplyr)
library(ggplot2)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
arhg_variants <- read_xlsx(file.path(BASE_DIR, "Variants/arhg_variants_clean_hg19.xlsx"))

arhg_types <- arhg_variants %>%
  mutate(
    Gene = str_trim(Hugo_Symbol),
    Class = case_when(
      str_detect(Hugo_Symbol, regex("^ARHGAP", ignore_case = TRUE)) ~ "GAP",
      str_detect(Hugo_Symbol, regex("^ARHGEF", ignore_case = TRUE)) ~ "GEF",
      TRUE ~ "Other"
    )
  )

gaps <- arhg_types %>% filter(Class == "GAP")
gefs <- arhg_types %>% filter(Class == "GEF")

domain_colors <- c(
  "RhoGAP" = "#D55E00",
  "RhoGEF" = "#E69F00",
  "CH" = "#009E73",
  "PH" = "#0072B2",
  "FF" = "#CC79A7",
  "No Domain" = "grey"
)

counts <- gaps %>%
  mutate(Domain.hg19 = ifelse(is.na(Domain.hg19) | Domain.hg19 == "",
                              "No Domain", Domain.hg19)) %>%
  count(Domain.hg19, name = "n") %>%
  mutate(
    pct   = n / sum(n),
    label = scales::percent(pct, accuracy = 0.1)
  )

pie_variant_domains <- function(gaps, title = "GAP variant domains") {
  ggplot(counts, aes(x = "", y = n, fill = Domain.hg19)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = percent(pct, accuracy = 0.1)),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = domain_colors) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

pie_gap <- pie_variant_domains(gaps, title = "GAP variant domains")
pie_gap

counts_gef <- gefs %>%
  mutate(Domain.hg19 = ifelse(is.na(Domain.hg19) | Domain.hg19 == "",
                              "No Domain", Domain.hg19)) %>%
  count(Domain.hg19, name = "n") %>%
  mutate(
    pct   = n / sum(n),
    label = scales::percent(pct, accuracy = 0.1)
  )

pie_variant_domains_gef <- function(gefs, title = "GEF variant domains") {
  ggplot(counts_gef, aes(x = "", y = n, fill = Domain.hg19)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = percent(pct, accuracy = 0.1)),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = domain_colors) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

pie_gef <- pie_variant_domains_gef(gefs, title = "GEF variant domains")
pie_gef
