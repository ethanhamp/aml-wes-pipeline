library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
exomes_effect_pred <- read_excel(file.path(BASE_DIR, "Variants/arhg_exomes_pred_effects.xlsx"))

arhg_typed <- exomes_effect_pred %>%
  mutate(
    gene = str_trim(gene),
    Class = case_when(
      str_detect(gene, regex("^ARHGAP", ignore_case = TRUE)) ~ "GAP",
      str_detect(gene, regex("^ARHGEF", ignore_case = TRUE)) ~ "GEF",
      TRUE ~ "Other"
    )
  )

gaps <- arhg_typed %>% filter(Class == "GAP")
gefs <- arhg_typed %>% filter(Class == "GEF")

effect_colors <- c(
  "Deleterious" = "#D55E00",
  "Possibly Damaging" = "#E69F00",
  "Benign" = "#009E73",
  "No Prediction" = "#0072B2",
  "Tolerated" = "#CC79A7"
)

pie_predicted_effect <- function(gaps, title = "Predicted effects") {
  # count + percent
  counts <- gaps %>%
    filter(!is.na(siftpred2) & siftpred2 != "") %>%
    mutate(siftpred2 = str_to_title(siftpred2)) %>%
    count(siftpred2, name = "n") %>%
    mutate(pct = n / sum(n))
  
  ggplot(counts, aes(x = "", y = n, fill = siftpred2)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = percent(pct, accuracy = 0.1)),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = effect_colors) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

pie_gap <- pie_predicted_effect(gaps, title = "ARHGAPs: Predicted Effects")
pie_gap

pie_predicted_effect <- function(gefs, title = "Predicted effects") {
  # count + percent
  counts <- gefs %>%
    filter(!is.na(siftpred2) & siftpred2 != "") %>%
    mutate(siftpred2 = str_to_title(siftpred2)) %>%
    count(siftpred2, name = "n") %>%
    mutate(pct = n / sum(n))
  
  ggplot(counts, aes(x = "", y = n, fill = siftpred2)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = percent(pct, accuracy = 0.1)),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = effect_colors) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

pie_gef <- pie_predicted_effect(gefs, title = "ARHGEFs: Predicted Effects")
pie_gef
