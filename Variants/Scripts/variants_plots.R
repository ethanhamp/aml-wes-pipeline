BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
variants <- read_excel(file.path(BASE_DIR, "Variants/arhg_variants.xlsx"))

totals <- variants %>%
  count(gene, name = "total_per_gene")

ggplot(
  variants %>% left_join(totals, by = "gene") %>%
    mutate(gene = fct_reorder(gene, total_per_gene, .desc = TRUE)),
  aes(x = gene, fill = effect_simple)
) +
  geom_bar() +
  coord_flip() +
  labs(x = "Gene", y = "Variant count", fill = "Mutation Type") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(expand = c(0, 0))
  
data_perc <- variants %>%
    mutate(percentage = effect_simple / sum(effect_simple))
  
ggplot(variants, aes(x= type, fill = effect_simple)) + 
    geom_bar(position = "stack") +
    labs(fill = "Effect")

counts <- variants %>%
  count(type, effect_simple, name = "n") %>%
  group_by(type) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

ggplot(counts, aes(x = type, y = pct, fill = effect_simple)) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  geom_text(
    aes(label = percent(pct, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    color = "white", size = 3
  ) +
  labs(x = "Class", y = "Percent", fill = "Mutation Type")

missense <- variants %>% filter(effect_simple == "missense")

colors <- c("ambiguous" = "blue", "benign" = "green", "pathogenic" = "red")

ggplot(
  missense %>% left_join(totals, by = "gene") %>%
    mutate(gene = fct_reorder(gene, total_per_gene, .desc = TRUE)),
  aes(x = gene, fill = alphamissense_pred)
) +
  geom_bar() +
  coord_flip() +
  labs(x = "Gene", y = "Variant count", fill = "Predicted effect") +
  theme_minimal(base_size = 12) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = colors, na.value = "grey80", drop = FALSE)

