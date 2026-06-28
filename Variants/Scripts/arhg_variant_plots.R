library(ggplot2)
library(readxl)
install.packages("ggrepel")
library(ggrepel)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
arhg <- read_excel(file.path(BASE_DIR, "Variants/arhg_sizes.xlsx"))

ggplot(arhg, aes(x=size.hg38, y=protein.size, color = class)) +
  geom_point() +
  theme_bw() +
  geom_text_repel(
    data = subset(arhg, variants.in.domain),
    aes(label = genes),
    max.overlaps = 100,
    na.rm = TRUE
  ) +
  labs(x="gene size (nt)", y="protein size (aa)", 
       title = "ARHG family genes: variants found in domain")

model <- lm(arhg$protein.size ~ arhg$size.hg38)
summary(model) # linear model not appropriate for this data

r2 <- summary(model)$r.squared

ggplot(arhg, aes(x=chromosome)) +
  geom_bar() +
  theme_bw() +
  labs(x="chromosomes", y="gene size (nt)",
       title = "ARHG family genes: gene size by chromosome")
  