library(ggplot2)
library(grid)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))

# Canvas dimensions
W <- 14
H <- 9

# Helper to draw a rounded-corner box
box <- function(xmin, xmax, ymin, ymax, fill, color, alpha = 1) {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           fill = fill, color = color, linewidth = 0.6, alpha = alpha)
}

# Helper for centered text
lbl <- function(x, y, label, size, fontface = "plain", color = "black", hjust = 0.5, vjust = 0.5) {
  annotate("text", x = x, y = y, label = label, size = size,
           fontface = fontface, color = color, hjust = hjust, vjust = vjust,
           lineheight = 1.1)
}

# Arrow segment
arr <- function(x, xend, y, yend, color = "#555555") {
  annotate("segment", x = x, xend = xend, y = y, yend = yend,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           color = color, linewidth = 0.6)
}

# ---- Colors ----
red_fill   <- "#FDECEA"
red_border <- "#C0392B"
red_head   <- "#FADBD8"

blue_fill   <- "#EBF5FB"
blue_border <- "#2471A3"
blue_head   <- "#D6EAF8"

navy_fill  <- "#1A2A4A"

# ---- Layout: x zones ----
# Left column center: x = 3.5, right column center: x = 10.5
# Figure runs x: 0.5 - 13.5

p <- ggplot() +
  xlim(0, W) + ylim(0, H) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# ===== TITLE =====
p <- p +
  lbl(7, 8.65, "ARHG Mutations: Independent or Convergent Oncogenic Events?",
      size = 6.5, fontface = "bold", color = "black")

# ===== LEFT COLUMN =====
lx1 <- 0.6; lx2 <- 6.6; lcx <- (lx1 + lx2) / 2

# Header box (red)
p <- p +
  box(lx1, lx2, 7.6, 8.35, red_head, red_border) +
  lbl(lcx, 8.1, "52% of ARHG-mutated patients", size = 5, fontface = "bold", color = red_border) +
  lbl(lcx, 7.82, "Co-mutated with NRAS / KRAS / KIT", size = 4, color = "#7B241C")

# Box A: ARHG mutation
p <- p +
  box(lx1 + 0.4, lx2 - 0.4, 6.6, 7.3, red_fill, red_border) +
  lbl(lcx, 7.05, "ARHG mutation", size = 4.5, fontface = "bold", color = red_border) +
  lbl(lcx, 6.78, "GEF gain-of-function  or  GAP loss-of-function", size = 3.3, color = "#7B241C")

# Arrow
p <- p + arr(lcx, lcx, 6.6, 6.35, red_border)

# Box B: NRAS / KRAS
p <- p +
  box(lx1 + 0.4, lx2 - 0.4, 5.6, 6.3, red_fill, red_border) +
  lbl(lcx, 6.03, "NRAS / KRAS", size = 4.5, fontface = "bold", color = red_border) +
  lbl(lcx, 5.76, "Direct RAS activation", size = 3.3, color = "#7B241C")

# Arrow
p <- p + arr(lcx, lcx, 5.6, 5.35, red_border)

# Box C: Rho GTPase
p <- p +
  box(lx1 + 0.4, lx2 - 0.4, 4.6, 5.3, red_fill, red_border) +
  lbl(lcx, 5.03, "\u2193 Rho GTPase activity", size = 4.5, fontface = "bold", color = red_border) +
  lbl(lcx, 4.76, "RhoA / Rac1 / Cdc42", size = 3.3, color = "#7B241C")

# Arrow
p <- p + arr(lcx, lcx, 4.6, 4.35, red_border)

# Box D: RAS activity
p <- p +
  box(lx1 + 0.4, lx2 - 0.4, 3.6, 4.3, red_fill, red_border) +
  lbl(lcx, 4.03, "\u2193 RAS activity", size = 4.5, fontface = "bold", color = red_border) +
  lbl(lcx, 3.76, "GTP-bound RAS", size = 3.3, color = "#7B241C")

# Arrow
p <- p + arr(lcx, lcx, 3.6, 3.35, red_border)

# Box E: PI3K / MAPK
p <- p +
  box(lx1 + 0.4, lx2 - 0.4, 2.6, 3.3, red_fill, red_border) +
  lbl(lcx, 3.03, "PI3K / MAPK activation", size = 4.5, fontface = "bold", color = red_border) +
  lbl(lcx, 2.76, "Convergent downstream signaling", size = 3.3, color = "#7B241C")

# ===== RIGHT COLUMN =====
rx1 <- 7.4; rx2 <- 13.4; rcx <- (rx1 + rx2) / 2

# Header box (blue)
p <- p +
  box(rx1, rx2, 7.6, 8.35, blue_head, blue_border) +
  lbl(rcx, 8.1, "48% of ARHG-mutated patients", size = 5, fontface = "bold", color = blue_border) +
  lbl(rcx, 7.82, "No co-mutations in canonical RAS pathway genes", size = 4, color = "#1A5276")

# Box A: ARHG mutation
p <- p +
  box(rx1 + 0.4, rx2 - 0.4, 6.6, 7.3, blue_fill, blue_border) +
  lbl(rcx, 7.05, "ARHG mutation", size = 4.5, fontface = "bold", color = blue_border) +
  lbl(rcx, 6.78, "GEF gain-of-function  or  GAP loss-of-function", size = 3.3, color = "#1A5276")

# Arrow
p <- p + arr(rcx, rcx, 6.6, 6.35, blue_border)

# Box B: Rho GTPase
p <- p +
  box(rx1 + 0.4, rx2 - 0.4, 5.6, 6.3, blue_fill, blue_border) +
  lbl(rcx, 6.03, "\u2193 Rho GTPase activity", size = 4.5, fontface = "bold", color = blue_border) +
  lbl(rcx, 5.76, "RhoA / Rac1 / Cdc42", size = 3.3, color = "#1A5276")

# Arrow
p <- p + arr(rcx, rcx, 5.6, 5.35, blue_border)

# Box C: Cross talk
p <- p +
  box(rx1 + 0.4, rx2 - 0.4, 4.6, 5.3, blue_fill, blue_border) +
  lbl(rcx - 0.3, 5.03, "Cross talk \u2192 PI3K / MAPK", size = 4.5, fontface = "bold", color = blue_border) +
  lbl(rcx - 0.3, 4.76, "Intrinsic pathway signaling", size = 3.3, color = "#1A5276") +
  # small note on right side
  lbl(rx2 - 0.5, 5.08, "via PI3K\nRalGDS\nNF-\u03baB", size = 2.5, color = "#1A5276", hjust = 0.5)

# Arrow from cross talk down to bottom box
p <- p + arr(rcx, rcx, 4.6, 4.35, blue_border)

# ===== BOTTOM CONVERGENCE BOX =====
# Both arrows converge to this
p <- p +
  arr(lcx, lcx, 2.6, 2.35, "#666666") +
  arr(rcx, rcx, 4.6, 4.35, blue_border)

# Redraw final arrow for left since we already drew above
# Bottom wide box
p <- p +
  box(1.2, 12.8, 1.55, 2.3, "#F0F0F0", "#666666") +
  lbl(7, 2.05, "Proliferation  \u00b7  Survival  \u00b7  Treatment resistance",
      size = 5, fontface = "bold", color = "#333333") +
  lbl(7, 1.76, "Same downstream phenotype — independent of direct RAS mutation",
      size = 3.3, color = "#555555")

# Arrows from column boxes to bottom box
p <- p +
  arr(lcx, lcx, 2.6, 2.3, "#666666") +
  arr(rcx, rcx, 4.6, 4.3, "#666666")

# ===== HYPOTHESIS BAR =====
p <- p +
  box(0.5, 13.5, 0.1, 1.35, navy_fill, navy_fill) +
  lbl(7, 0.88,
      "Hypothesis: ARHG mutations dysregulate Rho GTPase activity, feeding into RAS-adjacent signaling even without direct RAS mutations",
      size = 4, fontface = "bold.italic", color = "white") +
  lbl(7, 0.45,
      "ARHG mutations may be independent oncogenic events \u2014 not simply passengers to NRAS or KRAS",
      size = 3.3, color = "#BDC3C7")

# ===== SAVE =====
ggsave(
  file.path(BASE_DIR, "Presentations/arhg_ras_connection.png"),
  plot = p,
  width = W, height = H, dpi = 180, bg = "white"
)

cat("Saved arhg_ras_connection.png\n")
