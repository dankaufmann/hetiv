#!/usr/bin/env Rscript

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(hexSticker)

fill <- "#164F86"
border <- "#14275B"

x <- seq(-4.5, 4.5, length.out = 500)
thin_density <- dnorm(x, sd = 1.1)
fat_tail_density <- dt(x / 1.7, df = 2.8) / 1.7

scale_density <- function(y) y / max(y)

curve_data <- rbind(
  data.frame(
    x = x,
    y = 0.5 + 0.32 * scale_density(fat_tail_density),
    distribution = "fat tail"
  ),
  data.frame(
    x = x,
    y = 0.5 + 0.35 * scale_density(thin_density),
    distribution = "thin tail"
  )
)

distribution_plot <- ggplot(curve_data, aes(x, y, colour = distribution)) +
  geom_line(linewidth = 1.4, lineend = "round") +
  scale_colour_manual(
    values = c("fat tail" = "#FFFFFF", "thin tail" = border),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-4.1, 4.1), ylim = c(0.43, 0.88), expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

write_sticker <- function(filename) {
  hexSticker::sticker(
    subplot = distribution_plot,
    package = "hetiv",
    filename = filename,
    s_x = 1,
    s_y = 1.1,
    s_width = 1.8,
    s_height = 1.7,
    p_x = 1,
    p_y = 0.98,
    p_size = 28,
    p_family = "sans",
    p_fontface = "plain",
    p_color = "#FFFFFF",
    h_fill = fill,
    h_color = border,
    # h_size = 2.8,
    # dpi = 600,
    white_around_sticker = FALSE
  )
}

write_sticker("man/figures/logo.png")