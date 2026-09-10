# Fig. S1d-h: Regional cell body quantification
# Author: Ruohan Zhong

# Packages
library(tidyverse)
library(stringr)

# Load data
median_size <- read_csv("Aggregated.csv")
density <- read_csv("DensityDistribution.csv")
cells <- read_csv("AreaResults.csv")

# Define region order and labels
region_order <- c("mouth", "ring", "tentacles")

region_labels <- c(
  "mouth" = "Mouth",
  "ring" = "Ring",
  "tentacles" = "Tentacle"
)

region_colors <- c(
  "mouth" = "#4DBFD3",
  "ring" = "#dd4c2a",
  "tentacles" = "#7d49c9"
)

# Define biological replicate order
replicate_order <- c(
  "Data/projections\\2_juvenile5_oral001_Max.tif",
  "Data/projections\\3_juvenile6_oral001_Max.tif",
  "Data/projections\\1_juvenile1_oral_Max.tif",
  "Data/projections\\7_juvenile2_oral001_MAX.tif",
  "Data/projections\\9_animal2_oral_MAX.tif",
  "Data/projections\\4_R5-3_juvenile 1_oral_Max.tif",
  "Data/projections\\5_animal1_oral_MAX.tif",
  "Data/projections\\6_animal2_oral_MAX.tif",
  "Data/projections\\8_juvenile3_oral001_MAX.tif"
)

replicate_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")

median_size <- median_size %>%
  mutate(
    replicate = factor(fname, levels = replicate_order, labels = replicate_labels),
    state = factor(state, levels = region_order)
  )

density <- density %>%
  mutate(
    replicate = factor(fname, levels = replicate_order, labels = replicate_labels),
    state = factor(state, levels = region_order)
  )

cells <- cells %>%
  mutate(
    replicate = factor(fname, levels = replicate_order, labels = replicate_labels),
    state = factor(state, levels = region_order)
  )

# Fig. S1d - Individual cell body area distributions
p_cell_area_violin <- ggplot(
  cells,
  aes(x = replicate, y = area, fill = state)
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    linewidth = 0.25,
    alpha = 0.75
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.2,
    color = "black"
  ) +
  facet_wrap(
    ~ state,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(region_labels)
  ) +
  scale_fill_manual(values = region_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = "Biological replicate",
    y = expression("Cell body area ("*mu*"m"^2*")")
  ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_cell_area_violin

# Fig. S1e - Median cell body area by region
p_median_size <- ggplot(
  median_size,
  aes(x = state, y = area, group = replicate)
) +
  geom_line(color = "gray65", linewidth = 0.6, alpha = 0.8) +
  geom_point(
    aes(fill = state),
    shape = 21,
    size = 2,
    color = "black",
    stroke = 0.3
  ) +
  stat_summary(
    aes(group = state),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.3,
    color = "black"
  ) +
  scale_x_discrete(labels = region_labels) +
  scale_fill_manual(values = region_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = expression("Median cell body area ("*mu*"m"^2*")")
  ) +
  theme(
    legend.position = "none"
  )

p_median_size

# Fig. S1f - Region area by region
p_region_area <- ggplot(
  density,
  aes(x = state, y = region_area, group = replicate)
) +
  geom_line(color = "gray65", linewidth = 0.6, alpha = 0.8) +
  geom_point(
    aes(fill = state),
    shape = 21,
    size = 2,
    color = "black",
    stroke = 0.3
  ) +
  stat_summary(
    aes(group = state),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.3,
    color = "black"
  ) +
  scale_x_discrete(labels = region_labels) +
  scale_fill_manual(values = region_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = expression("Region area ("*mu*"m"^2*")")
  ) +
  theme(
    legend.position = "none"
  )

p_region_area

# Fig. S1g - Cell count density by region
density <- density %>%
  mutate(
    count_density_per_10000um2 = count_density * 10000
  )

p_count_density <- ggplot(
  density,
  aes(x = state, y = count_density_per_10000um2, group = replicate)
) +
  geom_line(color = "gray65", linewidth = 0.6, alpha = 0.8) +
  geom_point(
    aes(fill = state),
    shape = 21,
    size = 2,
    color = "black",
    stroke = 0.3
  ) +
  stat_summary(
    aes(group = state),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.3,
    color = "black"
  ) +
  scale_x_discrete(labels = region_labels) +
  scale_fill_manual(values = region_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = expression("Density by cell count per 10,000 "*mu*"m"^2)
  ) +
  theme(
    legend.position = "none"
  )

p_count_density

# Fig. S1h - Cell body area density by region
density <- density %>%
  mutate(
    area_density_percent = area_density * 100
  )

p_area_density <- ggplot(
  density,
  aes(x = state, y = area_density_percent, group = replicate)
) +
  geom_line(color = "gray65", linewidth = 0.6, alpha = 0.8) +
  geom_point(
    aes(fill = state),
    shape = 21,
    size = 2,
    color = "black",
    stroke = 0.3
  ) +
  stat_summary(
    aes(group = state),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.3,
    color = "black"
  ) +
  scale_x_discrete(labels = region_labels) +
  scale_fill_manual(values = region_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Density by cell body area (%)"
  ) +
  theme(
    legend.position = "none"
  )

p_area_density

# Session information
sessionInfo()
