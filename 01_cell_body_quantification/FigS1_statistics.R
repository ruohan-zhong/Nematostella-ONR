# Fig. S1e, g, h: Regional cell body quantification statistics
# Author: Ruohan Zhong

# Packages
library(tidyverse)
library(rstatix)

# Load data
median_size <- read_csv("Aggregated.csv")
density <- read_csv("DensityDistribution.csv")

# Define region and biological replicate order
region_order <- c("mouth", "ring", "tentacles")

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
    state = factor(state, levels = region_order),
    count_density_per_10000um2 = count_density * 10000,
    area_density_percent = area_density * 100
  )

# Friedman tests
f_e <- median_size %>% friedman_test(area ~ state | replicate)
f_g <- density %>% friedman_test(count_density_per_10000um2 ~ state | replicate)
f_h <- density %>% friedman_test(area_density_percent ~ state | replicate)

global_results <- bind_rows(
  mutate(f_e, metric = "median_size"),
  mutate(f_g, metric = "count_density"),
  mutate(f_h, metric = "area_density")
)

# Paired Wilcoxon signed-rank tests with Holm correction
pwc_e <- median_size %>% wilcox_test(area ~ state, paired = TRUE, p.adjust.method = "holm")
pwc_g <- density %>% wilcox_test(count_density_per_10000um2 ~ state, paired = TRUE, p.adjust.method = "holm")
pwc_h <- density %>% wilcox_test(area_density_percent ~ state, paired = TRUE, p.adjust.method = "holm")

pairwise_results <- bind_rows(
  mutate(pwc_e, metric = "median_size"),
  mutate(pwc_g, metric = "count_density"),
  mutate(pwc_h, metric = "area_density")
)

# Export results
write_csv(global_results, "supplementary_global_stats.csv")
write_csv(pairwise_results, "supplementary_pairwise_stats.csv")

# Session information
sessionInfo()
