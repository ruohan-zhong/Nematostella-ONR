# Fig. 6e-f: Feeding behavior
# Author: Ruohan Zhong

# Packages
library(dplyr)
library(ggplot2)

# Plot colors
my_cols <- c(WT = "#2166ac", MUT = "#f4881c")

# Load and prepare grasping data
grasp <- read.csv("grasp_data.csv") |>
  as_tibble() |>
  filter(method == "pipette") |>
  mutate(
    genotype   = factor(genotype, levels = c("WT", "MUT")),
    animalID   = factor(animalID),
    trial      = factor(trial, levels = c(1, 2, 3, 4)),
    trial_num  = as.numeric(as.character(trial)),
    method     = factor(method),
    starvation = factor(starvation),
    animal     = animalID,
    time       = Tgrasp
  )

# Load and prepare swallowing data
swallow <- read.csv("swallow_data.csv") |>
  as_tibble() |>
  select(-starts_with("Unnamed"), -starts_with("..."), everything()) |>
  mutate(
    genotype   = factor(genotype, levels = c("WT", "MUT")),
    animalID   = factor(animalID),
    trial      = factor(trial, levels = c("1", "2", "3", "4")),
    trial_num  = as.numeric(as.character(trial)),
    method     = factor(method),
    starvation = factor(starvation),
    animal     = animalID,
    time       = Tswallow
  )

# Fig. 6e - Grasping across feeding trials
grasp_summary <- grasp |>
  group_by(genotype, trial_num) |>
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    se_time   = sd(time, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups   = "drop"
  )

p_grasp <- ggplot() +
  geom_line(
    data = grasp,
    aes(x = trial_num, y = time, group = animal, color = genotype),
    alpha = 0.3,
    linewidth = 0.6
  ) +
  geom_line(
    data = grasp_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    linewidth = 1.2
  ) +
  geom_point(
    data = grasp_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    size = 2.2
  ) +
  geom_errorbar(
    data = grasp_summary,
    aes(
      x = trial_num,
      ymin = mean_time - se_time,
      ymax = mean_time + se_time,
      color = genotype
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_color_manual(values = my_cols) +
  scale_x_continuous(breaks = 1:4, labels = 1:4) +
  scale_y_continuous(limits = c(1, 6), breaks = c(1, 3, 5)) +
  labs(
    x = NULL,
    y = "Grasping (s)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.position = "none"
  )

p_grasp

# Fig. 6f - Swallowing across feeding trials
swallow_summary <- swallow |>
  group_by(genotype, trial_num) |>
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    se_time   = sd(time, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups   = "drop"
  )

p_swallow <- ggplot() +
  geom_line(
    data = swallow,
    aes(x = trial_num, y = time, group = animal, color = genotype),
    alpha = 0.3,
    linewidth = 0.6
  ) +
  geom_line(
    data = swallow_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    linewidth = 1.2
  ) +
  geom_point(
    data = swallow_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    size = 2.2
  ) +
  geom_errorbar(
    data = swallow_summary,
    aes(
      x = trial_num,
      ymin = mean_time - se_time,
      ymax = mean_time + se_time,
      color = genotype
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_color_manual(values = my_cols) +
  scale_x_continuous(breaks = 1:4, labels = 1:4) +
  scale_y_continuous(
    limits = c(40, 200),
    breaks = c(50, 100, 150, 200)
  ) +
  labs(
    x = "Feeding trial",
    y = "Swallowing (s)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.position = "none"
  )

p_swallow

# Session information
sessionInfo()
