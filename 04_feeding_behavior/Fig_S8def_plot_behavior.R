# Fig. S8d-f: Feeding behavior
# Author: Ruohan Zhong

# Packages
library(tidyverse)

# Plot colors
my_cols <- c(WT = "#2166ac", MUT = "#f4881c")

# Load and prepare 10-trial grasping data
grasp_10 <- read.csv("grasp_pipette_10.csv") |>
  as_tibble() |>
  filter(method == "pipette") |>
  mutate(
    genotype   = factor(genotype, levels = c("WT", "MUT")),
    animalID   = factor(animalID),
    trial      = factor(trial, levels = 1:10),
    trial_num  = as.numeric(as.character(trial)),
    method     = factor(method),
    starvation = factor(starvation),
    animal     = animalID,
    time       = Tgrasp
  )

# Load and prepare 10-trial swallowing data
swallow_10 <- read.csv("swallow_pipette_10.csv") |>
  as_tibble() |>
  filter(method == "pipette") |>
  mutate(
    genotype   = factor(genotype, levels = c("WT", "MUT")),
    animalID   = factor(animalID),
    trial      = factor(trial, levels = 1:10),
    trial_num  = as.numeric(as.character(trial)),
    method     = factor(method),
    starvation = factor(starvation),
    animal     = animalID,
    time       = Tswallow
  )

# Fig. S8d - Grasping across 10 feeding trials
grasp_10_summary <- grasp_10 |>
  group_by(genotype, trial_num) |>
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    se_time   = sd(time, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups   = "drop"
  )

p_grasp_10 <- ggplot() +
  geom_line(
    data = grasp_10,
    aes(x = trial_num, y = time, group = animal, color = genotype),
    alpha = 0.3,
    linewidth = 0.6
  ) +
  geom_line(
    data = grasp_10_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    linewidth = 1.2
  ) +
  geom_point(
    data = grasp_10_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    size = 2.2
  ) +
  geom_errorbar(
    data = grasp_10_summary,
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
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  scale_y_continuous(limits = c(1, 9), breaks = c(1, 3, 5, 7, 9)) +
  labs(
    x = "Feeding trial",
    y = "Grasping (s)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

p_grasp_10

# Fig. S8e - Swallowing across 10 feeding trials
swallow_10_summary <- swallow_10 |>
  group_by(genotype, trial_num) |>
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    se_time   = sd(time, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups   = "drop"
  )

p_swallow_10 <- ggplot() +
  geom_line(
    data = swallow_10,
    aes(x = trial_num, y = time, group = animal, color = genotype),
    alpha = 0.3,
    linewidth = 0.6
  ) +
  geom_line(
    data = swallow_10_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    linewidth = 1.2
  ) +
  geom_point(
    data = swallow_10_summary,
    aes(x = trial_num, y = mean_time, color = genotype),
    size = 2.2
  ) +
  geom_errorbar(
    data = swallow_10_summary,
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
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  scale_y_continuous(
    limits = c(45, 155),
    breaks = c(50, 75, 100, 125, 150)
  ) +
  labs(
    x = "Feeding trial",
    y = "Swallowing (s)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.position = "none"
  )

p_swallow_10

# Fig. S8f - Grasping with pipette and whisker delivery
df <- read.csv("grasp_pipette_vs_whisker.csv")

df_long <- df |>
  pivot_longer(
    cols = everything(),
    names_to = "group",
    values_to = "grasp_time"
  ) |>
  separate(group, into = c("method", "genotype"), sep = "_") |>
  mutate(
    genotype = recode(genotype, wt = "WT", mut = "MUT"),
    method = recode(method, pipette = "Pipette", whisker = "Whisker"),
    group_order = factor(
      paste(genotype, method, sep = "_"),
      levels = c(
        "WT_Pipette",
        "MUT_Pipette",
        "WT_Whisker",
        "MUT_Whisker"
      )
    ),
    genotype = factor(genotype, levels = c("WT", "MUT"))
  )

p_box <- ggplot(
  df_long,
  aes(x = group_order, y = grasp_time, fill = genotype)
) +
  geom_boxplot(
    alpha = 0.35,
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.3,
    alpha = 0.4,
    size = 1.2,
    color = "black"
  ) +
  scale_fill_manual(values = my_cols) +
  scale_x_discrete(
    labels = c(
      "WT_Pipette"  = "WT\npipette",
      "MUT_Pipette" = "MUT\npipette",
      "WT_Whisker"  = "WT\nwhisker",
      "MUT_Whisker" = "MUT\nwhisker"
    )
  ) +
  scale_y_continuous(limits = c(1, 20)) +
  labs(
    x = NULL,
    y = "Grasping (s)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text.x = element_text(vjust = 0.8),
    legend.position = "none"
  )

p_box

# Session information
sessionInfo()
