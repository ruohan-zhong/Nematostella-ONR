# Fig. 6e-f: Feeding behavior statistics
# Author: Ruohan Zhong

# Packages
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

# Load and prepare grasping data
grasp <- read.csv("grasp_data.csv") |>
  as_tibble() |>
  filter(method == "pipette") |>
  mutate(
    genotype   = factor(genotype, levels = c("WT", "MUT")),
    animalID   = factor(animalID),
    trial      = factor(trial, levels = c(1, 2, 3, 4)),
    method     = factor(method),
    starvation = factor(starvation)
  )

# Load and prepare swallowing data
swallow <- read.csv("swallow_data.csv") |>
  as_tibble() |>
  select(-starts_with("Unnamed"), -starts_with("..."), everything()) |>
  mutate(
    genotype   = factor(genotype, levels = c("WT", "MUT")),
    animalID   = factor(animalID),
    trial      = factor(trial, levels = c("1", "2", "3", "4")),
    method     = factor(method),
    starvation = factor(starvation)
  )

# Create output directory
dir.create("stats_output", showWarnings = FALSE)

# Fig. 6e - Grasping
mA_g <- lmer(
  Tgrasp ~ genotype * trial + (1 | animalID),
  data = grasp
)

sum_mA_g <- summary(mA_g)
coef_mA_g_df <- as.data.frame(sum_mA_g$coefficients)

write.csv(
  coef_mA_g_df,
  file = "stats_output/grasp_mA_summary_coefficients.csv",
  row.names = TRUE
)

anova_mA_g <- anova(mA_g)
anova_mA_g_df <- as.data.frame(anova_mA_g)

write.csv(
  anova_mA_g_df,
  file = "stats_output/grasp_mA_anova_table.csv",
  row.names = TRUE
)

emmA_g <- emmeans(mA_g, ~ genotype | trial)
emmA_g_df <- as.data.frame(emmA_g)

write.csv(
  emmA_g_df,
  file = "stats_output/grasp_mA_emmeans_genotype_by_trial.csv",
  row.names = FALSE
)

conA_g <- contrast(emmA_g, method = "revpairwise")
conA_g_df <- as.data.frame(conA_g)

write.csv(
  conA_g_df,
  file = "stats_output/grasp_mA_contrasts_genotype_by_trial.csv",
  row.names = FALSE
)

# Fig. 6f - Swallowing: Model A (genotype x trial)
mA <- lmer(
  Tswallow ~ genotype * trial + (1 | animalID),
  data = swallow
)

sum_mA <- summary(mA)
coef_mA_df <- as.data.frame(sum_mA$coefficients)

write.csv(
  coef_mA_df,
  file = "stats_output/swallow_mA_summary_coefficients.csv",
  row.names = TRUE
)

anova_mA <- anova(mA)
anova_mA_df <- as.data.frame(anova_mA)

write.csv(
  anova_mA_df,
  file = "stats_output/swallow_mA_anova_table.csv",
  row.names = TRUE
)

emmA <- emmeans(mA, ~ genotype | trial)
emmA_df <- as.data.frame(emmA)

write.csv(
  emmA_df,
  file = "stats_output/swallow_mA_emmeans_genotype_by_trial.csv",
  row.names = FALSE
)

conA <- contrast(emmA, method = "revpairwise")
conA_df <- as.data.frame(conA)

write.csv(
  conA_df,
  file = "stats_output/swallow_mA_contrasts_genotype_by_trial.csv",
  row.names = FALSE
)

# Fig. 6f - Swallowing: Model B (+ delivery method)
mB <- lmer(
  Tswallow ~ genotype * trial + method + (1 | animalID),
  data = swallow
)

sum_mB <- summary(mB)
coef_mB_df <- as.data.frame(sum_mB$coefficients)

write.csv(
  coef_mB_df,
  file = "stats_output/swallow_mB_summary_coefficients.csv",
  row.names = TRUE
)

anova_mB <- anova(mB)
anova_mB_df <- as.data.frame(anova_mB)

write.csv(
  anova_mB_df,
  file = "stats_output/swallow_mB_anova_table.csv",
  row.names = TRUE
)

emmB <- emmeans(mB, ~ genotype | trial)
emmB_df <- as.data.frame(emmB)

write.csv(
  emmB_df,
  file = "stats_output/swallow_mB_emmeans_genotype_by_trial.csv",
  row.names = FALSE
)

conB <- contrast(emmB, method = "revpairwise")
conB_df <- as.data.frame(conB)

write.csv(
  conB_df,
  file = "stats_output/swallow_mB_contrasts_genotype_by_trial.csv",
  row.names = FALSE
)

# Fig. 6f - Swallowing: Model C (+ delivery method and starvation)
mC <- lmer(
  Tswallow ~ genotype * trial + method + starvation + (1 | animalID),
  data = swallow
)

sum_mC <- summary(mC)
coef_mC_df <- as.data.frame(sum_mC$coefficients)

write.csv(
  coef_mC_df,
  file = "stats_output/swallow_mC_summary_coefficients.csv",
  row.names = TRUE
)

anova_mC <- anova(mC)
anova_mC_df <- as.data.frame(anova_mC)

write.csv(
  anova_mC_df,
  file = "stats_output/swallow_mC_anova_table.csv",
  row.names = TRUE
)

emmC <- emmeans(mC, ~ genotype | trial)
emmC_df <- as.data.frame(emmC)

write.csv(
  emmC_df,
  file = "stats_output/swallow_mC_emmeans_genotype_by_trial.csv",
  row.names = FALSE
)

conC <- contrast(emmC, method = "revpairwise")
conC_df <- as.data.frame(conC)

write.csv(
  conC_df,
  file = "stats_output/swallow_mC_contrasts_genotype_by_trial.csv",
  row.names = FALSE
)

# Session information
sessionInfo()
