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

# Type III ANOVA with Satterthwaite approximation
anova_mA_g <- anova(
  mA_g,
  type = 3,
  ddf = "Satterthwaite"
)

source_fig6e_global_stats <- as.data.frame(anova_mA_g) %>%
  tibble::rownames_to_column("Effect") %>%
  mutate(
    Panel = "Fig. 6e",
    Model = "genotype * trial"
  ) %>%
  relocate(
    Panel,
    Model,
    Effect
  )

# Estimated marginal means and genotype contrasts within each trial
emmA_g <- emmeans(
  mA_g,
  ~ genotype | trial,
  lmer.df = "kenward-roger"
)

conA_g <- contrast(
  emmA_g,
  method = "revpairwise",
  adjust = "none"
)

source_fig6e_pairwise_stats <- as.data.frame(
  summary(
    conA_g,
    infer = c(TRUE, TRUE),
    adjust = "none"
  )
) %>%
  mutate(
    Panel = "Fig. 6e",
    Model = "genotype * trial"
  ) %>%
  relocate(
    Panel,
    Model
  )

# Write Source Data
write.csv(
  source_fig6e_global_stats,
  "SourceData_Fig6e_global_stats.csv",
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  source_fig6e_pairwise_stats,
  "SourceData_Fig6e_pairwise_stats.csv",
  row.names = FALSE,
  quote = TRUE
)

# Fig. 6f - Swallowing

# Model A: genotype x trial
mA <- lmer(
  Tswallow ~ genotype * trial + (1 | animalID),
  data = swallow
)

# Model B: genotype x trial + delivery method
mB <- lmer(
  Tswallow ~ genotype * trial + method + (1 | animalID),
  data = swallow
)

# Model C: genotype x trial + delivery method + starvation
mC <- lmer(
  Tswallow ~ genotype * trial + method + starvation + (1 | animalID),
  data = swallow
)

# global Type III ANOVA

get_global_stats <- function(
    model,
    model_name
) {
  
  result <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  )
  
  as.data.frame(result) %>%
    tibble::rownames_to_column("Effect") %>%
    mutate(
      Panel = "Fig. 6f",
      Model = model_name
    ) %>%
    relocate(
      Panel,
      Model,
      Effect
    )
}

# WT vs mutant within each trial
get_pairwise_stats <- function(
    model,
    model_name
) {
  
  emm <- emmeans(
    model,
    ~ genotype | trial,
    lmer.df = "kenward-roger"
  )
  
  con <- contrast(
    emm,
    method = "revpairwise",
    adjust = "none"
  )
  
  as.data.frame(
    summary(
      con,
      infer = c(TRUE, TRUE),
      adjust = "none"
    )
  ) %>%
    mutate(
      Panel = "Fig. 6f",
      Model = model_name
    ) %>%
    relocate(
      Panel,
      Model
    )
}

# Global statistics for all swallowing models
source_fig6f_global_stats <- bind_rows(
  get_global_stats(
    mA,
    "Model A: genotype * trial"
  ),
  get_global_stats(
    mB,
    "Model B: genotype * trial + delivery method"
  ),
  get_global_stats(
    mC,
    "Model C: genotype * trial + delivery method + starvation"
  )
)

# Pairwise statistics for all swallowing models
source_fig6f_pairwise_stats <- bind_rows(
  get_pairwise_stats(
    mA,
    "Model A: genotype * trial"
  ),
  get_pairwise_stats(
    mB,
    "Model B: genotype * trial + delivery method"
  ),
  get_pairwise_stats(
    mC,
    "Model C: genotype * trial + delivery method + starvation"
  )
)

# Write Source Data
write.csv(
  source_fig6f_global_stats,
  "SourceData_Fig6f_global_stats.csv",
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  source_fig6f_pairwise_stats,
  "SourceData_Fig6f_pairwise_stats.csv",
  row.names = FALSE,
  quote = TRUE
)

# Session information
sessionInfo()
