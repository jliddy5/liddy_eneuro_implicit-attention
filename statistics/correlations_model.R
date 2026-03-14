# ==============================================================================
# correlations_exp1_exp2_model.R
# ==============================================================================
# Purpose: Bayesian regression of computational model parameters (retention,
#          error sensitivity) on RSVP accuracy for DT and DTF groups combined.
#          Models four regressions using Beta regression with RSVP phase
#          accuracy as predictor.
#
# Output:  Model summaries and Bayesian R² printed to console.
#
# ==============================================================================

library(bayesplot)
library(bridgesampling)
library(brms)
library(here)
library(openxlsx2)
library(tidybayes)
library(tidyverse)

# Load model data ------------------------------------------------------------ #
df_model <- read_xlsx(here("data", "data_modeling.xlsx")) |>
  select(id, group, a, b) |>
  mutate(
    id = as.factor(id),
    group = as.factor(group)
  ) |>
  filter(group %in% c("DT", "DTF")) |>
  mutate(
    id = droplevels(id),
    group = droplevels(group)
  )
head(df_model)

# Load RSVP data ------------------------------------------------------------- #
df_rsvp <- read_xlsx(here("data", "data_rsvp.xlsx")) |>
  mutate(
    id = as.factor(id),
    group = as.factor(group),
    phase = as.factor(phase)
  ) |>
  filter(group %in% c("DT", "DTF")) |>
  mutate(
    id = droplevels(id),
    group = droplevels(group)
  ) |>
  pivot_wider(
    names_from = phase,
    values_from = accuracy,
    names_prefix = "rsvp_"
  )

# Merge data ----------------------------------------------------------------- #
df <- df_model |>
  left_join(df_rsvp, by = c("id", "group"))

# Priors -------------------------------------------------------------------- #
priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "phi")
)

# Retention (a) ============================================================= #

# Model 1: Retention ~ Baseline RSVP
fit_a_baseline <- brm(
  formula = a ~ rsvp_baseline,
  data = df,
  family = Beta(),
  prior = priors,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

summary(fit_a_baseline, prob = 0.89)
bayes_R2(fit_a_baseline, summary = TRUE, probs = c(0.055, 0.945))

# Model 2: Retention ~ Late RSVP
fit_a_late <- brm(
  formula = a ~ rsvp_late,
  data = df,
  family = Beta(),
  prior = priors,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

summary(fit_a_late, prob = 0.89)
bayes_R2(fit_a_late, summary = TRUE, probs = c(0.055, 0.945))

# Error sensitivity (b) ===================================================== #

# Model 3: Error sensitivity ~ Baseline RSVP
fit_b_baseline <- brm(
  formula = b ~ rsvp_baseline,
  data = df,
  family = Beta(),
  prior = priors,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

summary(fit_b_baseline, prob = 0.89)
bayes_R2(fit_b_baseline, summary = TRUE, probs = c(0.055, 0.945))

# Model 4: Error sensitivity ~ Late RSVP
fit_b_late <- brm(
  formula = b ~ rsvp_late,
  data = df,
  family = Beta(),
  prior = priors,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

summary(fit_b_late, prob = 0.89)
bayes_R2(fit_b_late, summary = TRUE, probs = c(0.055, 0.945))
