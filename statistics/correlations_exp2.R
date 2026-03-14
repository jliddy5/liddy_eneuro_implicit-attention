# ==============================================================================
# correlations_exp2.R
# ==============================================================================
# Purpose: Bayesian correlations between RSVP accuracy and hand angle for the
#          DTF group (Experiment 2). Models four pairwise correlations using
#          bivariate Gaussian models with estimated residual correlation.
#
# Output:  Correlation summaries printed to console.
#
# ==============================================================================

library(brms)
library(here)
library(openxlsx2)
library(patchwork)
library(tidybayes)
library(tidyverse)

# Load reaching data --------------------------------------------------------- #
df_reach <- read_xlsx(here("data", "data_reaching.xlsx")) |>
  mutate(
    id = as.factor(id),
    group = as.factor(group)
  ) |>
  filter(group == "DTF") |>
  mutate(
    id = droplevels(id),
    group = droplevels(group)
    )

# Compute participant-level metrics
df_reach_metrics <- df_reach |>
  filter(cycle %in% c(13:17, 46:50)) |>
  mutate(window = case_when(
    cycle %in% 13:17 ~ "early",
    cycle %in% 46:50 ~ "late"
  )) |>
  group_by(id, window) |>
  summarise(
    ha_mean = mean(ha, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = window,
    values_from = ha_mean,
    names_prefix = "ha_"
  )

# Load RSVP data ------------------------------------------------------------- #
df_rsvp <- read_xlsx(here("data", "data_rsvp.xlsx")) |>
  mutate(
    id = as.factor(id),
    group = as.factor(group),
    phase = as.factor(phase)
  ) |>
  filter(group == "DTF") |>
  mutate(
    id = droplevels(id),
    group = droplevels(group)
  )

# Compute participant-level metrics
df_rsvp_metrics <- df_rsvp |>
  pivot_wider(
    names_from = phase,
    values_from = accuracy,
    names_prefix = "rsvp_"
  ) |>
  mutate(
    rsvp_drop = rsvp_baseline - rsvp_onset,
    rsvp_recovery = rsvp_early - rsvp_onset
  ) |>
  select(id, rsvp_drop, rsvp_recovery, 
         rsvp_baseline, rsvp_onset, rsvp_early, rsvp_late)

# Merge data ----------------------------------------------------------------- #
df_correlation <- df_reach_metrics |>
  left_join(df_rsvp_metrics, by = "id")

head(df_correlation)

# Model ---------------------------------------------------------------------- #

# Correlation 1: RSVP drop × Early HA
# Does initial disruption predict early learning magnitude?
fit_drop_early <- brm(
  bf(mvbind(ha_early, rsvp_drop) ~ 1) + set_rescor(TRUE),
  data = df_correlation,
  family = gaussian(),
  prior = c(
    prior(lkj(1), class = rescor)
  ),
  iter = 3500,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

# Correlation 2: RSVP recovery × Early HA
# Does recovery from disruption predict early learning magnitude?
fit_recovery_early <- brm(
  bf(mvbind(ha_early, rsvp_recovery) ~ 1) + set_rescor(TRUE),
  data = df_correlation,
  family = gaussian(),
  prior = c(
    prior(lkj(1), class = rescor)
  ),
  iter = 3500,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

# Correlation 3: Final RSVP level × Late HA
# Does final RSVP performance relate to late adaptation?
fit_late_late <- brm(
  bf(mvbind(ha_late, rsvp_late) ~ 1) + set_rescor(TRUE),
  data = df_correlation,
  family = gaussian(),
  prior = c(
    prior(lkj(1), class = rescor)
  ),
  iter = 3500,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

# Correlation 4: Baseline RSVP × Late HA
# Does baseline RSVP predict late adaptation? (expect no correlation)
fit_baseline_late <- brm(
  bf(mvbind(ha_late, rsvp_baseline) ~ 1) + set_rescor(TRUE),
  data = df_correlation,
  family = gaussian(),
  prior = c(
    prior(lkj(1), class = rescor)
  ),
  iter = 3500,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1501991256
)

# Check summaries
summary(fit_drop_early)
summary(fit_recovery_early)
summary(fit_late_late)
summary(fit_baseline_late)

# Extract correlation posteriors with 89% HDI
cor_drop_early <- fit_drop_early |>
  posterior::as_draws_df() |>
  select(starts_with("rescor")) |>
  pivot_longer(everything(), names_to = "param", values_to = "r") |>
  median_hdi(r, .width = 0.89) |>
  mutate(contrast = "rsvp_drop × ha_early")

cor_recovery_early <- fit_recovery_early |>
  posterior::as_draws_df() |>
  select(starts_with("rescor")) |>
  pivot_longer(everything(), names_to = "param", values_to = "r") |>
  median_hdi(r, .width = 0.89) |>
  mutate(contrast = "rsvp_recovery × ha_early")

cor_late_late <- fit_late_late |>
  posterior::as_draws_df() |>
  select(starts_with("rescor")) |>
  pivot_longer(everything(), names_to = "param", values_to = "r") |>
  median_hdi(r, .width = 0.89) |>
  mutate(contrast = "rsvp_late × ha_late")

cor_baseline_late <- fit_baseline_late |>
  posterior::as_draws_df() |>
  select(starts_with("rescor")) |>
  pivot_longer(everything(), names_to = "param", values_to = "r") |>
  median_hdi(r, .width = 0.89) |>
  mutate(contrast = "rsvp_baseline × ha_late")

# Combine and print
cor_summary <- bind_rows(cor_drop_early, cor_recovery_early, cor_late_late, cor_baseline_late) |>
  select(contrast, r, .lower, .upper) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))
cor_summary
