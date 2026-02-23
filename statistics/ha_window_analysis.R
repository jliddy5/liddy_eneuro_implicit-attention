# ==============================================================================
# ha_window_analysis.R
# ==============================================================================
# Purpose: Bayesian between-subjects comparison of hand angle using a robust
#          Student-t model. Compares group means (ST, DT, DTF) within a
#          specified cycle window.
#
# Usage:   Modify the window configuration (lines 24-25) to analyze different
#          phases of the experiment:
#
#          Window Name       | Cycles  | Purpose
#          ------------------|---------|----------------------------------------
#          "Baseline"        | 6:10    | Verify group-level equivalence (Methods)
#          "EarlyLearning"   | 13:17   | Early adaptation effects (Results)
#          "LateLearning"    | 46:50   | Late adaptation effects (Results)
#          "Aftereffect"     | 51      | No feedback (Results)
#
# Output:  Posterior summaries (group means, contrasts, effect sizes, Type S/M)
#          and exports posterior draws to results/ folder.
#
# ==============================================================================

library(bayesplot)
library(here)
library(openxlsx2)
library(rstan)
library(tidybayes)
library(tidyverse)

source(here("utils", "rstan_diagnostics.R"))
source(here("utils", "compute_type_s_error.R"))
source(here("utils", "compute_type_m_error.R"))

# Configure window analysis -------------------------------------------------- #
window_name <- "Baseline"
window_idx <- 6:10

# Load data ------------------------------------------------------------------ #
df <- read_xlsx(here("data", "data_reaching.xlsx"))

# Convert id and group to factors
df$id <- as.factor(df$id)
df$group <- factor(df$group, levels = c("ST", "DT", "DTF"))


# Filter cycles and compute mean per participant
df <- df |>
  filter(cycle %in% window_idx) |>
  select(id, group, cycle, ha)

df_means <- df |>
  group_by(id, group) |>
  summarise(
    ha = mean(ha, na.rm = TRUE),
    .groups = "drop"
  )

# Group-level summary statistics
group_summary <- df_means %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(ha, na.rm = TRUE),
    median = median(ha, na.rm = TRUE),
    sd = sd(ha, na.rm = TRUE),
    min = min(ha, na.rm = TRUE),
    max = max(ha, na.rm = TRUE),
    .groups = "drop"
  )
group_summary

# Model ---------------------------------------------------------------------- #
# Prepare data
data_stan <- list(
  N     = as.integer(length(df_means$ha)),            # Number of observations
  G     = as.integer(length(unique(df_means$group))), # Number of groups
  y     = df_means$ha,                                # Observations
  group = as.integer(df_means$group),                 # Group identifiers
  n_per_group = as.integer(table(df_means$group)),    # Number of observations per group
  meanY = mean(df_means$ha),                          # Full-sample mean
  sdY   = sd(df_means$ha)                             # Full-sample SD
)

# Compile model
model_stan <- stan_model(file = here("models", "model_tdist_betweensubjects.stan"))

# Fit model (seed for reproducibility: 1501991256)
model_fit <- sampling(model_stan,
                      data = data_stan,
                      iter = 3500, 
                      warmup = 1000, 
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 15),
                      seed = 1501991256)

# Run diagnostics
rstan_diagnostics(model_fit)

# Posterior predictive check
y_sim <- rstan::extract(model_fit, pars = "y_sim")[[1]]
y_obs <- data_stan$y
group_obs <- factor(data_stan$group, labels = levels(df_means$group))

ppc_ecdf_overlay_grouped(
  y    = y_obs,
  yrep = y_sim[1:500, ],
  group = group_obs
)

# Summarize results ---------------------------------------------------------- #

# Extract posterior draws
post_draws <- as_draws_df(model_fit)

# Group means and 89% HDIs
post_means <- post_draws |>
  select(starts_with("mu[")) |>
  set_names(levels(df$group)) |>
  mutate(draws = row_number()) |>
  pivot_longer(-draws, names_to = "group", values_to = "mu")

means_summary <- post_means |>
  group_by(group) |>
  median_hdi(mu, .width = 0.89) |>
  mutate(across(where(is.numeric), \(x) round(x, 1))) |>
  select(group, mu, .lower, .upper)
means_summary

# Directional probability (per contrast)
prob_DT_gt_ST <- mean(post_draws$`mu[2]` > post_draws$`mu[1]`)
prob_DTF_gt_ST <- mean(post_draws$`mu[3]` > post_draws$`mu[1]`)
prob_DTF_gt_DT <- mean(post_draws$`mu[3]` > post_draws$`mu[2]`)
cat("Pr(DT > ST) =", round(prob_DT_gt_ST, 2), "\n")
cat("Pr(DTF > ST) =", round(prob_DTF_gt_ST, 2), "\n")
cat("Pr(DTF > DT) =", round(prob_DTF_gt_DT, 2), "\n")

# Compute type S error (per contrast)
post_samples <- rstan::extract(model_fit)
type_s_DT_ST <- compute_type_s_error(post_samples, 2, 1)
type_s_DTF_ST <- compute_type_s_error(post_samples, 3, 1)
type_s_DTF_DT <- compute_type_s_error(post_samples, 3, 2)
cat("Type S (DT - ST) =", round(type_s_DT_ST, 2), "\n")
cat("Type S (DTF - ST) =", round(type_s_DTF_ST, 2), "\n")
cat("Type S (DTF - DT) =", round(type_s_DTF_DT, 2), "\n")

# Compute type M error (per contrast)
type_m_DT_ST <- compute_type_m_error(post_samples, 2, 1)
type_m_DTF_ST <- compute_type_m_error(post_samples, 3, 1)
type_m_DTF_DT <- compute_type_m_error(post_samples, 3, 2)
cat("Type M (DT - ST) =", round(type_m_DT_ST$median_type_m, 2), 
    "[", round(type_m_DT_ST$type_m_HDI[1], 2), ",", round(type_m_DT_ST$type_m_HDI[2], 2), "]\n")
cat("Type M (DTF - ST) =", round(type_m_DTF_ST$median_type_m, 2), 
    "[", round(type_m_DTF_ST$type_m_HDI[1], 2), ",", round(type_m_DTF_ST$type_m_HDI[2], 2), "]\n")
cat("Type M (DTF - DT) =", round(type_m_DTF_DT$median_type_m, 2), 
    "[", round(type_m_DTF_DT$type_m_HDI[1], 2), ",", round(type_m_DTF_DT$type_m_HDI[2], 2), "]\n")

# Group mean differences and 89% HDIs
post_contrasts <- post_draws |>
  transmute(
    draws = row_number(),
    `DT - ST` = `mu[2]` - `mu[1]`,
    `DTF - ST` = `mu[3]` - `mu[1]`,
    `DTF - DT` = `mu[3]` - `mu[2]`
  ) |>
  pivot_longer(-draws, names_to = "contrast", values_to = "diff")

contrasts_summary <- post_contrasts |>
  group_by(contrast) |>
  median_hdi(diff, .width = 0.89) |>
  as_tibble() |>
  mutate(across(where(is.numeric), \(x) round(x, 1))) |>
  select(contrast, diff, .lower, .upper)
contrasts_summary

# Effect sizes (Cohen's d)
# Pooled SD per draw
n_per_group <- table(df_means$group)

post_d <- post_draws |>
  transmute(
    draws = row_number(),
    # Pooled SD for each contrast
    pooled_DT_ST = sqrt(((n_per_group[1] - 1) * `sigma[1]`^2 + 
                           (n_per_group[2] - 1) * `sigma[2]`^2) / 
                          (n_per_group[1] + n_per_group[2] - 2)),
    pooled_DTF_ST = sqrt(((n_per_group[1] - 1) * `sigma[1]`^2 + 
                            (n_per_group[3] - 1) * `sigma[3]`^2) / 
                           (n_per_group[1] + n_per_group[3] - 2)),
    pooled_DTF_DT = sqrt(((n_per_group[2] - 1) * `sigma[2]`^2 + 
                            (n_per_group[3] - 1) * `sigma[3]`^2) / 
                           (n_per_group[2] + n_per_group[3] - 2)),
    # Cohen's d
    `DT - ST` = (`mu[2]` - `mu[1]`) / pooled_DT_ST,
    `DTF - ST` = (`mu[3]` - `mu[1]`) / pooled_DTF_ST,
    `DTF - DT` = (`mu[3]` - `mu[2]`) / pooled_DTF_DT
  ) |>
  select(draws, starts_with("DT")) |>
  pivot_longer(-draws, names_to = "contrast", values_to = "d")

d_summary <- post_d |>
  group_by(contrast) |>
  median_hdi(d, .width = 0.89) |>
  as_tibble() |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  select(contrast, d, .lower, .upper)
d_summary

# Export posterior draws ----------------------------------------------------- #

# Create results directory if it doesn't exist
if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}

# Create file name dynamically
output_file_name <- here("results", paste0("posterior_ha_", tolower(window_name), ".xlsx"))

# Write file
write_xlsx(
  list(
    data        = df_means,
    mean        = post_means,
    mean_diff   = post_contrasts,
    effect_size = post_d
  ),
  file = output_file_name
)
