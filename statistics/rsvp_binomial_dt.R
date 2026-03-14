# ==============================================================================
# rsvp_binomial_dt.R
# ==============================================================================
# Purpose: Bayesian binomial regression of RSVP accuracy for the DT and DTF
#          groups. Models phase-level accuracy (baseline, onset, early, late)
#          with a group interaction and exports posterior draws for figure
#          plotting.
#
# Output:  Posterior summaries and exports posterior draws to results/ folder.
#
# ==============================================================================

library(bayesplot)
library(brms)
library(emmeans)
library(here)
library(openxlsx2)
library(tidybayes)
library(tidyverse)

source(here("utils", "brms_diagnostics.R"))

# Load data ------------------------------------------------------------------ #
df <- read_xlsx(here("data", "data_rsvp.xlsx"))

# Convert id and group to factors
df$id <- as.factor(df$id)
df$group <- as.factor(df$group)

# Prep data
df <- df|>
  filter(group != "ST") |>
  mutate(
    n_trials = case_when(
      phase == "baseline" ~ 20,
      phase == "early"    ~ 20,
      phase == "late"     ~ 20,
      phase == "onset"     ~ 8
    ),
    successes = round(accuracy * n_trials),
    phase = factor(phase, levels = c("baseline", "onset", "early", "late"))
  ) |>
  droplevels()

# Quick visualization
ggplot(df, aes(x = phase, y = accuracy, color = group)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(y = "Pr(correct)", x = "Phase") +
  theme_minimal()

# Model ---------------------------------------------------------------------- #

# Seed (sample.int(1e6, 1))
model_seed <- 652398

# Priors
priors <- c(
  prior(normal(0.7, 0.5), class = "b"), # Each phase: ≈ .55-.78 on prob scale
  prior(exponential(2), class = "sd")   # Random effects ≈ ±.15–.20 on prob scale
)

# Prior predictive check
fit_score_prior <- brm(
  successes | trials(n_trials) ~ phase*group + (1 | id),
  data = df,
  family = binomial(link = "logit"),
  prior = priors,
  sample_prior = "only", 
  chains = 4, iter = 3500, warmup = 1000, cores = 4,
  seed = model_seed
)

yrep <- posterior_predict(fit_score_prior, ndraws = 10000)
yrep_prop <- sweep(yrep, 2, df$n_trials, "/")
ppc_ecdf_overlay_grouped(
  y    = df$successes / df$n_trials,
  yrep = yrep_prop[1:500, ],
  group = df$phase
)

# Fit model
fit_score <- brm(
  successes | trials(n_trials) ~ phase*group + (1 | id),
  data = df,
  family = binomial(link = "logit"),
  prior = priors,
  chains = 4, iter = 3500, warmup = 1000, cores = 4,
  control = list(max_treedepth = 15, adapt_delta = 0.95),
  seed =  model_seed
)
summary(fit_score)

# Run diagnostics
brms_diagnostics(fit_score)

# Posterior predictive check
yrep_post <- posterior_predict(fit_score, ndraws = 10000)
yrep_post_prop <- sweep(yrep_post, 2, df$n_trials, "/")
df$phase_group <- interaction(df$phase, df$group, sep = "_")
ppc_ecdf_overlay_grouped(
  y    = df$successes / df$n_trials,
  yrep = yrep_post_prop[1:500, ],
  group = df$phase_group
)

# Summarize results ---------------------------------------------------------- #

# Posterior draws: accuracy by phase and group (probability correct)
emm <- emmeans(fit_score, ~ phase | group)
emm_resp <- regrid(emm, transform = "response")
post_prob <- emm_resp |>
  gather_emmeans_draws() |>
  rename(draws = .draw, prob = .value) |>
  select(draws, phase, group, prob)

prob_summary <- post_prob |>
  group_by(phase, group) |>
  median_hdi(prob, .width = 0.89) |>
  mutate(across(where(is.numeric), ~round(., 2)))
prob_summary

# Posterior draws: change from baseline (Δ probability correct)
post_delta <- post_prob |>
  group_by(group, draws) |>
  mutate(delta = prob - prob[phase == "baseline"]) |>
  filter(phase != "baseline") |>
  ungroup()

delta_summary <- post_delta |>
  group_by(phase, group) |>
  median_hdi(delta, .width = 0.89) |>
  mutate(across(where(is.numeric), ~round(., 2)))
delta_summary

# Posterior draws: odds ratios to baseline (effect size)
post_or <- contrast(emm, method = "revpairwise", ref = "baseline", adjust = "none") |>
  gather_emmeans_draws() |>
  rename(draws = .draw, log_or = .value) |>
  mutate(or = exp(log_or)) |>
  filter(grepl(" - baseline$", contrast)) |>
  mutate(phase = sub(" - baseline", "", contrast)) |>
  ungroup() |>
  select(draws, phase, group, or)

or_summary <- post_or |>
  group_by(phase, group) |>
  median_hdi(or, .width = 0.89) |>
  mutate(across(where(is.numeric), ~round(., 2)))
or_summary

# Posterior draws: group difference at each phase (DTF - DT)
post_group_diff <- post_prob |>
  pivot_wider(names_from = group, values_from = prob) |>
  mutate(group_diff = DTF - DT) |>
  select(draws, phase, group_diff)

group_diff_summary <- post_group_diff |>
  group_by(phase) |>
  median_hdi(group_diff, .width = 0.89) |>
  mutate(across(where(is.numeric), ~round(., 2)))
group_diff_summary

# Export posterior draws ----------------------------------------------------- #
write_xlsx(
  list(
    PhaseAccuracy      = post_prob,
    DeltaFromBase      = post_delta,
    OddsRatioFromBase  = post_or,
    GroupDifference    = post_group_diff
  ),
  file = here("results", "posterior_rsvp_dts.xlsx")
)
