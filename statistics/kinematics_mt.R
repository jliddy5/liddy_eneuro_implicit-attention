# ==============================================================================
# kinematics_mt.R
# ==============================================================================
# Purpose: Bayesian lognormal regression of movement time over learning cycles
#          for ST, DT, and DTF groups. Fits two models (with/without group) and
#          compares via LOO. Exports posterior group-level fits as SVG panels
#          used to assemble Figure 5-1 in Illustrator.
#
# Note:    Figure 5-1 was assembled manually in Illustrator from the SVG outputs
#          of kinematics_rt.R, kinematics_mt.R, and kinematics_tad.R.
#
# Output:  Plots displayed in R. Uncomment ggsave() calls to export SVG panels.
#
# ==============================================================================

library(bayesplot)
library(brms)
library(emmeans)
library(here)
library(openxlsx2)
library(svglite)
library(tidybayes)
library(tidyverse)

# Load data ------------------------------------------------------------------ #
df <- read_xlsx(here("data", "data_reaching.xlsx"))

df$id <- as.factor(df$id)
df$group <- as.factor(df$group)

df <- df %>%
  filter(cycle %in% 11:50) %>%
  mutate(
    cycle_01 = (cycle - min(cycle)) / (max(cycle) - min(cycle)),
    group = fct_relevel(group, "ST")
  )

# Group colors (ST, DT, DTF)
group_colors <- c("#3182bd", "#c51b8a", "#ff8c12")

# Data ----------------------------------------------------------------------- #

df_summ <- df %>%
  group_by(group, cycle) %>%
  summarise(
    mean_mt = mean(mt),
    se_mt   = sd(mt) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_mt - 2 * se_mt,
    upper = mean_mt + 2 * se_mt
  )

theme_set(theme_bw())
p1 <- ggplot(df, aes(x = cycle, y = mt)) +
  geom_point(
    aes(color = group),
    position = position_jitter(width = .2, height = 0, seed = 1),
    alpha = 0.75, size = 0.5, shape = 16
  ) +
  geom_line(
    data = df_summ, aes(x = cycle, y = mean_mt, color = group),
    linewidth = 0.5, inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = df_summ,
    aes(x = cycle, ymin = lower, ymax = upper, fill = group),
    alpha = 0.25, inherit.aes = FALSE
  ) +
  scale_color_manual(values = group_colors, labels = c("ST", "DT", expression(DT[F]))) +
  scale_fill_manual(values = group_colors, labels = c("ST", "DT", expression(DT[F]))) +
  scale_x_continuous(name = "cycle", breaks = c(11,20,30,40,50)) +
  scale_y_continuous(name = "movement time (ms)", limits = c(100, 400)) +
  labs(title = "data: mean ± 2 SE") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )
p1

# ggsave(here("figures", "MT_Data.svg"), p1, width = 5, height = 2.5, units = "in", device = "svg")

# Model ---------------------------------------------------------------------- #

# Seed (sample.int(1e6, 1))
seed <- 973424

# Priors
priors_mt <- c(
  prior(normal(5.3, 0.25), class = "Intercept"),  # MT ~ 200 [95, 425] ms
  prior(normal(0, 0.1),    class = "b"),           # ±10% effects
  prior(exponential(8),    class = "sigma"),       # RE SDs
  prior(exponential(8),    class = "sd")           # residual SD
)

# Prior predictive check
fit_mt_priorpc <- brm(
  mt ~ cycle_01 * group + (1 + cycle_01 | id),
  data    = df,
  family  = lognormal(),
  prior   = priors_mt,
  sample_prior = "only",
  iter    = 2000, warmup = 1000, chains = 5,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  backend = "rstan",
  seed = seed
)
pp_check(fit_mt_priorpc, ndraws = 500) + scale_x_log10()
pp_check(fit_mt_priorpc, type = "stat", stat = "sd",  ndraws = 500)
pp_check(fit_mt_priorpc, type = "stat_grouped", group = "cycle_01", stat = "mean", ndraws = 500)

# Model 1: Separate groups
fit_mt_group <- brm(
  mt ~ cycle_01 * group + (1 + cycle_01 | id),
  data    = df,
  family  = lognormal(),
  prior   = priors_mt,
  iter    = 3500, warmup = 1000, chains = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  backend = "rstan",
  seed = seed
)
pp_check(fit_mt_group, ndraws = 1000) + scale_x_log10()
pp_check(fit_mt_group, type = "stat", stat = "sd",  ndraws = 1000)
pp_check(fit_mt_group, type = "stat_grouped", group = "cycle_01", stat = "mean", ndraws = 1000)
bayes_R2(fit_mt_group, summary = TRUE, probs = c(0.055, 0.945))

# Model 2: Combine groups
fit_mt <- brm(
  mt ~ cycle_01 + (1 + cycle_01 | id),
  data    = df,
  family  = lognormal(),
  prior   = priors_mt,
  iter    = 3500, warmup = 1000, chains = 5,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  backend = "rstan",
  seed = seed
)
pp_check(fit_mt, ndraws = 1000) + scale_x_log10()
pp_check(fit_mt, type = "stat", stat = "sd",  ndraws = 1000)
pp_check(fit_mt, type = "stat_grouped", group = "cycle_01", stat = "mean", ndraws = 1000)
bayes_R2(fit_mt, summary = TRUE, probs = c(0.055, 0.945))

em <- emmeans(fit_mt, ~ 1, type = "response", tran = "log")
summary(em, type = "response", level = 0.89)

loo_group   <- loo(fit_mt_group)
loo_nogroup <- loo(fit_mt)
loo_compare(loo_group, loo_nogroup)

# Posterior predicted slopes for each person --------------------------------- #
newdat_id_cycle <- df %>%
  distinct(id, group, cycle, cycle_01) %>%
  arrange(id, cycle)

epred_id <- newdat_id_cycle %>%
  add_epred_draws(fit_mt_group, ndraws = 1000, re_formula = NULL) %>%
  group_by(id, group, cycle) %>%
  median_hdi(.epred, .width = 0.89)

ggplot() +
  geom_point(
    data = df,
    aes(x = cycle, y = mt, color = group),
    alpha = 0.55, size = 0.8, shape = 16
  ) +
  geom_ribbon(
    data = epred_id,
    aes(x = cycle, ymin = .lower, ymax = .upper, fill = group),
    alpha = 0.18, inherit.aes = FALSE
  ) +
  geom_line(
    data = epred_id,
    aes(x = cycle, y = .epred, color = group),
    linewidth = 0.8, inherit.aes = FALSE
  ) +
  facet_wrap(~ id, scales = "free_y") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(y = "movement time (ms)", x = "cycle", title = "posterior: individual-level fits") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

# Posterior predicted slopes for each group ---------------------------------- #
newdat_cycle <- df %>%
  distinct(group, cycle, cycle_01) %>%
  arrange(group, cycle)

epred_group <- newdat_cycle %>%
  add_epred_draws(fit_mt_group, ndraws = 1000, re_formula = NA) %>%
  group_by(group, cycle) %>%
  median_hdi(.epred, .width = 0.89)

p2 <- ggplot() +
  geom_point(
    data = df,
    aes(x = cycle, y = mt, color = group),
    position = position_jitter(width = .2, height = 0, seed = 1),
    alpha = 0.75, size = 0.5, shape = 16
  ) +
  geom_line(
    data = epred_group,
    aes(x = cycle, y = .epred, color = group),
    linewidth = 0.5,
  ) +
  geom_ribbon(
    data = epred_group,
    aes(x = cycle, ymin = .lower, ymax = .upper, fill = group),
    alpha = 0.25,
  ) +
  scale_color_manual(values = group_colors, labels = c("ST", "DT", expression(DT[F]))) +
  scale_fill_manual(values = group_colors, labels = c("ST", "DT", expression(DT[F]))) +
  scale_x_continuous(name = "cycle", breaks = c(11,20,30,40,50)) +
  scale_y_continuous(name = "movement time (ms)", limits = c(100, 400)) +
  labs(title = "posterior: median ± 89% HDI") +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )
p2

# ggsave(here("figures", "MT_Posterior.svg"), p2, width = 5, height = 2.5, units = "in", device = "svg")

# Extract posterior summaries ------------------------------------------------ #
newdat_edges <- expand.grid(
  cycle_01 = c(0, 1),
  group = factor(c("ST","DT","DTF"), levels = c("ST","DT","DTF"))
)

pred_edges <- posterior_epred(fit_mt_group, newdata = newdat_edges, re_formula = NA)

idx0 <- newdat_edges$cycle_01 == 0
intercepts <- apply(pred_edges[, idx0], 2, quantile, probs = c(.5, .055, .945))

idx1 <- newdat_edges$cycle_01 == 1
changes <- pred_edges[, idx1] - pred_edges[, idx0]
changes_summ <- apply(changes, 2, quantile, probs = c(.5, .055, .945))

mt_summary <- tibble(
  group = levels(newdat_edges$group),
  intercept_med = intercepts[1,],
  intercept_lb  = intercepts[2,],
  intercept_ub  = intercepts[3,],
  change_med    = changes_summ[1,],
  change_lb     = changes_summ[2,],
  change_ub     = changes_summ[3,]
)
mt_summary
