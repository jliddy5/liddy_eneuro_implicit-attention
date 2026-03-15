# ==============================================================================
# model_parameters.R
# ==============================================================================
# Purpose: Fits a bivariate Student-t model (via Stan) to logit-transformed
#          state-space parameters (retention A, error sensitivity b) for ST, DT,
#          and DTF groups. Exports posterior group means to:
#            statistics/results/posterior_model_learningparameters.xlsx
#          Also produces the A vs b posterior scatter + correlation density plot
#          used in Figure 7.
#
# Output:  statistics/results/posterior_model_learningparameters.xlsx
#          figures/Figure7.tif
#
# ==============================================================================

library(cowplot)
library(HDInterval)
library(here)
library(openxlsx2)
library(patchwork)
library(purrr)
library(rstan)
library(tidybayes)
library(tidyverse)

source(here("utils", "rstan_diagnostics.R"))

# Load data ------------------------------------------------------------------ #
data <- read_xlsx(here("data", "data_modeling.xlsx"))

data <- data %>%
  mutate(
    id    = as.factor(id),
    group = factor(group, levels = c("ST", "DT", "DTF")),
    Logit_A = qlogis(a),
    Logit_b = qlogis(b)
  )

# Group colors (ST, DT, DTF)
group_colors <- c("#3182bd", "#c51b8a", "#ff8c12")

# Exploratory plots ---------------------------------------------------------- #
p1 <- ggplot(data, aes(x = group, y = Logit_A)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.9) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  labs(y = "logit(A)", title = "Logit A") +
  theme(legend.position = "none")

p2 <- ggplot(data, aes(x = group, y = Logit_b)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.9) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  labs(y = "logit(b)", title = "Logit b") +
  theme(legend.position = "none")

p1 + p2

# Model ---------------------------------------------------------------------- #
Y_mat <- cbind(data$Logit_A, data$Logit_b)
data_model <- list(
  N       = nrow(data),
  G       = length(levels(data$group)),
  group   = as.integer(data$group),
  Y       = Y_mat,
  Y_mean  = apply(Y_mat, 2, mean),
  Y_sd    = apply(Y_mat, 2, sd)
)

stan_model <- stan_model(file = here("models", "model_statespaceparam.stan"))

# Fit model (seed for reproducibility: 190943950)
fit <- rstan::sampling(stan_model, data = data_model,
                       iter = 3500, warmup = 1000, chains = 4,
                       control = list(adapt_delta = 0.9, max_treedepth = 15),
                       seed = 190943950)

# Diagnostics
rstan_diagnostics(fit)

# Extract posterior samples -------------------------------------------------- #
posterior_samples <- rstan::extract(fit, permuted = TRUE)

# Retention, A --------------------------------------------------------------- #
mu_A <- data.frame(
  ST  = posterior_samples$mu_prob[,1,1],
  DT  = posterior_samples$mu_prob[,2,1],
  DTF = posterior_samples$mu_prob[,3,1]
)

mu_A_long <- mu_A %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "mu") %>%
  mutate(group = factor(group, levels = c("ST", "DT", "DTF")))

ggplot(mu_A_long, aes(x = mu, fill = group)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  labs(title = "posterior density of the mean", x = "retention, A", y = "density", fill = NULL) +
  theme_minimal()

# Summarize
mu_A %>%
  pivot_longer(everything(), names_to = "group", values_to = "mu") %>%
  group_by(group) %>%
  median_hdi(mu, .width = 0.89) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))

# Contrasts
mu_A_diff <- data.frame(
  DTF_ST = mu_A$DTF - mu_A$ST,
  DTF_DT = mu_A$DTF - mu_A$DT,
  DT_ST  = mu_A$DT  - mu_A$ST
)

mu_A_diff %>%
  pivot_longer(everything(), names_to = "contrast", values_to = "diff") %>%
  group_by(contrast) %>%
  median_hdi(diff, .width = 0.89) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))

# Error sensitivity, b ------------------------------------------------------- #
mu_b <- data.frame(
  ST  = posterior_samples$mu_prob[,1,2],
  DT  = posterior_samples$mu_prob[,2,2],
  DTF = posterior_samples$mu_prob[,3,2]
)

mu_b_long <- mu_b %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "mu") %>%
  mutate(group = factor(group, levels = c("ST", "DT", "DTF")))

ggplot(mu_b_long, aes(x = mu, fill = group)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  labs(title = "posterior density of the mean", x = "error sensitivity, b", y = "density", fill = NULL) +
  theme_minimal()

# Summarize
mu_b %>%
  pivot_longer(everything(), names_to = "group", values_to = "mu") %>%
  group_by(group) %>%
  median_hdi(mu, .width = 0.89) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))

# Directional probabilities
cat("Pr(DTF > ST) =", round(mean(mu_b$DTF > mu_b$ST), 2), "\n")
cat("Pr(DT > ST) =",  round(mean(mu_b$DT  > mu_b$ST), 2), "\n")
cat("Pr(DTF > DT) =", round(mean(mu_b$DTF > mu_b$DT), 2), "\n")
cat("Pr(ST < DT < DTF) =", round(mean(mu_b$ST < mu_b$DT & mu_b$DT < mu_b$DTF), 2), "\n")

# Contrasts
mu_b_diff <- data.frame(
  DTF_ST = mu_b$DTF - mu_b$ST,
  DTF_DT = mu_b$DTF - mu_b$DT,
  DT_ST  = mu_b$DT  - mu_b$ST
)

mu_b_diff %>%
  pivot_longer(everything(), names_to = "contrast", values_to = "diff") %>%
  group_by(contrast) %>%
  median_hdi(diff, .width = 0.89) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))

# Export posterior draws ----------------------------------------------------- #
wb <- wb_workbook()
wb <- wb_add_worksheet(wb, "Retention")
wb <- wb_add_data(wb, sheet = "Retention", x = data.frame(
  Mean_ST  = mu_A$ST,
  Mean_DT  = mu_A$DT,
  Mean_DTF = mu_A$DTF
))
wb <- wb_add_worksheet(wb, "ErrorSensitivity")
wb <- wb_add_data(wb, sheet = "ErrorSensitivity", x = data.frame(
  Mean_ST  = mu_b$ST,
  Mean_DT  = mu_b$DT,
  Mean_DTF = mu_b$DTF
))
wb_save(wb, here("results", "posterior_model_parameters.xlsx"))

# Export logit-scale group summary (used by degeneracy test) ----------------- #
logit_params <- map_dfr(1:3, function(g) {
  tibble(
    Group   = c("ST", "DT", "DTF")[g],
    Mu_A    = median(posterior_samples$mu[, g, 1]),
    Sigma_A = median(posterior_samples$sigma[, g, 1]),
    Mu_b    = median(posterior_samples$mu[, g, 2]),
    Sigma_b = median(posterior_samples$sigma[, g, 2]),
    Nu      = median(posterior_samples$nu[, g])
  )
})
wb_logit <- wb_workbook()
wb_logit <- wb_add_worksheet(wb_logit, "LogitParams")
wb_logit <- wb_add_data(wb_logit, sheet = "LogitParams", x = logit_params)
wb_save(wb_logit, here("results", "posterior_model_logitparameters.xlsx"))

# Figure 7 ------------------------------------------------------------------- #
group_labels <- c("ST", "DT", "DT[F]")
group_colors_named <- c("ST" = "#3182bd", "DT" = "#c51b8a", "DT[F]" = "#ff8c12")

mu_draws   <- posterior_samples$mu_prob
corr_draws <- posterior_samples$corr_A_b
n_draws    <- dim(mu_draws)[1]

mu_df <- map_dfr(1:3, function(g) {
  tibble(
    Retention        = mu_draws[, g, 1],
    ErrorSensitivity = mu_draws[, g, 2],
    Group            = group_labels[g]
  )
}) %>% mutate(Group = factor(Group, levels = group_labels))

main_plot <- ggplot(mu_df, aes(x = Retention, y = ErrorSensitivity, color = Group)) +
  geom_point(alpha = 0.25, size = 0.5, stroke = 0) +
  facet_wrap(~Group, labeller = label_parsed) +
  scale_color_manual(values = group_colors_named) +
  scale_x_continuous(limits = c(0.82, 0.97), breaks = seq(0.85, 0.95, 0.05)) +
  scale_y_continuous(limits = c(0.01, 0.05), breaks = seq(0.01, 0.05, 0.01)) +
  labs(
    x = expression("retention, " * italic(A)),
    y = expression("error sensitivity, " * italic(b))
  ) +
  theme_minimal(base_family = "sans", base_size = 8) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    axis.text        = element_text(color = "black", size = 5),
    axis.title       = element_text(size = 10),
    strip.text       = element_text(size = 10, face = "bold"),
    axis.line        = element_line(color = "black", linewidth = 0.15),
    axis.ticks       = element_line(color = "black", linewidth = 0.15),
    axis.ticks.length = unit(0.075, "cm"),
    legend.position  = "none"
  )

corr_df <- tibble(
  Correlation = as.numeric(corr_draws),
  Group       = rep(group_labels, each = n_draws)
)

r_plot_list <- map(group_labels, function(g) {
  ggplot(filter(corr_df, Group == g), aes(x = Correlation)) +
    geom_density(fill = group_colors_named[g], color = "black", alpha = 0.6, linewidth = 0.2) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(x = expression(italic(r)), y = NULL) +
    theme_minimal(base_size = 6) +
    theme(
      panel.grid       = element_blank(),
      axis.title.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      strip.text       = element_blank(),
      axis.title.x     = element_text(size = 5),
      axis.text.x      = element_text(size = 5),
      axis.line.x      = element_line(color = "black", linewidth = 0.1),
      plot.background  = element_blank(),
      panel.background = element_blank(),
      panel.border     = element_blank(),
      legend.position  = "none"
    )
})
names(r_plot_list) <- group_labels

final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(r_plot_list[["ST"]],    x = 0.08, y = 0.57, width = 0.15, height = 0.25) +
  draw_plot(r_plot_list[["DT"]],    x = 0.38, y = 0.15, width = 0.15, height = 0.25) +
  draw_plot(r_plot_list[["DT[F]"]], x = 0.69, y = 0.15, width = 0.15, height = 0.25)

final_plot

ggsave(
  filename = here("..", "figures", "figure7.tif"),
  plot = final_plot, width = 6.5, height = 2.25, units = "in", dpi = 600, bg = "white"
)
