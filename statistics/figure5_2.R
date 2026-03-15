# ==============================================================================
# figure5_2.R
# ==============================================================================
# Purpose: Bootstrap analysis of hand angle learning curves for ST, DT, and DTF
#          groups. Computes bootstrap HDIs for group means and pairwise
#          differences per cycle. Produces a two-panel figure (group curves +
#          contrast differences) saved as Figure5-2.tif.
#
# Output:  figures/Figure5-2.tif
#
# ==============================================================================

library(HDInterval)
library(here)
library(openxlsx2)
library(patchwork)
library(purrr)
library(tidyverse)
library(zoo)

# Supporting functions ------------------------------------------------------- #

# Function: Compute bootstrap HDI of a mean
boot_hdi <- function(x, nboot = 10000, credMass = 0.89) {
  boots <- replicate(nboot, mean(sample(x, replace = TRUE), na.rm = TRUE))
  as.numeric(HDInterval::hdi(boots, credMass = credMass))
}

# Function: Compute bootstrap HDI of a mean difference
boot_diff_hdi <- function(x1, x2, nboot = 10000, credMass = 0.89) {
  boots <- replicate(nboot, {
    m1 <- mean(sample(x1, replace = TRUE), na.rm = TRUE)
    m2 <- mean(sample(x2, replace = TRUE), na.rm = TRUE)
    m1 - m2
  })
  c(mean = mean(boots), HDInterval::hdi(boots, credMass = credMass))
}

# Load data ------------------------------------------------------------------ #
df <- read_xlsx(here("data", "data_reaching.xlsx"))

df$id <- as.factor(df$id)
df$group <- factor(df$group, levels = c("ST", "DT", "DTF"))

# Compute moving average per participant
df_smooth <- df |>
  filter(cycle >= 8, cycle <= 50) |>
  arrange(id, cycle) |>
  group_by(id) |>
  mutate(
    ha_smooth = rollapply(ha, width = 5, mean, align = "center", fill = NA),
  ) |>
  ungroup()

# Group means ---------------------------------------------------------------- #
df_summary <- df_smooth |>
  group_by(group, cycle) |>
  summarise(
    mean = mean(ha_smooth, na.rm = TRUE),
    hdi = list(boot_hdi(ha_smooth))
  ) |>
  mutate(
    hdi_lb = map_dbl(hdi, ~ .x[1]),
    hdi_ub = map_dbl(hdi, ~ .x[2])
  ) |>
  select(-hdi)

# Group mean differences ----------------------------------------------------- #

# Compute pairwise group differences per cycle
pairs <- list(
  "bold(DT)[F] - bold(ST)" = c("DTF", "ST"),
  "bold(DT) - bold(ST)"    = c("DT", "ST"),
  "bold(DT)[F] - bold(DT)" = c("DTF", "DT")
)

df_diff <- map_dfr(pairs, ~ {
  df_smooth %>%
    filter(group %in% .x) %>%
    group_by(cycle) %>%
    summarise(
      group1 = .x[1],
      group2 = .x[2],
      diff_stats = list(boot_diff_hdi(ha_smooth[group == .x[1]],
                                      ha_smooth[group == .x[2]]))
    ) %>%
    mutate(mean_diff = map_dbl(diff_stats, 1),
           hdi_lb    = map_dbl(diff_stats, 2),
           hdi_ub    = map_dbl(diff_stats, 3))
}, .id = "contrast_label") %>%
  select(-diff_stats)

# Create figure -------------------------------------------------------------- #
group_colors <- c("#3182bd", "#c51b8a", "#ff8c12")

# Panel a
p1 <- ggplot(df_summary, aes(x = cycle, y = mean, color = group, fill = group)) +
  geom_smooth(method = "loess", span = 0.33, se = FALSE, linewidth = 1) +
  geom_ribbon(aes(ymin = hdi_lb, ymax = hdi_ub), alpha = 0.2, color = NA) +
  scale_color_manual(values = group_colors, labels = c("ST", "DT", expression(DT[F]))) +
  scale_fill_manual(values = group_colors, labels = c("ST", "DT", expression(DT[F]))) +
  labs(x = "cycle", y = "hand angle (°)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, -10, 0)
  )

# Find cycles where HDI > 0
df_diff_pos <- df_diff %>%
  filter(hdi_lb > 0) %>%
  mutate(y = -1.75)

# Panel b
p2 <- ggplot(df_diff, aes(x = cycle, y = mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_ribbon(aes(ymin = hdi_lb, ymax = hdi_ub), fill = "grey80", alpha = 0.6, color = NA) +
  geom_smooth(method = "loess", span = 0.33, se = FALSE, linewidth = 1, color = "grey30") +
  geom_text(data = df_diff_pos, aes(x = cycle, y = y), label = "*", size = 3, color = "black") +
  facet_wrap(~ contrast_label, ncol = 1, scales = "fixed", labeller = label_parsed) +
  coord_cartesian(xlim = c(10, 50), ylim = c(-2, 6)) +
  scale_y_continuous(breaks = c(0, 3, 6)) +
  labs(x = "cycle", y = expression(Delta~"hand angle (°)")) +
  theme_classic(base_size = 10) +
  theme(
    legend.position   = "none",
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold"),
    panel.grid.major  = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey95", linewidth = 0.2)
  )

# Combine panels
p_combined <- (p1 + p2) +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.02, 0.98)
  )
p_combined

ggsave(
  filename = here("..", "figures", "figure5_2.tif"),
  plot = p_combined, width = 6.5, height = 3.25, units = "in", dpi = 600
)
