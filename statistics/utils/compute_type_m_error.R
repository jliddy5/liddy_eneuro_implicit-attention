# ------------------------------------------------------------------------------
# compute_type_m_error.R
# ------------------------------------------------------------------------------
# Compute Type M (magnitude/exaggeration) error for a between-group contrast.
#
# Unlike the traditional approach (Gelman & Carlin, 2014), which assumes a
# known true effect size, this implementation marginalizes over posterior
# uncertainty. For each posterior draw, the sampled group means are treated
# as "true" parameters, and a hypothetical replication is simulated. Type M
# is the ratio of the simulated effect magnitude to the true magnitude.
#
# NOTE: Requires posterior_samples from model_tdist_betweensubjects.stan, which
#   generates y_sim and group_sim in its generated quantities block.
#
# Arguments:
#   posterior_samples  Extracted Stan samples (must include mu, y_sim, group_sim)
#   group1, group2     Group indices for the contrast (default: 1 vs 2)
#   threshold          Minimum |true_diff| to include (default: 0, i.e., all draws)
#
# Returns:
#   List with mean_type_m, median_type_m, type_m_HDI, and type_m_ratios
# ------------------------------------------------------------------------------
compute_type_m_error <- function(posterior_samples, group1 = 1, group2 = 2, threshold = 0) {
  
  n_draws <- dim(posterior_samples$mu)[1]
  type_m_ratios <- rep(NA, n_draws)
  
  for (d in 1:n_draws) {
    # True mean difference for this posterior draw
    true_diff <- posterior_samples$mu[d, group1] - posterior_samples$mu[d, group2]
    
    # Simulated study
    y_sim_d <- posterior_samples$y_sim[d, ]
    group_sim_d <- posterior_samples$group_sim[d, ]
    
    # Simulated group means
    sim_mu1 <- mean(y_sim_d[group_sim_d == group1])
    sim_mu2 <- mean(y_sim_d[group_sim_d == group2])
    
    # Simulated mean difference
    sim_diff <- sim_mu1 - sim_mu2
    
    # Compute Type M ratio only if true_diff is above threshold
    if (abs(true_diff) >= threshold) {
      type_m_ratios[d] <- abs(sim_diff) / abs(true_diff)
    } else {
      type_m_ratios[d] <- NA
    }
  }
  
  # Remove NAs before summarizing
  valid_ratios <- type_m_ratios[!is.na(type_m_ratios)]
  
  list(
    mean_type_m = mean(valid_ratios),
    median_type_m = median(valid_ratios),
    type_m_HDI = quantile(valid_ratios, probs = c(0.055, 0.945)),
    type_m_ratios = valid_ratios
  )
}
