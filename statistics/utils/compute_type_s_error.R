# ------------------------------------------------------------------------------
# compute_type_s_error.R
# ------------------------------------------------------------------------------
# Compute Type S (sign) error rate for a between-group contrast.
#
# Unlike the traditional approach (Gelman & Carlin, 2014), which assumes a
# known true effect size, this implementation marginalizes over posterior
# uncertainty. For each posterior draw, the sampled group means are treated
# as "true" parameters, and a hypothetical replication is simulated. Type S
# error is scored as 1 if the simulated effect has the opposite sign.
#
# NOTE: Requires posterior_samples from model_tdist_betweensubjects.stan, which
#   generates y_sim and group_sim in its generated quantities block.
#
# Arguments:
#   posterior_samples  Extracted Stan samples (must include mu, y_sim, group_sim)
#   group1, group2     Group indices for the contrast (default: 1 vs 2)
#
# Returns:
#   Scalar: estimated Type S error rate (proportion of sign disagreements)
# ------------------------------------------------------------------------------
compute_type_s_error <- function(posterior_samples, group1 = 1, group2 = 2) {
  
  n_draws <- dim(posterior_samples$mu)[1]
  type_s_errors <- rep(NA, n_draws)
  
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
    
    # Check if sign flipped (safely ignore cases where true_diff == 0)
    if (true_diff != 0) {
      type_s_errors[d] <- sign(true_diff) != sign(sim_diff)
    } else {
      type_s_errors[d] <- NA
    }
  }
  
  # Return estimated Type S error rate (excluding NAs if any)
  mean(type_s_errors, na.rm = TRUE)
}
