# ------------------------------------------------------------------------------
# rstan_diagnostics.R
# ------------------------------------------------------------------------------
# Function: rstan_diagnostics()
# Purpose:  Evaluate key computational diagnostics for an rstan model fit.
# References: Gelman et al. (2020); Baribault & Collins (2025)
# ------------------------------------------------------------------------------
# Diagnostics reported:
#   • Rhat (≤ 1.01)          → convergence across chains
#   • ESS (≥ 400)            → effective sample size for bulk and tails
#   • Divergent transitions  → geometric instability
#   • Max treedepth hits     → inefficient integration or step-size limits
#   • E-BFMI (≥ 0.2)         → energy exploration adequacy
# ------------------------------------------------------------------------------
rstan_diagnostics <- function(fit) {
  library(posterior)
  library(rstan)
  
  # Extract sampler parameters (per chain)
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  # Get max_treedepth from the model (or default to 10)
  max_td <- fit@stan_args[[1]]$control$max_treedepth %||% 10
  
  # Sampler-derived diagnostics
  diag <- list(
    # Divergences
    divergences = sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
    
    # Max treedepth hits
    treedepth = sum(sapply(
      sp,
      \(x) sum(x[, "treedepth__"] >= max_td)
    )),
    
    # Energy-BFMI
    ebfmi = sapply(sp, function(x) {
      e <- x[, "energy__"]
      var(e) / mean(diff(e)^2)
    })
  )
  
  # Posterior draw summaries for model-derived diagnostics (Rhat & ESS)
  s <- posterior::summarize_draws(posterior::as_draws(fit))
  
  # Print diagnostic summary
  cat(sprintf(
    "\nR̂ > 1.01: %d | ESS < 400: %d | Divergences: %d | Treedepth hits: %d | E-BFMI < 0.2: %d\n",
    sum(s$rhat > 1.01, na.rm = TRUE),
    sum(s$ess_bulk < 400 | s$ess_tail < 400, na.rm = TRUE),
    diag$divergences,
    diag$treedepth,
    sum(diag$ebfmi < 0.2)
  ))
  
  # Return structured results invisibly
  invisible(list(
    sampler = diag,
    summary = s
  ))
}
