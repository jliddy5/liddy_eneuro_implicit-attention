# ------------------------------------------------------------------------------
# brms_diagnostics.R
# ------------------------------------------------------------------------------
# Function: brms_diagnostics()
# Purpose:  Evaluate key computational diagnostics for a brms model fit.
# References: Gelman et al. (2020); Baribault & Collins (2025)
# ------------------------------------------------------------------------------
# Diagnostics reported:
#   • Rhat (≤ 1.01)          → convergence across chains
#   • ESS (≥ 400)            → effective sample size for bulk and tails
#   • Divergent transitions  → geometric instability
#   • Max treedepth hits     → inefficient integration or step-size limits
#   • E-BFMI (≥ 0.2)         → energy exploration adequacy
# ------------------------------------------------------------------------------

brms_diagnostics <- function(fit) {
  library(posterior)
  library(rstan)
  
  # Extract sampler parameters (per chain)
  sp <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
  
  # Sampler-derived diagnostics
  diag <- list(
    # Divergences
    divergences = sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
    
    # Max treedepth hits
    treedepth = sum(sapply(
      sp,
      \(x) sum(x[, "treedepth__"] >= attr(sp[[1]], "max_depth") %||% 10)
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
