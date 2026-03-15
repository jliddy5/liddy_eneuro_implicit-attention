data {
  int<lower=1> N;                  // Number of observations
  int<lower=1> G;                  // Number of groups (1, 2, ..., G)
  int<lower=1, upper=G> group[N];  // Group ID for each observation
  vector[2] Y[N];                  // [logit(A), logit(b)] per observation
  vector[2] Y_mean;                // Full-sample logit mean per group
  vector[2] Y_sd;                  // Full-sample logit SD per group
}

parameters {
  vector[2] mu[G];                 // Group-specific means
  vector<lower=0>[2] sigma[G];     // Group-specific SDs
  cholesky_factor_corr[2] L[G];    // Group-specific correlations
  real<lower=0> nu_minus_one[G];   // Group-specific DOF
}

transformed parameters {
  cov_matrix[2] Sigma[G];
  real<lower=1> nu[G];
  
  for (g in 1:G) {
    nu[g] = nu_minus_one[g] + 1;
    Sigma[g] = diag_pre_multiply(sigma[g], L[g]) * diag_pre_multiply(sigma[g], L[g])';
  }
}

model {
  // Priors
  for (g in 1:G) {
    nu_minus_one[g] ~ exponential(0.03448);
    mu[g] ~ normal(Y_mean, Y_sd);

    for (j in 1:2) {
      sigma[g, j]^2 ~ inv_gamma(2, square(Y_sd[j]));
    }

    L[g] ~ lkj_corr_cholesky(2);
  }

  // Likelihood
  for (n in 1:N) {
    Y[n] ~ multi_student_t(nu[group[n]], mu[group[n]], Sigma[group[n]]);
  }
}

generated quantities {
  vector[2] mu_prob[G];
  real corr_A_b[G];

  for (g in 1:G) {
    // Inverse logit to (0,1)
    mu_prob[g, 1] = inv_logit(mu[g, 1]); // Retention
    mu_prob[g, 2] = inv_logit(mu[g, 2]); // Error sensitivity
    
    // Correlation
    corr_A_b[g] = multiply_lower_tri_self_transpose(L[g])[1, 2];
  }
}
