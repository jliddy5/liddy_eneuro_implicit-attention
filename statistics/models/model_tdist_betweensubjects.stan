data {
  int<lower=1> N;                   // Number of observations
  int<lower=1> G;                   // Number of groups
  vector[N] y;                      // Observations
  int<lower=1, upper=G> group[N];   // Group identifiers
  int<lower=0> n_per_group[G];      // Number of observations per group
  real meanY;                       // Full-sample mean
  real sdY;                         // Full-sample SD
}

transformed data {
  // Prior hyperparameters
  real nu_lambda;
  real sigma_beta;

  nu_lambda = 0.03448;
  sigma_beta = square(sdY);
}

parameters {
  real<lower=0> nuMinusOne;     // Degrees of freedom (shared). Ensures nu > 1.
  vector[G] mu;                 // Group means (separate)
  vector<lower=0>[G] sigma;     // Group SDs (separate)
}

transformed parameters {
  real<lower=0> nu;
  nu = nuMinusOne + 1;
}

model {
  // Priors
  nuMinusOne ~ exponential(nu_lambda);   // Degrees of freedom (shared)
  mu ~ normal(meanY, 10);                // Group means (separate)
  sigma ~ inv_gamma(2, sigma_beta);      // Group SDs (separate)

  // Likelihood
  for (i in 1:N) {
    y[i] ~ student_t(nu, mu[group[i]], sigma[group[i]]);
  }
}

generated quantities {
  // Simulated observations and group labels
  array[N] int group_sim;
  vector[N] y_sim;

  int counter = 1;
  for (g in 1:G) {
    for (n in 1:n_per_group[g]) {
      group_sim[counter] = g;
      y_sim[counter] = student_t_rng(nu, mu[g], sigma[g]);
      counter += 1;
    }
  }
}
