functions {
  #include hsgp_functions.stan
}

data {
  int<lower=0> N;
  int<lower=0> N0;
  int<lower=0> N1;

  vector[N] t;        // time points
  array[N0] int tid0; // treatment indicator for the first segment
  array[N1] int tid1; // time points for the first segment

  vector[N] y;        // observations  
}

transformed data {
  real C = 1.5;
  int M = 20;
  
  // Boundary condition
  real<lower=0> L = C * max(abs(t));

  // Compute the eigenvalues
  vector[M] lambdas = eigenvalues(M, L);
  
  // Compute the eigenvectors
  matrix[N, M] PHI = eigenvectors(t, M, L, lambdas);
}

parameters {
  // Baseline parameters
  array[2] real alpha;

  // Observation noise
  real<lower=0> sigma_eps;

  // HSGP parameters
  array[2] real<lower=0> sigma;
  array[2] real<lower=0> ell;
  array[2] vector[M] z;
}

transformed parameters {
  vector[N] mu0 = alpha[1] + hsgp_matern52(t, sigma[1], ell[1], lambdas, PHI, z[1]);
  vector[N] mu1 = alpha[2] + hsgp_matern52(t, sigma[2], ell[2], lambdas, PHI, z[2]);
}

model {
  // Priors
  alpha ~ normal(0, 1);
  sigma_eps ~ inv_gamma(10, 10);

  // HSGP related priors
  sigma ~ inv_gamma(10, 10);
  ell ~ inv_gamma(10, 10);
  z[1] ~ normal(0, 1);
  z[2] ~ normal(0, 1);

  // Likelihood
  y[tid0] ~ normal(mu0[tid0], sigma_eps);
  y[tid1] ~ normal(mu1[tid1], sigma_eps);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] y_rep_cf; // Counterfactual predictions

  for (n in 1:N0) {
    log_lik[n] = normal_lpdf(y[n] | mu0[n], sigma_eps);
    y_rep[n] = normal_rng(mu0[n], sigma_eps);
    y_rep_cf[n] = normal_rng(mu0[n], sigma_eps);
  }

  for (n in N0+1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu1[n], sigma_eps);
    y_rep[n] = normal_rng(mu1[n], sigma_eps);
    y_rep_cf[n] = normal_rng(mu0[n], sigma_eps);
  }
}


