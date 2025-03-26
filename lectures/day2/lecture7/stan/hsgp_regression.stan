functions {
  // Spectral density functions
	vector spd_se(vector omega, real sigma, real ell) {
		return sigma^2 * sqrt(2 * pi()) * ell * exp(-0.5 * ell^2 * omega .* omega);
	}

	vector spd_matern32(vector omega, real sigma, real ell) {
		return sigma^2 * 12 * sqrt(3) / ell^3 * (3 / ell^2 + omega .* omega).^(-2);
	}

	vector spd_matern52(vector omega, real sigma, real ell) {
		return sigma^2 * 400 * sqrt(5) / (3 * ell^5) * (5 / ell^2 + omega .* omega).^(-3);
	}

  // Eigenvalues and eigenvectors
  vector eigenvalues(int M, real L) {
		vector[M] lambda;
		for (m in 1:M) {
			lambda[m] = (m * pi() / (2 * L))^2;
		}
		return lambda;
	}

	matrix eigenvectors(vector x, int M, real L, vector lambda) {
		int N = rows(x);
		matrix[N, M] PHI;
		for (m in 1:M) {
			PHI[,m] = sqrt(1 / L) * sin(sqrt(lambda[m]) * (x + L));
		}
		return PHI;
	}

  // HSGP functions
  vector hsgp_se(vector x, real sigma, real ell, vector lambdas, matrix PHI, vector z) {
		int N = rows(x);
		int M = cols(PHI);
		vector[N] f;
		matrix[M, M] Delta;

		// Spectral densities evaluated at the square root of the eigenvalues
		vector[M] spds = spd_se(sqrt(lambdas), sigma, ell);

		// Construct the diagonal matrix Delta
		Delta = diag_matrix(sqrt(spds));

		// Compute the HSGP sample
		f = PHI * Delta * z;

		return f;
	}

  vector hsgp_matern32(vector x, real sigma, real ell, vector lambdas, matrix PHI, vector z) {
		int N = rows(x);
		int M = cols(PHI);
		vector[N] f;
		matrix[M, M] Delta;

		// Spectral densities evaluated at the square root of the eigenvalues
		vector[M] spds = spd_matern32(sqrt(lambdas), sigma, ell);

		// Construct the diagonal matrix Delta
		Delta = diag_matrix(sqrt(spds));

		// Compute the HSGP sample
		f = PHI * Delta * z;

		return f;
	}

  vector hsgp_matern52(vector x, real sigma, real ell, vector lambdas, matrix PHI, vector z) {
		int N = rows(x);
		int M = cols(PHI);
		vector[N] f;
		matrix[M, M] Delta;

		// Spectral densities evaluated at the square root of the eigenvalues
		vector[M] spds = spd_matern52(sqrt(lambdas), sigma, ell);

		// Construct the diagonal matrix Delta
		Delta = diag_matrix(sqrt(spds));

		// Compute the HSGP sample
		f = PHI * Delta * z;

		return f;
	}
}

data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
  real<lower=0> C;
  int<lower=1> M;
}

transformed data {
  // Boundary condition
  real<lower=0> L = C * max(abs(x));
  
  // Compute the eigenvalues
  vector[M] lambdas = eigenvalues(M, L);

  // Compute the eigenvectors
  matrix[N, M] PHI = eigenvectors(x, M, L, lambdas);
}

parameters {
  real alpha;
  real<lower=0> sigma_eps;

  real<lower=0> sigma;
  real<lower=0> ell;
  vector[M] z;
}

transformed parameters {
  vector[N] f = hsgp_se(x, sigma, ell, lambdas, PHI, z);
  vector[N] mu = alpha + f;
}

model {
  // Priors
  alpha ~ normal(0, 1);
  sigma_eps ~ inv_gamma(5, 5);

  // GP related priors
  sigma ~ inv_gamma(5, 5);
  ell ~ inv_gamma(5, 5);
  z ~ normal(0, 1);

  // Likelihood
  y ~ normal(mu, sigma_eps);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_eps);
    y_rep[n] = normal_rng(mu[n], sigma_eps);
  }
}