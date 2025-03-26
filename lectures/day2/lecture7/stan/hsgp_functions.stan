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

vector hsgp_se(vector x, real sigma, real ell, vector lambdas, matrix PHI, vector z) {
  // HSGP functions
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