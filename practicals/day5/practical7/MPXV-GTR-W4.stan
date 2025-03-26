functions{
	real ctmc_scale_lpdf(real rate, vector blens){
		real total_tree_time = sum(blens);
		real log_normalization = 0.5 * log(total_tree_time) - 0.5723649263381958; //lgamma(0.5);
		real log_like = log_normalization - 0.5 * log(rate) - rate * total_tree_time;
		return log_like;
	}


	// transform node heights to proportion, except for the root
	vector transform(real[] p, real rootHeight, int[,] map, real[] lowers){
		int S = size(p)+2;
		int nodeCount = S*2-1;
		
		vector[nodeCount] heights;
		int j = 1;
		
		heights[map[1,1]] = rootHeight;
		for( i in 2:nodeCount ){
			// internal node: transform
			if(map[i,1] > S){
				heights[map[i,1]] = lowers[map[i,1]] + (heights[map[i,2]] - lowers[map[i,1]])*p[map[i,1]-S];
				j += 1;
			}
			else{
				heights[map[i,1]] = lowers[map[i,1]];
			}
		}
		return heights;
	}
	

	real oneOnX_lpdf(real x){
		return -log(x);
	}


	real constant_coalescent_lpdf(vector heights, real popSize, int[,] map){
		int nodeCount = rows(heights);
		int S = (nodeCount+1)/2;
		int intervalCount = nodeCount - 1;

		int indices[nodeCount] = sort_indices_asc(heights);
		vector[nodeCount] events = append_row(rep_vector(1, S), rep_vector(-1, S-1));
		vector[intervalCount] lineageCount = cumulative_sum(events[indices])[:intervalCount];
		vector[intervalCount] intervals = heights[indices][2:] - heights[indices][:intervalCount];

		return -sum(intervals .* ((lineageCount .* (lineageCount-1.0))/2.0)) / popSize - (S-1)*log(popSize);
	}
	

	matrix[] calculate_hky_p_matrices(vector freqs, real kappa, vector blens, vector rs){
		int C = rows(rs);
		int bcount = rows(blens);
		matrix[4,4] pmats[bcount*C]; // probability matrices

		matrix[4,4] Q; // rate matrix
		matrix[4,4] P2 = diag_matrix(sqrt(freqs));        // diagonal sqrt frequencies
		matrix[4,4] P2inv = diag_matrix(1.0 ./ sqrt(freqs)); // diagonal inverse sqrt frequencies
		matrix[4,4] A; // another symmetric matrix
		vector[4] eigenvalues;
		matrix[4,4] eigenvectors;
		matrix[4,4] m1;
		matrix[4,4] m2;
		// symmetric rate matrix
		matrix[4,4] R = [[0.0, 1.0, kappa, 1.0],
						 [1.0, 0.0, 1.0, kappa],
						 [kappa, 1.0, 0.0, 1.0],
						 [1.0, kappa, 1.0, 0.0]];
		real s = 0.0;
		int index = 1;

		Q = R * diag_matrix(freqs);
		for (i in 1:4) {
			Q[i,i] = 0.0;
			Q[i,i] = -sum(Q[i,1:4]);
			s -= Q[i,i] * freqs[i];
		}
		Q /= s;

		A = P2 * Q * P2inv;

		eigenvalues = eigenvalues_sym(A);
		eigenvectors = eigenvectors_sym(A);

		m1 = P2inv * eigenvectors;
		m2 = eigenvectors' * P2;

		for(c in 1:C){
			for( b in 1:bcount ){
				pmats[index] = m1 * diag_matrix(exp(eigenvalues*blens[b]*rs[c])) * m2;
				index += 1;
			}
		}

		return pmats;
	}

}

data{
	int <lower=0> L;                      // alignment length
	int <lower=0> S;                      // number of tips
	vector<lower=0,upper=1>[4] tipdata[S,L]; // alignment as partials
	int <lower=0,upper=2*S> peel[S-1,3];  // list of nodes for peeling
	real weights[L];
	int map[2*S-1,2];                     // list of node in preorder [node,parent]
	int C;
	real lower_root;
	real lowers[2*S-1]; // list of lower bounds for each internal node (for reparametrization)
	vector<lower=0>[4] frequencies_alpha; // parameters of the prior on frequencies
}

transformed data{
	int bcount = 2*S-2; // number of branches
	int nodeCount = 2*S-1; // number of nodes
	int pCount = S-2; // number of proportions
}

parameters{
	real<lower=0.1> wshape;
	real <lower=0,upper=1> props[pCount]; // proportions
	real <lower=0> rate;
	real <lower=lower_root> height; // root height
	real <lower=0> theta;
	real<lower=0> kappa;
	simplex[4] freqs;
}

transformed parameters{
	vector[C] ps = rep_vector(1.0/C, C);
	vector[C] rs;
	vector <lower=0> [2*S-1] heights;
	vector<lower=0> [bcount] blensUnscaled; // unscaled branch lengths

	
		{
			real m = 0;
			for(i in 1:C){
				rs[i] = pow(-log(1.0 - (2.0*(i-1)+1.0)/(2.0*C)), 1.0/wshape);
			}
			m = sum(rs)/C;
			for(i in 1:C){
				rs[i] /= m;		
			}
		}

	heights = transform(props, height, map, lowers);
	blensUnscaled[map[2:,1]] = heights[map[2:,2]] - heights[map[2:,1]];
}

model{
	real probs[C];
	vector[4] partials[C,S-1,L];   // partial probabilities for the S tips and S-1 internal nodes
	matrix[4,4] pmats[bcount*C]; // finite-time transition matrices for each branch
	vector [bcount] blens; // branch lengths
	int p1;
	int p2;
	int p3;

	wshape ~ exponential(1.0);
	rate ~ ctmc_scale(blensUnscaled);
	theta ~ oneOnX();
	heights ~ constant_coalescent(theta, map);
	kappa ~ lognormal(1.0,1.25);
	freqs ~ dirichlet(frequencies_alpha);

	blens = rate*blensUnscaled;
	pmats = calculate_hky_p_matrices(freqs, kappa, blens, rs);
	// calculate tree likelihood

    	for(c in 1:C){
		for( n in 1:(S-1) ) {
			p1 = peel[n,1];
			p2 = peel[n,2];
			p3 = peel[n,3] - S;
			if(peel[n,1] <= S && peel[n,2] <= S){
				for( i in 1:L ) {
					partials[c,p3,i] = (pmats[p1+(c-1)*bcount]*tipdata[p1,i]) .* (pmats[p2+(c-1)*bcount]*tipdata[p2,i]);
				}
			}
			else if(peel[n,1] <= S){
				for( i in 1:L ) {
					partials[c,p3,i] = (pmats[p1+(c-1)*bcount]*tipdata[p1,i]) .* (pmats[p2+(c-1)*bcount]*partials[c,p2-S,i]);
				}
			}
			else if(peel[n,2] <= S){
				for( i in 1:L ) {
					partials[c,p3,i] = (pmats[p1+(c-1)*bcount]*partials[c,p1-S,i]) .* (pmats[p2+(c-1)*bcount]*tipdata[p2,i]);
				}
			}
			else{
				for( i in 1:L ) {
					partials[c,p3,i] = (pmats[p1+(c-1)*bcount]*partials[c,p1-S,i]) .* (pmats[p2+(c-1)*bcount]*partials[c,p2-S,i]);
				}		
			}
		}
	}
	p3 = peel[S-1,3] - S;
	for( i in 1:L ) {
		for(c in 1:C){
			probs[c] = ps[c] * sum(partials[c,p3,i] .* freqs);
		}
		target += log(sum(probs))*weights[i];
	}
    
	
	// add log det jacobian
	for( i in 2:nodeCount ){
		// skip leaves
		if(map[i,1] > S ){
			target += log(heights[map[i,2]] - lowers[map[i,1]]);
		}
	}

}

