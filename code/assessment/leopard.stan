
functions {

	matrix tmatrix(real[] S_equilibrium, real H) {
	
		matrix[14,14] M;
		real S[14];
		
		# initialize
		for (i in 1:14)
			for (j in 1:14)
				M[i,j] <- 0.0;
				
		# apply harvest rate
		for (i in 1:14)
			S[i] <- S_equilibrium[i] * (1 - H);
		
		# cub
		M[2,1] <- S[1];

		# sub-adult
		M[3,2] <- 0.5 * S[2];    # juvenile to sub-adult female
		M[9,2] <- 0.5 * S[2];    # juvenile to sub-adult male

		# adult
		M[4,3]   <- S[3];        # sub-adult female
		M[5,4]   <- S[4];        # adult female (36-48)
		M[6,5]   <- S[5];        # adult female (48-60)
		M[7,6]   <- S[6];        # adult female (60-72)
		M[8,7]   <- S[7];        # adult female (72-84)
		M[8,8]   <- S[8];        # adult female (>84)
		M[10,9]  <- S[9];        # sub-adult male
		M[11,10] <- S[10];       # adult male (36-48)
		M[12,11] <- S[11];       # adult male (48-60)
		M[13,12] <- S[12];       # adult male (60-72)
		M[14,13] <- S[13];       # adult male (72-84)
		M[14,14] <- S[14];       # adult male (>84)
		
		return M;
	}
	
	real birth(vector N, real[] k, real h) {
	
		real B;
	
		vector[5] nf;
		vector[5] nm;
		vector[5] nc;
		real nh;
		
		nf <- N[4:8];
		nm <- N[10:14];
		for (a in 1:5)
			nc[a] <- nf[a] * k[a];
		
		nh <- sum(nf) / h;
		
		B <- sum(nc) / nh * 2 * sum(nm) * nh / (sum(nm) + nh);
				
		return B;
	}
}

data {
	int T;
	int A;
	real S[A];
	real k[5];
	real h;
	int kills[T];
}

transformed data {
}

parameters {
	real<lower=0,upper=1> H;
	real<lower=0,upper=1000> N0[A];
}

transformed parameters { 
}

model {
	vector[A] N;
	matrix[A,A] M;
	real K;
	
	for (a in 1:A) 
		N[a] <- N0[a];
	
	for (t in 2:T) {
		M <- tmatrix(S,H);
		N <- M * N;
		N[1] <- birth(N,k,h);
		K <- sum(N) * H;
		kills[t] ~ poisson(K);
		print("K:", K);
	}
}

generated quantities {
	matrix[A,A] Mout;
	
	Mout <- tmatrix(S,H);
	
}
