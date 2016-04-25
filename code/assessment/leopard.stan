
# to do:
# make harvest rate time variant
# remove estimation of h
# model is not fit to first data point:
# fit N0 to numbers[1], kills[1] etc. 
# add additional H parameter (for unaccounted mortality)
# to be used as a sensitivity

functions {

	matrix tmatrix(real[] S_equilibrium, real[] selectivity, real H) {
	
		matrix[6,6] M;
		real S[6];
		
		# initialize
		for (i in 1:6)
			for (j in 1:6)
				M[i,j] <- 0.0;
				
		# apply harvest rate
		for (a in 1:6)
			S[a] <- S_equilibrium[a] * (1 - H * selectivity[a]);
		
		# cub
		M[2,1] <- S[1];

		# sub-adult
		M[3,2] <- 0.5 * S[2];    # juvenile to sub-adult female
		M[5,2] <- 0.5 * S[2];    # juvenile to sub-adult male

		# adult
		M[4,3]   <- S[3];        # sub-adult female
		M[4,4]   <- S[4];        # adult female
		M[6,5]   <- S[5];        # sub-adult male
		M[6,6]   <- S[6];        # adult male
		
		return M;
	}
	
	real birth(vector N, real k, real h) {
	
		real B;
	
		real nf;
		real nm;
		real nc;
		real nh;
		
		nf <- N[4];
		nm <- N[6];
		nc <- nf * k;
		nh <- nf / h;
		
		B <- nc / nh * 2 * nm * nh / (nm + nh);
				
		return B;
	}
}

data {
	// dimensions
	int T;
	int A;
	// fixed input parameters
	real S[A];
	real k;
	// data
	int kills[T];
	real density[T];
	real numbers[T];
	vector[4] proportions[T];
}

transformed data {
}

parameters {
	real h;
	real logq;
	real<lower=0,upper=0.99> H;
	real<lower=0,upper=1000> N0[A];
	real<lower=0,upper=1.0> selectivity[A];
}

transformed parameters { 
	real q;
	q <- exp(logq);
}

model {
	vector[A] N;
	matrix[A,A] M;
	vector[4] theta;
	real lambda;
	
	// initial numbers
	for (a in 1:A) 
		N[a] <- N0[a];
	
	// projection
	for (t in 2:T) {
	
		// create transition matrix
		M <- tmatrix(S,selectivity,H);
		// project forward one step
		N <- M * N;
		// add cubs
		N[1] <- birth(N,k,h);

		// fit to kills
		if (kills[t] > 0) {
			lambda <- 0;
			for (a in 1:A)
				lambda <- lambda + N[a] * H * selectivity[a];
			kills[t] ~ poisson(lambda);
		}
		// fit to numbers estimates
		if (numbers[t] > 0)
			numbers[t] ~ lognormal(log(sum(N)), 0.1);
		// fit to density estimates
		if (density[t] > 0)
			density[t] ~ lognormal(log(q * sum(N)), 0.1);
		// fit to proportions
		theta[1] <- N[3] / (N[3] + N[4] + N[5] + N[6]);
		theta[2] <- N[4] / (N[3] + N[4] + N[5] + N[6]);
		theta[3] <- N[5] / (N[3] + N[4] + N[5] + N[6]);
		theta[4] <- N[6] / (N[3] + N[4] + N[5] + N[6]);
		if (sum(proportions[t]) > 0) {
			theta ~ normal(proportions[t], 0.1);
		}
	}
	
	// priors
	h ~ uniform(1,5);
	logq ~ uniform(-10.0,0.0);
	H ~ uniform(0,1);
	N0 ~ uniform(0,1000);
	selectivity ~ uniform(0,1);
	
	// Jacobian
	//increment_log_prob(-logq);
}

generated quantities {
	
	vector[A] Nhat;
	matrix[A,A] Mhat;
	
	real killsHat[T];
	real densityHat[T];
	real numbersHat[T];
	vector[4] thetaHat[T];
	real lambda;

	// initial numbers
	for (a in 1:A) 
		Nhat[a] <- N0[a];
		
	// initial kills
	lambda <- 0;
	for (a in 1:A)
		lambda <- lambda + Nhat[a] * H * selectivity[a];
	killsHat[1] <- lambda;
	// initial numbers
	numbersHat[1] <- sum(Nhat);
	// initial density
	densityHat[1] <- q * sum(Nhat);
	// initial proportions
	thetaHat[1][1] <- Nhat[3] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
	thetaHat[1][2] <- Nhat[4] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
	thetaHat[1][3] <- Nhat[5] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
	thetaHat[1][4] <- Nhat[6] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
	
	// projection
	for (t in 2:T) {
	
		// create transition matrix
		Mhat <- tmatrix(S,selectivity,H);
		// project forward one step
		Nhat <- Mhat * Nhat;
		// add cubs
		Nhat[1] <- birth(Nhat,k,h);

		// predict kills
		lambda <- 0;
		for (a in 1:A)
			lambda <- lambda + Nhat[a] * H * selectivity[a];
		killsHat[t] <- lambda;
		// predict numbers
		numbersHat[t] <- sum(Nhat);
		// predict density
		densityHat[t] <- q * sum(Nhat);
		// predict proportions
		thetaHat[t][1] <- Nhat[3] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
		thetaHat[t][2] <- Nhat[4] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
		thetaHat[t][3] <- Nhat[5] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
		thetaHat[t][4] <- Nhat[6] / (Nhat[3] + Nhat[4] + Nhat[5] + Nhat[6]);
	}
	
}
