
data {

	int T;
	int T1;
	int T2;
	real y[T];
	real f[T-1];
	int m[T1,T2+1];

}
parameters {

	real alpha1;
	real alphaa;
	real alphar;
	real alphal;
	real beta1;
	real betaa;
	real betar;
	real betal;
	real<lower=0.01,upper=10> sigy;
	
	real<lower=0,upper=1e5> N1[T];
	real<lower=0,upper=1e5> Na[T];

}
transformed parameters {

	real phi1[T-1];
	real phia[T-1];
	real rho[T-1];
	real lambda[T-1];
	
	vector<lower=0,upper=1>[T2+1] p[T1];
	
	#print("T1:",T1);
	#print("m[2]:", m[2]);
	#print("m[32]:", m[32]);
	#print("m[T1]:", m[T1]);
	
	# logistic regression equations
	for (t in 1:(T-1)) {
		
		# survivorship as a function of normalised frost
		phi1[t] = inv_logit(alpha1 + beta1 * f[t]);
		phia[t] = inv_logit(alphaa + betaa * f[t]);
		
		# productivity as a function of time
		rho[t] = exp(alphar + betar * t); 
		 
		# reporting rate as a function of time
		lambda[t] = inv_logit(alphal + betal * (t+1));
	}
	
	#print("phi1:", phi1);
	#print("phia:", phia);
	#print("rho:", rho);
	#print("lambda:", lambda);
	
	# recovery table probabilities
	for (t1 in 1:(T1-1)) {
	
		# calculate the offset diagonal
		# (release and recover the next year)
		p[t1, t1] = lambda[t1] * (1-phi1[t1]);
	
		# calculate value one above the diagonal
		# (survive as chicks but die as adults)
		p[t1, t1+1] = lambda[t1+1] * phi1[t1] * (1-phia[t1+1]);
	
		# calculate remaining terms above diagonal
		for(t2 in (t1+2):T2) {
		
			vector[T2] lphi;
			for(t in (t1+1):(t2-1)){
				lphi[t] = log(phia[t]);
			}
			
			p[t1,t2] = lambda[t2] * phi1[t1] * (1 - phia[t2]) * exp(sum(lphi[(t1+1):(t2-1)]));
		}
		
		for(t2 in 1:(t1-1)){
			# zero probabilities in lower triangle of table
			p[t1, t2] = 0;
		}
	
		# probability of an animal never being seen again
		p[t1, T2+1] = 1 - sum(p[t1, 1:T2]);
	}
	
	# Final row
	p[T1,T1] = lambda[T1] * (1 - phi1[T1]);
	for(t2 in 1:(T1-1)){
		p[T1,t2] = 0;
	}
	p[T1,T2+1] = 1 - p[T1,T2];
	
	#print("p[2]:", p[2]);
	#print("p[32]:", p[32]);
	#print("p[T1]:", p[T1]);
	
}
model {

	# logistic regression parameters
	alpha1 ~ normal(0,1);
	alphaa ~ normal(0,1);
	alphar ~ normal(0,1);
	alphal ~ normal(0,1);
	beta1 ~ normal(0,1);
	betaa ~ normal(0,1);
	betar ~ normal(0,1);
	betal ~ normal(0,1);
	
	# observation error prior
	sigy ~ inv_gamma(0.001,0.001);
	
	# initial population priors
	for(t in 1:2){
		N1[t] ~ normal(200, 0.000001);
		Na[t] ~ normal(1000, 0.000001);
	}
	
	# model dynamics using the Normal approximation
	# to Poisson/Binomial process
	for(t in 3:T){
	
		real mean1;
		real meana;
		real sig1;
		real siga;
		
		mean1 = rho[t-1] * phi1[t-1] * Na[t-1];
		meana = phia[t-1] * (N1[t-1] + Na[t-1]);
		
		sig1 = Na[t-1] * rho[t-1] * phi1[t-1];
		siga = (N1[t-1]+Na[t-1]) * phia[t-1] * (1-phia[t-1]);
		
		N1[t] ~ normal(mean1,sig1);
		Na[t] ~ normal(meana,siga);
	}
	
	# Define the observation process for the census/index data
	for(t in 3:T){
		y[t] ~ normal(Na[t],sigy);
	}
	
	# Define the recovery likelihood
	for(t in 1:T1){
		m[t] ~ multinomial(p[t]);
	}
}
generated quantities {

	vector[T] r;

	for (t in 1:3) {
		r[t] = -1; 
	}
	for (t in 4:T){
		r[t] = (Na[t]+N1[t]) / (Na[t-1]+N1[t-1]);
	}
}
