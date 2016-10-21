model{
  
  for (i in 1:n.stocks){
    
  # censoring mechanism  
  censored[i] ~ dinterval(a.time[i],ctime[i])
  
  # sampling model
  a.time[i] ~ dlnorm(mu[i],tau)
  
  # linear predictor
  log(mu[i]) <- betas[1:n.covs] %*% COVS[i,1:n.covs] + habitat[hab[i]] + taxonomy[tax[i]]
  
  CS[i] <- -log(1-plnorm(ctime[i],mu[i],tau))
  
  
  }
  
  # regression parameters
  for (k in 1:n.covs){
    betas[k] ~ dnorm(0,1e-6)
  }
  
  ###
  #random effects 
  ###
  
  for (j in 1:n.hab){
    habitat[j] <- hab.xi*hab.eta[j]
    hab.eta[j] ~ dnorm(0,hab.prec)
  }
  # half cauchy prior on random effect variance
  hab.xi ~ dnorm(0,0.04)
  hab.prec ~ dgamma(0.5,0.5)
  
  for (l in 1:n.tax){
    taxonomy[l] <- tax.xi*tax.eta[l]
    tax.eta[l] ~ dnorm(0,tax.prec)
  }
  # finite population variance
  fp.sd.tax <- sd(taxonomy)
  fp.sd.hab <- sd(habitat)
  
  # half cauchy on random effects
  tax.xi ~ dnorm(0,0.04)
  tax.prec ~ dgamma(0.5,0.5)
  
  tau ~ dgamma(0.1,0.1)
  
}