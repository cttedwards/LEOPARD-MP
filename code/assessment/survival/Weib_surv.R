model{
  
  for (i in 1:n.stocks){
    
  # censoring mechanism  
  censored[i] ~ dinterval(a.time[i],ctime[i])
  
  # sampling model
  a.time[i] ~ dweib(tau,mu[i])
  
  # linear predictor
  mu[i] <- exp(logmean+betas[1:n.covs] %*% COVS[i,1:n.covs] + region[reg[i]] + habitat[hab[i]] + classfx[class[i]] + familyfx[family[i]] + orderfx[order[i]])
  
  # cox-snell residuals - use pweib alone for predictions of P(t<T|X) i.e. for some predictor X what is the probability that a stock is assessed within T years
  CS[i] <- -log(1-pweib(ctime[i],tau,mu[i]))
 
  }
  
  # regression parameters
  for (k in 1:n.covs){
    betas[k] ~ dnorm(0,1e-6)
  }
  
  ###
  #random effects 
  ###
  
  for (j in 1:n.hab){
    hab.pmu[j] <- exp(logmean+ habitat[j] + classfx.pred + orderfx.pred + familyfx.pred + regfx.pred)
    for (t in tmin:tmax){
       hab.pred[j,t] <- pweib(t,tau,hab.pmu[j])
    }
    habitat[j] <- hab.xi*hab.eta[j]
    hab.eta[j] ~ dnorm(0,hab.prec)
  }
  # half cauchy prior on random effect variance
  habfx.pred <- hab.xi*hab.eta.pred
  hab.eta.pred ~ dnorm(0,hab.prec)
  hab.xi ~ dnorm(0,0.001)
  hab.prec ~ dgamma(0.5,0.5)
  
  for (l in 1:n.class){
    class.pmu[l] <- exp(logmean+ habfx.pred + classfx[l] + orderfx.pred + familyfx.pred + regfx.pred)
    for (t in tmin:tmax){
       class.pred[l,t] <- pweib(t,tau,class.pmu[l])
    }
    classfx[l] <- class.xi*class.eta[l]
    class.eta[l] ~ dnorm(0,class.prec)
  }
  # finite population variance
  fp.sd.class <- sd(classfx)
  fp.sd.habitat <- sd(habitat)
  
  # half cauchy on random effects
  classfx.pred <- class.xi*class.eta.pred
  class.eta.pred ~ dnorm(0,class.prec)
  class.xi ~ dnorm(0,0.001)
  class.prec ~ dgamma(0.5,0.5)
  
  for (l in 1:n.order){
    order.pmu[l] <- exp(logmean+ habfx.pred + classfx[classord[l]] + orderfx[l] + familyfx.pred + regfx.pred)
    for (t in tmin:tmax){
      order.pred[l,t] <- pweib(t,tau,order.pmu[l])
    }
    orderfx[l] <- order.xi*order.eta[l]
    order.eta[l] ~ dnorm(0,order.prec)
  }
  # finite population variance
  fp.sd.order <- sd(orderfx)
  
  # half cauchy on random effects
  orderfx.pred <- order.xi*order.eta.pred
  order.eta.pred ~ dnorm(0,order.prec)
  order.xi ~ dnorm(0,0.001)
  order.prec ~ dgamma(0.5,0.5)
  
  for (l in 1:n.family){
    family.pmu[l] <- exp(logmean+ habfx.pred + classfx[classfam[l]] + familyfx[l] + orderfx[orderfam[l]] + regfx.pred)
    for (t in tmin:tmax){
      family.pred[l,t] <- pweib(t,tau,family.pmu[l])
    }
    familyfx[l] <- family.xi*family.eta[l]
    family.eta[l] ~ dnorm(0,family.prec)
  }
  # finite population variance
  fp.sd.family <- sd(familyfx)
  
  # half cauchy on random effects
  familyfx.pred <- family.xi*family.eta.pred
  family.eta.pred ~ dnorm(0,family.prec)
  family.xi ~ dnorm(0,0.001)
  family.prec ~ dgamma(0.5,0.5)
  
  for (l in 1:n.region){
    region.pmu[l] <- exp(logmean+ habfx.pred + classfx.pred + region[l] + orderfx.pred + familyfx.pred)
    for (t in tmin:tmax){
      region.pred[l,t] <- pweib(t,tau,region.pmu[l])
    }
    region[l] <- region.xi*region.eta[l]
    region.eta[l] ~ dnorm(0,region.prec)
  }
  # finite population variance
  fp.sd.region <- sd(region)
  
  # half cauchy on random effects
  regfx.pred <- region.xi*region.eta.pred
  region.eta.pred ~ dnorm(0,region.prec)
  region.xi ~ dnorm(0,0.001)
  region.prec ~ dgamma(0.5,0.5)
  
  logmean ~ dnorm(0,0.00001)
  tau ~ dgamma(0.00001,0.00001)
  
}