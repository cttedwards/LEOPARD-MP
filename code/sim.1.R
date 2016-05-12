library(leopard)
## leopard version 0.0.0.9106 (2016-04-29 10:24:16)
library(reshape2)
library(ggplot2)

# survival rates per demographic category
ss <- c(0.3270833, 0.7197452, 0.9615385, 0.8882435, 0.9729730, 0.9382353, 0.9230769, 0.7219583, 0.9124424, 0.9642857, 1.0000000, 1.0000000, 0.9000000, 0.2857143)

# litter size per maternal age class
ll <- c(2, 2, 2, 2, 2)

# initial population numbers
x.initial <- c(nc  = 406,
               nj  = 201,
               saf = 92,
               f36 = 67,
               f48 = 63,
               f60 = 62,
               f72 = 93,
               f84 = 469,
               sam = 39,
               m36 = 62,
               m48 = 47,
               m60 = 78,
               m72 = 84,
               m84 = 246)

# number of iterations for simulation
niter <- 100


harvest.rate <- seq(0, 0.4, length = 101)
selectivity  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#selectivity  <- c(0, 0, rep(1,12))

prob.inf  <- array(0, dim = c(length(harvest.rate), niter, 10))
cub.surv  <- array(0, dim = c(length(harvest.rate), niter, 10))
male.surv <- array(0, dim = c(length(harvest.rate), niter, 10))
real.harem <- array(0, dim = c(length(harvest.rate), niter, 10))
real.cub.surv <- array(0, dim = c(length(harvest.rate), niter, 10))
prob.ext <- array(0, dim = c(length(harvest.rate), niter, 10))

# loop over harvest rates
for (i in 1:length(harvest.rate)) {
  
  # create new leopard object
  ss.harvest <- ss * (1 - selectivity * harvest.rate[i])
  xx <- leopard(x.initial, survival.rates = ss.harvest, litter.sizes = ll, harem.size.min = 1.5)
  
  # loop over iterations
  for (j in 1:niter) {
    
    # re-initialise
    yy <- xx
    
    # project forward to adjust numbers
    for (k in 1:100) {
      
      yy <- birth(survival(yy))
      
      # for every 10th year
      if (k %% 10 == 0) {
        # record probability of infanticide
        # and cub survival
        prob.inf[i, j, k/10]   <- yy@prob.infanticide
        cub.surv[i, j, k/10]   <- yy@expected.survival.rate[1] * (1 - yy@prob.infanticide)
        male.surv[i, j, k/10]  <- sum(yy@expected.survival.rate[9:14])/6
        real.harem[i, j, k/10] <- yy@realised.harem.size
        real.cub.surv[i, j, k/10]   <- yy@realised.survival.rate[1] * (1 - yy@prob.infanticide)
        prob.ext[i, j, k/10]  <- 1 - if(sum(yy@.Data) / sum(x.initial)>1){
                                     1
                                     } else{
                                       sum(yy@.Data) / sum(x.initial)
                                     }
        
        
      }
      
      yy <- transition(yy)
      
    }
  }
}

# plot median probability of infanticide
prob.inf <- apply(prob.inf, c(1, 3), median)
dimnames(prob.inf) <- list(H = harvest.rate, Year = seq(10, 100, length = 10))
ggplot(melt(prob.inf)) + geom_line(aes(H, value, col = as.factor(Year))) + 
  ggtitle('Probability of infanticide') + labs(y = '', col = 'Year')

# plot median cub survivorship
cub.surv <- apply(cub.surv, c(1, 3), median)
dimnames(cub.surv) <- list(H = harvest.rate, Year = seq(10, 100, length = 10))
ggplot(melt(cub.surv)) + geom_line(aes(H, value, col = as.factor(Year))) + 
  ggtitle('Expected cub survival') + labs(y = '', col = 'Year')

# plot median male survivorship
male.surv <- apply(male.surv, c(1, 3), median)
dimnames(male.surv) <- list(H = harvest.rate, Year = seq(10, 100, length = 10))
ggplot(melt(male.surv)) + geom_line(aes(H, value, col = as.factor(Year))) + 
  ggtitle('Expected male survival') + labs(y = '', col = 'Year')

# plot median harem size
real.harem <- apply(real.harem, c(1, 3), median)
dimnames(real.harem) <- list(H = harvest.rate, Year = seq(10, 100, length = 10))
ggplot(melt(real.harem)) + geom_line(aes(H, value, col = as.factor(Year))) + 
  ggtitle('Realised harem size') + labs(y = '', col = 'Year')

# plot median cub survivorship
real.cub.surv <- apply(real.cub.surv, c(1, 3), median)
dimnames(real.cub.surv) <- list(H = harvest.rate, Year = seq(10, 100, length = 10))
ggplot(melt(real.cub.surv)) + geom_line(aes(H, value, col = as.factor(Year))) + 
  ggtitle('Realised cub survival') + labs(y = '', col = 'Year')

# plot median pop size
prob.ext <- apply(prob.ext, c(1, 3), median)
dimnames(prob.ext) <- list(H = harvest.rate, Year = seq(10, 100, length = 10))
ggplot(melt(prob.ext)) + geom_line(aes(H, value, col = as.factor(Year))) + 
  ggtitle('Extinction probability') + labs(y = '', col = 'Year')








