## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(leopard)
library(reshape2)
library(ggplot2)

## ------------------------------------------------------------------------
# survival rates per demographic category
#ss <- c(0.3270833, 0.7197452, 0.9615385, 0.8882435, 0.9729730, 0.9382353, 0.9230769, 0.7219583, 0.9124424, 0.9642857, 1.0000000, 1.0000000, 0.9000000, 0.2857143)
ss <- c(0.4610291, 0.7197452, 0.9615385, 0.8882435, 0.9729730, 0.9382353, 0.9230769, 0.7219583, 0.9124424, 0.9642857, 1.0000000, 1.0000000, 0.9000000, 0.2857143)

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
niter <- 1e2

## ------------------------------------------------------------------------
xx <- leopard(x.initial, survival.rates = ss, litter.sizes = ll)
harem(xx)@harem.size

## ------------------------------------------------------------------------
xx <- leopard(x.initial, survival.rates = ss, litter.sizes = ll)

## ---- fig.width=6--------------------------------------------------------
simulated.survival <- numeric(1e4)
for (i in 1:length(simulated.survival)) 
    simulated.survival[i] <- infanticide(survival(xx))@realised.survival.rate[1]

hist(simulated.survival, prob = TRUE)
abline(v = 0.3270833, col = 2)

## ---- fig.width=6, warning=FALSE-----------------------------------------
harvest.rate <- seq(0, 0.4, length = 11)
selectivity  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)

prob.inf <- array(0, dim = c(length(harvest.rate), niter, 5))
cub.surv <- array(0, dim = c(length(harvest.rate), niter, 5))

# create new leopard object
xx <- leopard(x.initial, survival.rates = ss, litter.sizes = ll)
    
# loop over harvest rates
for (i in 1:length(harvest.rate)) {
    
    # loop over iterations
    for (j in 1:niter) {
        
        # re-initialise
        yy <- xx
        
        yy <- survival(yy)
        yy <- infanticide(yy)
        yy <- transition(yy)
        
        # project forward to adjust numbers
        for (k in 1:50) {
            
            yy <- birth(yy)
            
            zz <- harvest(yy, list(list(rate = harvest.rate[i], preference = selectivity)))
            
            yy <- survival(yy, zz[[1]]@kills)
            yy <- infanticide(yy)
            
            # for every 10th year
            if (k %% 10 == 0) {
                # record probability of infanticide
                # and cub survival
                prob.inf[i, j, k/10] <- yy@prob.infanticide
                cub.surv[i, j, k/10] <- yy@expected.survival.rate[1] * (1 - yy@prob.infanticide)
            }
                
            yy <- transition(yy)
            
        }
    }
}

# plot median probability of infanticide
prob.inf <- apply(prob.inf, c(1, 3), median)
dimnames(prob.inf) <- list(H = harvest.rate, Year = seq(10, 50, length = 5))
ggplot(melt(prob.inf)) + geom_line(aes(H, value, col = as.factor(Year))) + 
    ggtitle('Probability of infanticide') + labs(y = '', col = 'Year')

# plot median cub survivorship
cub.surv <- apply(cub.surv, c(1, 3), median)
dimnames(cub.surv) <- list(H = harvest.rate, Year = seq(10, 50, length = 5))
ggplot(melt(cub.surv)) + geom_line(aes(H, value, col = as.factor(Year))) + 
    ggtitle('Expected cub survival') + labs(y = '', col = 'Year')

## ------------------------------------------------------------------------
xx <- leopard(x.initial, survival.rates = ss, litter.sizes = ll)

## ---- fig.width=6--------------------------------------------------------
simulated.days <- numeric(niter)
for (i in 1:length(simulated.days)) 
    simulated.days[i] <- implementation(xx, list(list(size = 1)))[[1]]@days[]

hist(simulated.days, prob = TRUE)

## ---- fig.width=6--------------------------------------------------------
simulated.age.class <- numeric(niter)
for (i in 1:length(simulated.age.class)) 
    simulated.age.class[i] <- which(implementation(xx, list(list(size = 1, preference = "Braczkowski")))[[1]]@kills[] > 0)

hist(simulated.age.class, prob = TRUE)

