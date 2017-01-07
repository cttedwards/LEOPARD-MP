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

