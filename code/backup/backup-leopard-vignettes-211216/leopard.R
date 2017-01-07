## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(leopard)

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

# create new leopard object
xx <- leopard(x.initial, survival.rates = ss, litter.sizes = ll)

# examine leopard-class object
xx

## ------------------------------------------------------------------------
yy <- survival(xx)
yy@realised.survival.rate
yy <- infanticide(xx)
yy@realised.survival.rate
yy@.Data * yy@realised.survival.rate

## ------------------------------------------------------------------------
yy <- transition(survival(xx))
yy@.Data

## ------------------------------------------------------------------------
implementation(xx, list(trophy = list(size = 50)))

## ------------------------------------------------------------------------
implementation(xx, list(trophy = list(size = 50, preference = "Braczkowski"), problem_animal = list(size = 100)))

## ------------------------------------------------------------------------
hcr <- new('control')
hcr@.Data <- function(x) floor(sum(x * hcr@target_preference) * hcr@target_harvest_rate)
hcr@target_harvest_rate <- 0.1

hcr(xx)

