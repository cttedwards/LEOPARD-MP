###########################################################################################
# Load functions and libraries
###########################################################################################

rm(list = ls())
setwd("/Users/RossTyzackPitman/Documents/OneDrive/Data/GitHub/Databases/LEOPARD-MP/code")

library(leopard)
library(ggplot2)
library(reshape2)

source('utils/saver.r')
source('utils/loader.r')
source('utils/reader.r')
source('utils/pdfr.r')
source('params.R')
#source('params.kzn.R')
source('aging_error.r')

# extinction probability function
prob.ext.func <- function(x){
  prob.extinction <- 1 - mean(x > 0)
  return(prob.extinction)
}

###########################################################################################
# Model.1 Description
###########################################################################################

# Site parameters         - Sabi Sands
# Harvest scenario        - ≥7 years old
# Porportion removed      - 1 (100%)
# Problem animal control  - NA
# Non-compliance          - NA
# Aging error             - NA
# Recovery years          - NA

###########################################################################################
# Setup model objects
###########################################################################################

# unadjusted cub survival
#param[1] <- 0.4520594

# dimensions
## number of monte carlo samples
nreps <- 1000
## number of projection years
nyr.proj <- 60

# query objects
eigen.value <- matrix(ncol = nreps, nrow = nyr.proj)
harvest.value <- matrix(ncol = nreps, nrow = nyr.proj)

## matrix of paramter values
params <- matrix(param,nrow=length(param),ncol=nreps)  
rownames(params) <- names(param)                       
colnames(params) <- 1:nreps

# initial population size (Sabi Sands average population structure estimate from 2013-2015; using 70 leopard)
x.initial <- c(nc  = 14,
               nj  = 7,
               saf = 3,
               f36 = 2,
               f48 = 2,
               f60 = 2,
               f72 = 3,
               f84 = 17,
               sam = 1,
               m36 = 2,
               m48 = 2,
               m60 = 3,
               m72 = 3,
               m84 = 9)

# initial population size (Phinda average population structure estimate from 2002-2012; using 30 leopard)
#x.initial <- c(nc  = 7,
#nj  = 5,
#saf = 2,
#f36 = 1,
#f48 = 1,
#f60 = 1,
#f72 = 1,
#f84 = 3,
#sam = 3,
#m36 = 1,
#m48 = 1,
#m60 = 1,
#m72 = 1,
#m84 = 2)

x.initial <- x.initial * 64 

# population projection array
x <- array(x.initial,dim=c(length(x.initial),nreps,nyr.proj))
dimnames(x) <- list(age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)
x[,,2:nyr.proj] <- NA

# harvest rate
harvest.rate <- 1

# setup selectivity object
selectivity <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
only.males  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)

###########################################################################################
# Run model
###########################################################################################

# projection
for (i in 1:nreps) {
  
  # sample parameters
  param.sample <- params[,i]
  
  # create new object to hold leopard numbers
  # and vital rates
  xx <- leopard(numbers        = x.initial, 
                survival.rates = param.sample[1:14], 
                litter.sizes   = param.sample[15:19], 
                deterministic  = FALSE,
                harem.size     = 1.14)
  
  # assign multiplicative maternal effects
  xx@maternal.effect[] <- matrix(maternal.effects, nrow = 2, ncol = 5, byrow = T)
  
  # record initial numbers
  x[,i,1] <- xx
  
  # loop forward over years
  for (y in 2:nyr.proj) {
    
    # total numbers
    xx.t0 <- sum(xx)
    
    # correlated deviation in survival: 
    # log-normal with cv = 0.2
    # truncated at 0 and 1
    sdev <- rnorm(1)
    sigma <- sqrt(log(1+0.20^2))
    
    param.sample[1:14] <- exp(log(param.sample[1:14]) +  sigma * sdev - sigma^2/2)
    param.sample[1:14] <- vapply(vapply(param.sample[1:14],
                                        function(x) max(x,0),numeric(1)),function(x) min(x,1),numeric(1))
    
    # create list of sequential hunting scenarios (i.e., Braczkowski or selectivity)
    removals <- list(trophy         = list(rate = harvest.rate, preference = selectivity), 
                     problem_animal = list(rate = 0.00),
                     noncompliance  = list(rate = 0.00),
                     aging_error    = list(rate = 0.00, preference = only.males))
    
    removals <- harvest(xx, removals)
    
    # include trophy hunting aging error 
    #source('incorp.aging.error.final.r')
    
    total.removals <- removals$trophy@kills + removals$problem_animal@kills + removals$noncompliance@kills + removals$aging_error@kills
    
    # run initialisation
    source('initialisation.r')
    
    # add recovery years
    #source('two.years.recovery.r')
    #source('three.years.recovery.r')
    
    # calculate stochastic survival
    xx <- survival(xx, total.removals)
    
    # calculate infanticide
    xx <- infanticide(xx)
    
    eigen <- Re(eigen(tmatrix(xx))$values)[1]
    
    # step forward
    xx <- transition(xx)
    
    # calculate stochastic birth
    #xx <- birth(xx)
    
    # record numbers
    x[,i,y] <- xx
    
    # calculate eigen value
    xx.t1 <- sum(xx)
    
    eigen.value[y, i] <- xx.t1 / xx.t0
    harvest.value[y, i] <- sum(total.removals)
    
  }
  
}

dimnames(eigen.value) <- list(year = 1:nyr.proj, iter = 1:nreps)
dimnames(harvest.value) <- list(year = 1:nyr.proj, iter = 1:nreps)

saver(x,
      eigen.value,
      name = 'model.1')

###########################################################################################
# Plot
###########################################################################################

loader('model.1')

# average growth rate
mean(eigen.value, na.rm = TRUE)

# total popualtion size
x.tot <- apply(x, 2:3, sum)
par(bg = NA) 
#par(bg = "white") 
cairo_pdf(file = '/Users/RossTyzackPitman/Documents/OneDrive/Data/GitHub/Databases/PhD_Chapter3/MSE_Paper/figures/Population.Size.Model.1.pdf', 
    width = 8, height = 5)
boxplot(x.tot[,11:nyr.proj],
        ylab = "Population Size",
        xaxt = "n",
        xlab = "Year",
        ylim = c(0,10000),
        outline = FALSE)
title("Hunting = 100% males \u2265 7 yrs", line = -2)
axis(side = 1, at = 1:nyr.proj)
dev.off()

par(bg = NA) 
#par(bg = "white") 
pdf(file = '/Users/RossTyzackPitman/Documents/OneDrive/Data/GitHub/Databases/PhD_Chapter3/MSE_Paper/figures/Eigen.Value.Model.1.pdf', 
    width = 8, height = 5)
boxplot(t(eigen.value),
        ylab = expression(paste("Growth rate (", lambda, ")")),
        xaxt = "n",
        xlab = "Year",
        ylim = c(0.7,1.2),
        outline = FALSE)
axis(side = 1, at = 1:nyr.proj)
abline(h = 1, col = 2, lty = 1)
dev.off()

par(bg = NA) 
#par(bg = "white") 
pdf(file = '/Users/RossTyzackPitman/Documents/OneDrive/Data/GitHub/Databases/PhD_Chapter3/MSE_Paper/figures/Harvest.Value.Model.1.pdf', 
    width = 8, height = 5)
boxplot(t(harvest.value),
        ylab = "Number harvested",
        xaxt = "n",
        xlab = "Year",
        #ylim = c(0.7,1.2),
        outline = FALSE)
axis(side = 1, at = 1:nyr.proj)
#abline(h = 1, col = 2, lty = 1)
dev.off()
