
###########################################################################################
  # Load packages and utility functions
###########################################################################################

rm(list = ls())
setwd("/Users/RossTyzackPitman/Documents/OneDrive/PhD/Data/R_Database/R_PROJECTS/PhD_Chapter3/MSE_Paper/LEOPARD-MP/code")

library(leopard)

source('utils/saver.r')
source('utils/loader.r')
source('utils/reader.r')
source('params.R')

# dimensions

## number of monte carlo samples
nreps <- 100
## number of projection years
nyr.proj <- 20

## matrix of paramter values
params <- matrix(param,nrow=length(param),ncol=nreps)  
rownames(params) <- names(param)                       
colnames(params) <- 1:nreps

# initial population size (Sabi Sands average population structure estimate from 2009-2015; using 2010 leopard)
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

# population projection array
x <- array(x.initial,dim=c(length(x.initial),nreps,nyr.proj))
dimnames(x) <- list(age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)
x[,,2:nyr.proj] <- NA

# projection
for (i in 1:nreps) {
  
  # sample parameters
  param.sample <- params[,i]
  
  # create new object to hold leopard numbers
  # and vital rates
  xx <- leopard(x.initial, param.sample[1:14], param.sample[15:19], harem.size = 1.25)
  
  # assign multiplicative maternal effects
  xx@maternal.effect[] <- matrix(maternal.effects, nrow=2, ncol=5, byrow=T)
  
  # loop forward over years
  for (y in 2:nyr.proj) {
    
    # correlated deviation in survival: 
    # log-normal with cv = 0.2
    # truncated at 0 and 1
    sdev <- rnorm(1)
    sigma <- sqrt(log(1+0.20^2))
    
    # (check this: if you take 1000 values of sdev the right hand side should have a mean value
    # approximately equal to param.sample)
    param.sample[1:14] <- exp(log(param.sample[1:14]) +  sigma * sdev - sigma^2/2)
    param.sample[1:14] <- vapply(vapply(param.sample[1:14],function(x) max(x,0),numeric(1)),function(x) min(x,1),numeric(1))
    
    # calculate stochastic survival
    removals <- implementation(xx, 0, preference = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
    xx <- survival(xx, removals@kills)
    
    #xx <- survival(xx)
    
    # calculate stochastic birth
    xx <- birth(xx)
    #eigen <- Re(eigen(tmatrix(xx))$values)[1]
    #print(Re(eigen(tmatrix(xx))$values)[1])
    
    # step forward
    xx <- transition(xx)
    
    # record numbers
    x[,i,y] <- xx
    
  }
  
}

###########################################################################################
# PLOTS
###########################################################################################

x.tot <- apply(x, 2:3, sum)
boxplot(x.tot, ylab = "N", xaxt = "n", xlab = "Step", outline = FALSE, ylim = c(0, 4000))
axis(side = 1, at = 1:nyr.proj)

###########################################################################################
# END
###########################################################################################


#Re(eigen(tmatrix(xx))$values)[1]

#xx@.Data
#sum(xx@.Data[1:2]) / sum(xx@.Data[3:14])

#x.tot

#tmatrix(xx)
#eigen(tmatrix(xx))
#mean(eigen.df[,1])

###########################################################################################
# Adding problem leopard offtake
###########################################################################################

#xx - rep(1,14)
#xx@maternal.birth[1,] - 1:5





