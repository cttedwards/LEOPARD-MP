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

eigen.df <- data.frame()

# dimensions

## number of monte carlo samples
nreps <- 100
## number of projection years
nyr.proj <- 30

## matrix of paramter values
params <- matrix(param,nrow=length(param),ncol=nreps)  
rownames(params) <- names(param)                       
colnames(params) <- 1:nreps

# initial population size (Swanepoel et al. 2014; ±2000 leopard)
x.initial <- c(nc=694,
               nj=410,
               saf=117,
               f36=99,
               f48=80,
               f60=73,
               f72=45,
               f84=194,
               sam=84,
               m36=41,
               m48=24,
               m60=24,
               m72=28,
               m84=95)

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
  xx <- leopard(x.initial, param.sample[1:14], param.sample[15:19], harem.size = 5)
  
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
    removals <- implementation(xx, 50, preference = c(rep(0, 9), rep(1, 5)))

    xx <- survival(xx, removals@kills)
    #xx <- survival(xx)

    # calculate stochastic birth
    xx <- birth(xx)
    eigen <- Re(eigen(tmatrix(xx))$values)[1]
    #print(Re(eigen(tmatrix(xx))$values)[1])
    
    # step forward
    xx <- transition(xx)
    
    # record numbers
    x[,i,y] <- xx
   
    eigen.df <- rbind(eigen.df, eigen)
    
     
  }
  
}

###########################################################################################
# PLOTS
###########################################################################################

x.tot <- apply(x, 2:3, sum)
boxplot(x.tot, ylab = "N", xaxt = "n", xlab = "Step", outline = FALSE)
axis(side = 1, at = 1:nyr.proj)

###########################################################################################
# END
###########################################################################################


#Re(eigen(tmatrix(xx))$values)[1]

xx@.Data
sum(xx@.Data[1:2]) / sum(xx@.Data[3:14])
#x.tot

#tmatrix(xx)
#eigen(tmatrix(xx))


mean(eigen.df[,1])
