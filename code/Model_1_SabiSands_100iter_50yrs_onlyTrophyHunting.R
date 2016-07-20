###########################################################################################
# Model Run 1 | Sabi Sands - 100 iterations - 50 years - only trophy hunting
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

###########################################################################################
# Setup model objects
###########################################################################################

# dimensions
## number of monte carlo samples
nreps <- 100
## number of projection years
nyr.proj <- 50

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
#               nj  = 5,
#               saf = 2,
#               f36 = 1,
#               f48 = 1,
#               f60 = 1,
#               f72 = 1,
#               f84 = 3,
#               sam = 3,
#               m36 = 1,
#               m48 = 1,
#               m60 = 1,
#               m72 = 1,
#               m84 = 2)

# population projection array
x <- array(x.initial,dim=c(length(x.initial),nreps,nyr.proj))
dimnames(x) <- list(age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)
x[,,2:nyr.proj] <- NA

# harvest rate
harvest.rate <- seq(0, 0.6, length = 101)

# setup selectivity object
selectivity <- matrix(data = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,  # all males
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,  # males ≥ 3
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  # males ≥ 6
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  # males ≥ 7
                               0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,  # males & females ≥ 6
                               0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,  # males & females ≥ 7
                               0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  # excl dependent young
                               0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,  # only adults
                               0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1), # all males with females ≥ 7
                      nrow = 14) 

selectivity <- t(selectivity) ; rownames(selectivity) = c("all.males", "males≥3", "males≥6", "males≥7", "male.female≥6",
                                                          "male.female≥7", "excl.dep.young", "only.adults", "all.males.w.females≥7")

# create population size object
population.size <- array(0, dim = c(length(harvest.rate), nreps, 5))

# create extinction prob object
extinction.probability <-  as.data.frame(NA)

# query objects
#eigen.df <- data.frame()
total.removed <- data.frame()
#pop.size <- data.frame()

prob.inf     <- array(0, dim = c(length(harvest.rate), nreps, 5))
prob.inf.df  <- as.data.frame(NA)

cub.surv    <- array(0, dim = c(length(harvest.rate), nreps, 5))
cub.surv.df <- as.data.frame(NA)

# query number killed
total.harvested    <- array(0, dim = c(length(harvest.rate), nreps, 5))
total.harvested.df <- as.data.frame(NA)

###########################################################################################
# Functions
###########################################################################################

prob.ext.func <- function(x){
  prob.extinction <- 1 - mean(x > 0)
  return(prob.extinction)
}

###########################################################################################
# Run model
###########################################################################################

# loop over range of selectivities
for(z in 1:nrow(selectivity)){
  
  # loop over harvest rates
  for (k in 1:length(harvest.rate)) { 
    
    # projection
    for (i in 1:nreps) {
      
      # sample parameters
      param.sample <- params[,i]
      
      # create new object to hold leopard numbers
      # and vital rates
      xx <- leopard(x.initial, param.sample[1:14], param.sample[15:19])
      
      # assign multiplicative maternal effects
      xx@maternal.effect[] <- matrix(maternal.effects, nrow=2, ncol=5, byrow=T)
      
      # loop forward over years
      for (y in 1:nyr.proj) {
        
        # correlated deviation in survival: 
        # log-normal with cv = 0.2
        # truncated at 0 and 1
        sdev <- rnorm(1)
        sigma <- sqrt(log(1+0.20^2))
        
        param.sample[1:14] <- exp(log(param.sample[1:14]) +  sigma * sdev - sigma^2/2)
        param.sample[1:14] <- vapply(vapply(param.sample[1:14],function(x) max(x,0),numeric(1)),function(x) min(x,1),numeric(1))
        
        # create list of sequential hunting scenarios
        removals <- list(trophy = list(rate = harvest.rate[k], preference = selectivity[z,]), 
                         problem_animal = list(rate = 0.0))
        
        removals <- harvest(xx, removals)
        
        # include trophy hunting aging error 
        #source('incorp.aging.error.final.r')
        
        total.removals <- removals$trophy@kills + removals$problem_animal@kills
        
        # add recovery years
        #source('two.years.recovery.r')
        #source('three.years.recovery.r')
        
        # calculate stochastic survival
        xx <- survival(xx, total.removals)
        
        # calculate stochastic birth
        xx <- birth(xx)
        
        # for every 10th year
        if (y %% 10 == 0) {
          # record probability of infanticide
          # and cub survival
          prob.inf[k, i, y/10]   <- xx@prob.infanticide
          cub.surv[k, i, y/10]   <- xx@expected.survival.rate[1] * (1 - xx@prob.infanticide)
          population.size[k, i, y/10]  <- sum(xx@.Data)
          total.harvested[k, i, y/10]  <- sum(total.removals)
          
        }
        
        #eigen <- Re(eigen(tmatrix(xx))$values)[1]
        #print(Re(eigen(tmatrix(xx))$values)[1])
        
        # step forward
        xx <- transition(xx)
        
        # record numbers
        x[,i,y] <- xx
        
        #eigen.df <- rbind(eigen.df, eigen)
        #total.removed   <- rbind(total.removed, total.removals)
        #pop.size <- rbind(pop.size, xx@.Data)
        
      }
      
    }
    
  }
  
  extinction.probability <- cbind(extinction.probability, apply(population.size, c(1, 3), prob.ext.func))
  total.harvested.df <- cbind(total.harvested.df, apply(total.harvested, c(1, 3), mean))
  cub.surv.df <- cbind(cub.surv.df, apply(cub.surv, c(1, 3), median))
  prob.inf.df <- cbind(prob.inf.df, apply(prob.inf, c(1, 3), median))
  
}


# save output
saver(extinction.probability,
      total.harvested.df,
      cub.surv.df,
      prob.inf.df,
      name = 'model_run_1')
