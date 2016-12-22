###########################################################################################
# Setup model objects
###########################################################################################

rm(list = ls())
setwd("..")

library(leopard)
library(ggplot2)
library(reshape2)

source('utils/saver.r')
source('utils/loader.r')
source('params.R')
source('aging_error.r')

# extinction probability function
prob.ext.func <- function(x){
  prob.extinction <- 1 - mean(x > 0)
  return(prob.extinction)
}

# unadjusted cub survival
param[1] <- 0.4610291

# dimensions
## number of monte carlo samples
nreps <- 1000
## number of projection years
nyr.proj <- 60

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

x.initial.multiplier <- c(1, 10, 20, 60)

# population projection array
x <- array(NA,dim=c(length(x.initial.multiplier), length(x.initial),nreps,nyr.proj))
dimnames(x) <- list(multiplier = x.initial.multiplier, age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)

# harvest rate
harvest.rate <- 0

# setup selectivity object
selectivity <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)

###########################################################################################
# Run model
###########################################################################################

for (m in 1:length(x.initial.multiplier)) {
  
  # projection
  for (i in 1:nreps) {
    
    # sample parameters
    param.sample <- params[,i]
    
    # create new object to hold leopard numbers
    # and vital rates
    xx <- leopard(x.initial * x.initial.multiplier[m], param.sample[1:14], param.sample[15:19], deterministic = FALSE)
    
    # assign multiplicative maternal effects
    xx@maternal.effect[] <- matrix(maternal.effects, nrow = 2, ncol = 5, byrow = T)
    
    # record numbers
    x[m,,i,1] <- xx
    
    # loop forward over years
    for (y in 2:nyr.proj) {
      
      # correlated deviation in survival: 
      # log-normal with cv = 0.2
      # truncated at 0 and 1
      sdev <- rnorm(1)
      sigma <- sqrt(log(1+0.20^2))
      
      param.sample[1:14] <- exp(log(param.sample[1:14]) +  sigma * sdev - sigma^2/2)
      param.sample[1:14] <- vapply(vapply(param.sample[1:14],
                                          function(x) max(x,0),numeric(1)),function(x) min(x,1),numeric(1))
      
      # create list of sequential hunting scenarios
      removals <- list(trophy = list(rate = harvest.rate, preference = selectivity), 
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
      
      # calculate infanticide
      xx <- infanticide(xx)
      
      # step forward
      xx <- transition(xx)
      
      # calculate stochastic birth
      xx <- birth(xx)
      
      # record numbers
      x[m,,i,y] <- xx
      
    }
    
  }
}


saver(x,
      name = 'SSGR_SensitivityAnalysis')

###########################################################################################
# Plots
###########################################################################################

loader("SSGR_SensitivityAnalysis")

# sum across age classes
x <- apply(x, c(1, 3, 4), sum)

# total population size
x.tot <- apply(x, c(1,3), mean)

x.tot <- melt(x.tot)

#levels(as.factor(x.tot$multiplier))

x.tot$multiplier <- factor(as.factor(x.tot$multiplier), levels = rev(levels(as.factor(x.tot$multiplier))))

#levels(x.tot$multiplier)

par(bg = NA) 
#par(bg = "white") 
pdf(file = '/Users/RossTyzackPitman/Documents/OneDrive/Data/GitHub/Databases/PhD_Chapter3/MSE_Paper/figures/Population.Size.Sensitivity.pdf', 
    width = 8, height = 5)

xlabels <- c("0", "10", "20", "30", "40", "50")

pop.size.sensitivity <- 
  ggplot(subset(x.tot, year > 10)) + 
  geom_line(aes(year, value, col = as.factor(multiplier)), linetype = "solid") +
  xlab("Year\n") + ylab("Population size\n") +
  theme_bw() + theme(strip.background = element_rect(fill = "white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  
  scale_x_continuous(labels = xlabels) +

  
  #theme(legend.position = c(.5, .5)) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_rect(colour = NA)) +
  scale_colour_grey() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.8))

pop.size.sensitivity

dev.off()





# probability of extinction

#x.tot <- apply(x, c(1,3), prob.ext.func)

#ggplot(melt(x.tot)) + geom_line(aes(year, value, col = as.factor(multiplier)))


# change in population size

#x.tot <- sweep(x, 1, sum(x.initial) * x.initial.multiplier, '/')
#x.tot <- apply(x.tot, c(1,3), mean)

#ggplot(melt(x.tot)) + geom_line(aes(year, value, col = as.factor(multiplier)))

