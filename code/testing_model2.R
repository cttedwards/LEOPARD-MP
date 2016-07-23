###########################################################################################
# Setup model objects
###########################################################################################

rm(list = ls())

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

# dimensions
## number of monte carlo samples
nreps <- 100
## number of projection years
nyr.proj <- 50

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

# alternative numbers
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

# unadjusted cub survival
param[1] <- 0.4520594

# population projection array
x <- array(x.initial,dim=c(length(x.initial),nreps,nyr.proj))
dimnames(x) <- list(age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)
x[,,2:nyr.proj] <- NA

# harem size array
hh <- array(1.14, dim=c(nreps,nyr.proj))
dimnames(hh) <- list(rep = 1:nreps, year = 1:nyr.proj)
hh[,2:nyr.proj] <- NA

# survival array
ss <- array(param[1:14], dim=c(length(x.initial), nreps, nyr.proj))
dimnames(ss) <- list(age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)
ss[,,2:nyr.proj] <- NA

###########################################################################################
# Run model
###########################################################################################

# projection
for (i in 1:nreps) {
  
  # create new object to hold leopard numbers
  # and vital rates
  xx <- leopard(x.initial, param[1:14], param[15:19])
  
  # assign multiplicative maternal effects
  xx@maternal.effect[] <- matrix(maternal.effects, nrow = 2, ncol = 5, byrow = T)
  
  # loop forward over years
  for (y in 1:nyr.proj) {
    
    # record harem size
    hh[i, y] <- harem(xx)@harem.size[]
      
    # calculate stochastic birth
    xx <- birth(xx)
    
    # calculate stochastic survival
    xx <- survival(xx)
    
    # record numbers
    ss[,i,y] <- xx@realised.survival.rate
    
    # step forward
    xx <- transition(xx)
    
    # record numbers
    x[,i,y] <- xx
    
  }
  
}

# numbers dynamics
x.tot <- apply(x, 2:3, sum)
boxplot(x.tot,
        ylab="N",
        xaxt="n",
        xlab="Year",
        outline=FALSE)
axis(side=1, at=1:nyr.proj)

# numbers dynamics
x[x == 0] <- NA
ggplot() + 
    geom_line(data = melt(apply(x, c(1,3), mean, na.rm = TRUE)), aes(year, value, col = age.class)) + 
    facet_wrap(~age.class)

# survivorship dynamics
ss[ss == 0] <- NA
ggplot() + 
    geom_line(data = melt(apply(ss, c(1,3), mean, na.rm = TRUE)), aes(year, value, col = age.class)) + 
    facet_wrap(~age.class)

# harem size
hh[hh == 0] <- NA
boxplot(hh,
        ylab="harem size",
        xaxt="n",
        xlab="Year",
        outline=FALSE)
axis(side=1, at=1:nyr.proj)
abline(h = 1.14, col = 4)




