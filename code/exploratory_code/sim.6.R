
###########################################################################################
# Load packages and utility functions
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
nreps <- 1000
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
#harvest.rate <- seq(0, 1, length = 101)
harvest.rate <- 0

# setup selectivity object
selectivity <- matrix(data = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,  # all males
                               #0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,  # males ≥ 3
                               #0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  # males ≥ 6
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  # males ≥ 7
                               #0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,  # males & females ≥ 6
                               #0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,  # males & females ≥ 7
                               #0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  # excl dependent young
                               #0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,  # only adults
                               0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1), # all males with females ≥ 7
                      nrow = 14) 

selectivity <- t(selectivity) ; rownames(selectivity) = c("all.males", 
                                                          #"males≥3", 
                                                          #"males≥6", 
                                                          "males≥7", 
                                                          #"male.female≥6",
                                                          #"male.female≥7", 
                                                          #"excl.dep.young", 
                                                          #"only.adults", 
                                                          "all.males.w.females≥7")

# create population size object
population.size <- array(0, dim = c(length(harvest.rate), nreps, 1))

# create extinction prob object
extinction.probability <-  as.data.frame(NA)

total.removed <- data.frame()
#pop.size <- data.frame()

# query number killed
total.harvested    <- array(0, dim = c(length(harvest.rate), nreps, 1))
total.harvested.df <- as.data.frame(NA)

total.population.size <- data.frame()

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
    
    # projection
    for (i in 1:nreps) {
      
      # sample parameters
      param.sample <- params[,i]
      
      # create new object to hold leopard numbers
      # and vital rates
      xx <- leopard(x.initial, param.sample[1:14], param.sample[15:19])
      
      # assign multiplicative maternal effects
      xx@maternal.effect[] <- matrix(maternal.effects, nrow = 2, ncol = 5, byrow = T)
      
      # loop forward over years
      for (y in 1:nyr.proj) {
        
        # correlated deviation in survival: 
        # log-normal with cv = 0.2
        # truncated at 0 and 1
        sdev <- rnorm(1)
        sigma <- sqrt(log(1+0.20^2))
        
        param.sample[1:14] <- exp(log(param.sample[1:14]) +  sigma * sdev - sigma^2/2)
        param.sample[1:14] <- vapply(vapply(param.sample[1:14],
                                            function(x) max(x,0),numeric(1)),function(x) min(x,1),numeric(1))
        
        # create list of sequential hunting scenarios
        removals <- list(trophy = list(rate = harvest.rate, preference = selectivity[z,]), 
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
        
        # -----
        
        population.size[,i,]  <- sum(xx@.Data)
        #total.harvested[,i,]  <- sum(total.removals)
 
        # -----
        
        # step forward
        xx <- transition(xx)
        
        # record numbers
        x[,i,y] <- xx
         
      }
      
      total.population.size <- rbind(total.population.size, mean(population.size))
      
      #extinction.probability <- cbind(extinction.probability, apply(population.size, c(1, 3), prob.ext.func))
      #total.harvested.df     <- cbind(total.harvested.df,     apply(total.harvested, c(1, 3), mean))
      
    }
    
}

x.tot <- apply(x, 2:3, mean)
boxplot(x.tot,
        ylab="Extinction Probability",
        xaxt="n",
        xlab="Year",
        outline=FALSE)
axis(side=1, at=1:nyr.proj)

# save output
#saver(extinction.probability,
#      total.harvested.df,
#      cub.surv.df,
#      prob.inf.df,
#      name = 'model_run_2')
#
#loader('model_run_1')
#loader('model_run_2')


extinction.probability[1] <- NULL

all.males <- extinction.probability[1,1:100]
all.males$group <- rownames(selectivity)[1]

males7 <- extinction.probability[1,101:200]
males7$group <- rownames(selectivity)[2]

all.males.w.females7 <- extinction.probability[1,201:300]
all.males.w.females7$group <- rownames(selectivity)[3]

all.data <- rbind(all.males,
                  males7,
                  all.males.w.females7)

colnames(all.data)[1:100] <- seq(1,100,1)

all.data.melt <- melt(all.data, id = c("group"))

colnames(all.data.melt)[1:3] <- c("Scenario", "Year", "ExtinctionProb")

###########################################################################################
# Plots
###########################################################################################


# plot extinction probability for multiple ages of harvest
cbPalette <- c("#000000", "#999999", "#E69F00") #, "#56B4E9", "#009E73") #, "#F0E442", "#0072B2", "#39e600", "#CC79A7", "#D55E00")

ext.prob.plot <- ggplot(all.data.melt) + 
  geom_point(aes(Year, ExtinctionProb, group = "Scenario")) + 
  #scale_color_manual(values = cbPalette) +
  #facet_wrap("Scenario") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Extinction probability') + 
  labs(y = '', col = 'Year')
ext.prob.plot





















# -----Extinction probability (101)--------

# clean up output and prepare for plotting
extinction.probability[1] <- NULL

all.males             <- extinction.probability[1:101,1:5]  ; all.males$group <- rownames(selectivity)[1]
#males3                <- extinction.probability[1:101,6:10] ; males3$group    <- rownames(selectivity)[2]
#males6                <- extinction.probability[1:101,11:15] ; males6$group    <- rownames(selectivity)[3]
#males7                <- extinction.probability[1:101,16:20] ; males7$group    <- rownames(selectivity)[4]
#male.female6          <- extinction.probability[1:101,21:25] ; male.female6$group  <- rownames(selectivity)[5]
#male.female7          <- extinction.probability[1:101,26:30] ; male.female7$group  <- rownames(selectivity)[6]
#excl.dep.young        <- extinction.probability[1:101,31:35] ; excl.dep.young$group  <- rownames(selectivity)[7]
#only.adults           <- extinction.probability[1:101,36:40] ; only.adults$group  <- rownames(selectivity)[8]
#all.males.w.females7  <- extinction.probability[1:101,41:45] ; all.males.w.females7$group  <- rownames(selectivity)[9]
all.males.w.females7  <- extinction.probability[1:101,6:10] ; all.males.w.females7$group  <- rownames(selectivity)[2]

all.data <- rbind(all.males, 
                  #males3, 
                  #males6, 
                  #males7, 
                  #male.female6, 
                  #male.female7, 
                  #excl.dep.young, 
                  #only.adults, 
                  all.males.w.females7)

colnames(all.data)[1:5] <- c(10, 20, 30, 40, 50)

all.data$H <- c(harvest.rate, 
                #harvest.rate, 
                #harvest.rate, 
                #harvest.rate, 
                #harvest.rate, 
                #harvest.rate, 
                #harvest.rate, 
                #harvest.rate, 
                harvest.rate)

# reshape and change factor levels
all.data.melt <- melt(all.data, id = c("group", "H"))
all.data.melt$group <- factor(all.data.melt$group, levels = c(#"excl.dep.young", 
                                                              #"only.adults",
                                                              "all.males", 
                                                              #"males≥3", 
                                                              #"males≥6", 
                                                              #"male.female≥6",
                                                              #"males≥7",
                                                              #"male.female≥7",
                                                              "all.males.w.females≥7"))

# plot extinction probability for multiple ages of harvest
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73") #, "#F0E442", "#0072B2", "#39e600", "#CC79A7", "#D55E00")
ext.prob.plot <- ggplot(all.data.melt) + 
  geom_line(aes(H, value, color = variable)) + 
  scale_color_manual(values = cbPalette) +
  facet_wrap("group") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Extinction probability') + 
  labs(y = '', col = 'Year')
ext.prob.plot

pdfr(ext.prob.plot, width = 10, name = 'extinction probability')

###########################################################################################
# End
###########################################################################################



