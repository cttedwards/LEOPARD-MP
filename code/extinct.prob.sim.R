
library(leopard)
library(reshape2)
library(ggplot2)

# survival rates per demographic category
ss <- c(0.3270833, 0.7197452, 0.9615385, 0.8882435, 0.9729730, 0.9382353, 0.9230769, 0.7219583, 0.9124424, 0.9642857, 1.0000000, 1.0000000, 0.9000000, 0.2857143)

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
niter <- 10

# harvest rate
harvest.rate <- seq(0, 0.6, length = 101)

# setup selectivity object
selectivity <- matrix(data = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,  # all males
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,  # males ≥ 3
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  # males ≥ 6
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  # males ≥ 7
                               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  # all ages
                               0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,  # males and females ≥ 7
                               0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  # excl dependent young
                               0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,  # only adults
                               0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1), # all males with females ≥ 7
                      nrow = 14) 

selectivity <- t(selectivity) ; rownames(selectivity) = c("all.males", "males≥3", "males≥6", "males≥7", "all.ages",
                                                          "male.female≥7", "excl.dep.young", "only.adults", "all.males.w.females≥7")

# create extinction probability object
prob.ext <- array(0, dim = c(length(harvest.rate), niter, 10))
# create final extinction prob object
all.prob.ext <- as.data.frame(NA)


# loop over range of selectivities
for(z in 1:nrow(selectivity)){
  
  # loop over harvest rates
  for (i in 1:length(harvest.rate)) {
    
    # create new leopard object
    ss.harvest <- ss * (1 - selectivity[z,] * harvest.rate[i])
    xx <- leopard(x.initial, survival.rates = ss.harvest, litter.sizes = ll, harem.size.min = 1.5)
    
    # loop over iterations
    for (j in 1:niter) {
      
      # re-initialise
      yy <- xx
      
      # project forward to adjust numbers
      for (k in 1:100) {
        
        yy <- birth(survival(yy))
        
        # for every 10th year
        if (k %% 10 == 0) {
          # record extinction probability
          prob.ext[i, j, k/10]  <- 1 - if(sum(yy@.Data) / sum(x.initial) > 1){
            1
          } else{
            sum(yy@.Data) / sum(x.initial)
          }
          
        }
        
        yy <- transition(yy)
        
      }
      
    }
    
  }
  
  all.prob.ext <- cbind(all.prob.ext, apply(prob.ext, c(1, 3), median))
  
}

# clean up output and prepare for plotting
all.prob.ext[1] <- NULL

all.males             <- all.prob.ext[1:101,1:10]  ; all.males$group <- rownames(selectivity)[1]
males3                <- all.prob.ext[1:101,11:20] ; males3$group    <- rownames(selectivity)[2]
males6                <- all.prob.ext[1:101,21:30] ; males6$group    <- rownames(selectivity)[3]
males7                <- all.prob.ext[1:101,31:40] ; males7$group    <- rownames(selectivity)[4]
all.ages              <- all.prob.ext[1:101,41:50] ; all.ages$group  <- rownames(selectivity)[5]
male.female7          <- all.prob.ext[1:101,51:60] ; male.female7$group  <- rownames(selectivity)[6]
excl.dep.young        <- all.prob.ext[1:101,61:70] ; excl.dep.young$group  <- rownames(selectivity)[7]
only.adults           <- all.prob.ext[1:101,71:80] ; only.adults$group  <- rownames(selectivity)[8]
all.males.w.females7  <- all.prob.ext[1:101,81:90] ; all.males.w.females7$group  <- rownames(selectivity)[9]

all.data <- rbind(all.males, males3, males6, males7, all.ages, male.female7, excl.dep.young, only.adults, all.males.w.females7)
colnames(all.data)[1:10] <- seq(10, 100, length = 10)
all.data$H <- c(harvest.rate, harvest.rate, harvest.rate, harvest.rate, harvest.rate, harvest.rate, harvest.rate, harvest.rate, harvest.rate)

# reshape and change factor levels
all.data.melt <- melt(all.data, id = c("group", "H"))
all.data.melt$group <- factor(all.data.melt$group, levels = c("all.ages", "excl.dep.young", "only.adults",
                                                              "all.males", "all.males.w.females≥7", 
                                                              "males≥3", "males≥6", "males≥7",
                                                              "male.female≥7"))

# plot extinction probability for multiple ages of harvest
ggplot(all.data.melt) + geom_line(aes(H, value, col = as.factor(variable))) + 
  facet_wrap("group") +
  ggtitle('Extinction probability') + labs(y = '', col = 'Year')

