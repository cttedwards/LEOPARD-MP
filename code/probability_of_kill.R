
library(reshape2)
library(ggplot2)

rm(list=ls())

# for all age classes
age.structure <- c(0.35,0.2,0.06,0.05,0.04,0.04,0.02,0.1,0.04,0.02,0.01,0.01,0.01,0.05)

# hunter selectivity (Braczkowski et al 2015)
# the probability of shooting a leopard at each age class
hunter.preference <- c(0.256,0.244,0.712,0.462,0.462,0.462,0.359,0.726,0.712,0.897,0.974,0.949,1.00,0.99)

ilogit <- function(x) exp(x) / (exp(x) + 1)

# probablity of kill per day
prob.kill.leopard <- function(total.number = 2000, catchability = 0.005, theta = 0.4) {
    
    # total numbers at age in hunting area
    numbers.at.age <- total.number * age.structure
    
    # assuming independent encounters for each age class
    # the probability of ecounter at age is a function of
    # the numbers at age and the catchability
    prob.encounter.at.age <- ilogit(log(catchability) + theta * log(numbers.at.age))
    
    # probability of not encountering
    prob.not.encounter.at.age <- 1 - prob.encounter.at.age
    
    # the probability of encountering at least one individual
    prob.encounter <- 1 - prod(prob.not.encounter.at.age)
    
    # the probability of a kill is the probability of encouter multiplied
    # by the conditional probability of killing (i.e. the probability that 
    # a hunter will kill an individual of that age given that it has been 
    # encountered)
    prob.kill.at.age <- hunter.preference * prob.encounter.at.age
    
    # probability of not killing
    prob.not.kill.at.age <- 1 - prob.kill.at.age
    
    # the probability of one kill, of any age, is the sum
    # of the probabilities of killing an individual of each
    # age and not killing any of the other ages
    prob.kill <- numeric(length(prob.kill.at.age))
    
    for(age in 1:length(prob.kill))
        prob.kill[age] <- prob.kill.at.age[age] * prod(prob.not.kill.at.age[-age])
    
    prob.kill <- sum(prob.kill)
    
    # return
    return(list(prob.kill = prob.kill, prob.encounter = prob.encounter))
}

# examine influence of catchability

pk <- pe <- nd <- numeric(101)
qv <- seq(0, 0.005, length = 101)
for(i in 1:101) {
    
    out   <- prob.kill.leopard(catchability = qv[i])
    pk[i] <- out$prob.kill
    pe[i] <- out$prob.encounter
    nd[i] <- ifelse(pk[i] == 0, NA, 1/pk[i])
    
}

dfr <- data.frame(catchability = qv, melt(data.frame(prob.kill = pk, prob.encounter = pe, number.days = nd)))
gg <- ggplot(dfr) + geom_line(aes(catchability, value)) + facet_wrap(~variable, scales = 'free_y')
print(gg)
dfr

#plot(pk~tn, xlab='total numbers', ylab='probability of kill', type = 'l')
#plot(pe~tn, xlab='total numbers', ylab='probability of encounter', type = 'l')
#plot(nd~tn, xlab='total numbers', ylab='number of days', type = 'l')

# influence of numbers
pk <- pe <- nd <- numeric(101)
tn <- seq(0, 2000, length = 101)
for(i in 1:101) {
    
    out   <- prob.kill.leopard(total.number = tn[i], catchability = 0.005)
    pk[i] <- out$prob.kill
    pe[i] <- out$prob.encounter
    nd[i] <- ifelse(pk[i] == 0, NA, 1/pk[i])
    
}

dfr <- data.frame(number.leopard = tn, melt(data.frame(prob.kill = pk, prob.encounter = pe, number.days = nd)))
gg <- ggplot(dfr) + geom_line(aes(number.leopard, value)) + facet_wrap(~variable, scales = 'free_y')
print(gg)
dfr

#plot(pk~tn, xlab='total numbers', ylab='probability of kill', type = 'l')
#plot(pe~tn, xlab='total numbers', ylab='probability of encounter', type = 'l')
#plot(nd~tn, xlab='total numbers', ylab='number of days', type = 'l')


