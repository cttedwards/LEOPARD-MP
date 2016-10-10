
library(reshape2)
library(ggplot2)

rm(list = ls())

# simple exploration

# binomial probability i.e. same p for each age class

p <- seq(0.01, 0.99, length = 100)

# at higher densities the probability of killing only one individual decreases
plot(p,dbinom(1, 14, p), type = 'l', ylab = 'binomial probability', xlab = 'probability of killing one')
#lines(p,dbinom(2, 14, p))
#lines(p,dbinom(4, 14, p))
#lines(p,dbinom(6, 14, p))

# but the probability of killing at least one individual continues to increase
plot(p, 1 - pbinom(1, 14, p), type = 'l', ylab = 'binomial probability', xlab = 'probability of killing from at least one age class')
#lines(p,1 - pbinom(2, 14, p))
#lines(p,1 - pbinom(4, 14, p))
#lines(p,1 - pbinom(6, 14, p))



# probability of kill for 14 age classes
prob.kill.at.age <- rbeta(14, 20, 2)

# probability of not kill
prob.not.kill.at.age <- 1 - prob.kill.at.age

# probability of killing at age and not killing all other ages
binomial.prob.kill.at.age <- numeric(14)

for (i in 1:14) 
    binomial.prob.kill.at.age[i] <- prob.kill.at.age[i] * prod(prob.not.kill.at.age[-i])

# probability of killing one individual of any age on a given day
prob.kill <- sum(binomial.prob.kill.at.age)

# probabiliy of killing at least one on a given day
#prob.kill <- 1 - prod(prob.not.kill.at.age)

# expected waiting time from a geometric distribution
1 / prob.kill

# more detailed exploration

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
    
    # the probability of at least one kill, of any age
    prob.kill <- 1 - prod(prob.not.kill.at.age)
    
    # return
    return(list(prob.kill = prob.kill, prob.encounter = prob.encounter))
}

# influence of numbers
pk <- pe <- nd <- numeric(101)
tn <- seq(0, 2000, length = 101)
for (i in 1:101) {
    
    out   <- prob.kill.leopard(total.number = tn[i], catchability = 0.007)
    pk[i] <- out$prob.kill
    pe[i] <- out$prob.encounter
    nd[i] <- ifelse(pk[i] == 0, NA, 1/pk[i])
    
}

dfr <- data.frame(number.leopard = tn, melt(data.frame(prob.kill = pk, prob.encounter = pe, number.days = nd)))
gg <- ggplot(dfr) + geom_line(aes(number.leopard, value)) + facet_wrap(~variable, scales = 'free_y')
print(gg)
dfr

# alternative probablity of kill per day
prob.kill.leopard <- function(total.number = 2000, catchability = 0.005, theta = 0.4) {
    
    # total numbers at age in hunting area
    numbers.at.age <- total.number * age.structure
    
    # assuming independent encounters for each age class
    # the probability of ecounter at age is a function of
    # the numbers at age and the catchability
    prob.encounter.at.age <- ilogit(log(catchability) + theta * log(numbers.at.age))
    
    # the probability of a kill is the probability of encouter multiplied
    # by the conditional probability of killing (i.e. the probability that 
    # a hunter will kill an individual of that age given that it has been 
    # encountered)
    prob.kill.at.age <- hunter.preference * prob.encounter.at.age
    
    # for each age calculate the expected waiting time
    waiting.time <- numeric(length(prob.kill.at.age))
    
    for (age in 1:length(waiting.time))
        waiting.time[age] <- 1 / prob.kill.at.age[age]
    
    # realised waiting time is the minimum across ages
    waiting.time <- min(waiting.time)
    
    # daily probability of a kill (for comparative purposes only)
    prob.kill <- 1/waiting.time
    
    # return
    return(list(prob.kill = prob.kill, waiting.time = waiting.time))
}

# influence of numbers
pk <- pe <- nd <- numeric(101)
tn <- seq(0, 2000, length = 101)
for (i in 1:101) {
    
    out   <- prob.kill.leopard(total.number = tn[i], catchability = 0.007)
    pk[i] <- out$prob.kill
    nd[i] <- out$waiting.time
    
}

dfr <- data.frame(number.leopard = tn, melt(data.frame(prob.kill = pk, number.days = nd)))
gg <- ggplot(dfr) + geom_line(aes(number.leopard, value)) + facet_wrap(~variable, scales = 'free_y')
print(gg)

# examine influence of catchability

pk <- pe <- nd <- numeric(101)
qv <- seq(0, 0.005, length = 101)
for(i in 1:101) {
    
    out   <- prob.kill.leopard(catchability = qv[i])
    pk[i] <- out$prob.kill
    nd[i] <- out$waiting.time
    
}

dfr <- data.frame(catchability = qv, melt(data.frame(prob.kill = pk, number.days = nd)))
gg <- ggplot(dfr) + geom_line(aes(catchability, value)) + facet_wrap(~variable, scales = 'free_y')
print(gg)

# test library
library(leopard)


# for all age classes
age.structure <- c(0.35,0.2,0.06,0.05,0.04,0.04,0.02,0.1,0.04,0.02,0.01,0.01,0.01,0.05)

xx <- newleopard(2000 * age.structure)


for (i in 1:1000) {
    xx <- newleopard(rbeta(1,  20, 2) * 2000 * age.structure)
    days[i] <- kill(xx)$waiting.time
}

# geometric distribution
hist(days)


tmp <- numeric(50)
days <- numeric(1000)

for (i in 1:1000) {
    for (j in 1:50) {
        tmp[j] <- kill(xx)$waiting.time
    }
    days[i] <- mean(tmp)
}

# normal distribution (from CLT)
hist(days)
    


