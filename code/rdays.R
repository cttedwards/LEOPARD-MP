##########################################
# RANDOM GENERATION OF WAITING TIME DATA #
# (TIME UNTIL FIRST KILL)                #
##########################################

###########
# 1. Explore relationship between CHANGING PROBABILITY
# of encounter and distribution of waiting times

# simulate random days
rdays <- function(p) {
    if (length(p) > 1) {
        for (i in 1:length(p)) {
            if (runif(1) < p[i]) {
                days <- i
                break
            }
        }
    } else {
        days <- 0
        while (runif(1) > p) days <- days + 1
    }
    return(days) 
}


# increasing probability over time
pp <- seq(0.01, 0.99, length.out = 30)

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp)
hist(sim, prob = TRUE, main = 'increasing prob.')

# constant probability over time
pp <- 0.2

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp)
hist(sim, prob = TRUE, main = 'constant prob.')

###########
# 2. Explore relationship between CHANGING CATCHABILITY
# and distribution of waiting times
n     <- 200
theta <- 0.1

ilogit <- function(x) exp(x) / (exp(x) + 1)

# linear catchability
catchability <- seq(0.001,0.4, length.out = 100)

logitpp <- log(catchability) + theta * log(n)

pp <- ilogit(logitpp)

plot(1:length(pp), pp, ylim = c(0,1), xlab = 'sequential days', ylab = 'prob. encounter', type = 'l')

plot(catchability, pp, ylim = c(0,1), xlab = 'catchability', ylab = 'prob. encounter', type = 'l')

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp)
hist(sim, prob = TRUE, main = 'linear increasing prob.')

# constant catchability
catchability2 <- 0.2

logitpp <- log(catchability2) + theta * log(n)

pp2 <- ilogit(logitpp)

plot(1:length(pp2), pp2, ylim = c(0,1), xlab = 'sequential days', ylab = 'prob. encounter', type = 'p')

plot(catchability2, pp2, ylim = c(0,1), xlab = 'catchability', ylab = 'prob. encounter', type = 'p')

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp)
hist(sim, prob = TRUE, main = 'constant prob.')


# concave catchability
#catchability3 <- 1/(10000 - 9995/99 * (1:100 - 1))
catchability3 <- seq(0.01,2, length.out = 100) ^ 2

#plot(catchability)
#lines(catchability3)

logitpp <- log(catchability3) + theta * log(n)

pp3 <- ilogit(logitpp)

plot(1:length(pp), pp, ylim = c(0,1), xlab = 'sequential days', ylab = 'prob. encounter', type = 'l')
lines(1:length(pp3), pp3, col = 2)

plot(catchability, pp, ylim = c(0,1), xlab = 'catchability', ylab = 'prob. encounter', type = 'l')
lines(catchability3, pp3, col = 2)

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp3)
hist(sim, prob = TRUE, , main = 'non-linear increasing prob.')

###########
# 3. Simulate waiting time data

# probablity of kill per day since hunting began
prob.encounter <- function(d, beta = c(0.01, 0.02, 2), n = 200, theta = 0.1) {
    
    ilogit <- function(x) exp(x) / (exp(x) + 1)
    
    catchability <-  (beta[1] + (d - 1) * beta[2]) ^ beta[3]

    return(list(catchability = catchability, prob = ilogit(log(catchability) + theta * log(n))))
}

# simulate random days
rdays <- function(p) {
    success <- FALSE
    for (i in 1:length(p)) {
        if (runif(1) < p[i]) {
            days <- i
            success <- TRUE
            break
        }
    }
    ifelse(success, return(days), return(NA)) 
}

# non-linear increase in catchability over time
pp <- prob.encounter(1:100)[[2]]

plot(pp, type = 'l', ylim = c(0, 1), xlab = 'sequential days', ylab = 'prob. encounter')

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp)
hist(sim, prob = TRUE, main = 'days')

# or linear increase in catchability 
pp <- prob.encounter(1:100, beta = c(0.001,0.004,1), n = 200, theta = 0.1)[[2]]

plot(pp * 0.3, type = 'l', ylim = c(0, 1), xlab = 'sequential days', ylab = 'prob. encounter')

sim <- numeric(1000)
for (i in 1:length(sim)) {
    tmp <- c()
    for (j in 1:20) tmp <- c(tmp, rdays(pp * 0.3))
    sim[i] <- min(tmp) #rdays(pp * 0.2)
}
hist(sim, prob = TRUE, main = 'days')

# or constant catchability 
pp <- prob.encounter(1:200, beta = c(0.05,0.00,1), n = 200, theta = 0.2)[[2]]

plot(pp * 0.3, type = 'l', ylim = c(0, 1), xlab = 'sequential days', ylab = 'prob. encounter')

sim <- numeric(1000)
for (i in 1:length(sim)) {
    tmp <- c()
    for (j in 1:20) tmp <- c(tmp, rgeom(1, pp * 0.3)) #tmp <- c(tmp, rdays(pp * 0.3))
    sim[i] <- sample(tmp, 1) #rdays(pp * 0.2)
}
hist(sim, prob = TRUE, main = 'days')

# success rate
1 - mean(is.na(sim))


# CHECK LIBRARY
library(leopard)

age.structure <- c(0.35,0.2,0.06,0.05,0.04,0.04,0.02,0.1,0.04,0.02,0.01,0.01,0.01,0.05)

xx <- newleopard(2000 * age.structure)

# approximate negative binomial distribution
days <- numeric(1000)
for (i in 1:1000) {
    days[i] <- kill(xx,beta = c(0.0001, 0.0004, 1))$waiting.time
}

hist(days, main = 'linear increase in catchability')

# geometric distribution
days <- numeric(1000)
for (i in 1:1000) {
    days[i] <- kill(xx, beta = c(0.1, 0, 1), theta = 0.4)$waiting.time
}

hist(days, main = 'constant catchability')



#########################################################################
# the waiting time will decrease as the number of age classes increases #
# => need alternative formulation
# => population age/day matrix using conditional probabilities

# probablity of kill per day since hunting began
prob.kill <- function(d = 1:100, beta = c(0.1, 0.0, 1), n = 200, theta = 0.4, preference = 0.3) {
    
    ilogit <- function(x) exp(x) / (exp(x) + 1)
    
    catchability <-  (beta[1] + (d - 1) * beta[2]) ^ beta[3]
    
    prob.enc <- ilogit(log(catchability) + theta * log(n))
    
    prob.kill <- prob.enc * preference
    
    return(prob.kill)
}

age.structure <- c(0.35,0.2,0.06,0.05,0.04,0.04,0.02,0.1,0.04,0.02,0.01,0.01,0.01,0.05)

n <- 2000 * age.structure

prob.kill.at.age <- matrix(0, 14, 100)
for (i in 1:14)
    prob.kill.at.age[i,] <- prob.kill(n = n[i], beta = c(0.001, 0.001, 1), theta = 1)

prob.kill.at.age.at.time <- matrix(0, 14, 100)
prob.kill.at.age.at.time[,1] <- prob.kill.at.age[,1]
for (i in 2:100)
    prob.kill.at.age.at.time[,i] <- prob.kill.at.age[,i] * prod(1 - prob.kill.at.age[,1:(i-1)])

plot(prob.kill.at.age.at.time[1,], type = 'l', ylim = c(0,0.2))
for (i in 2:14) lines(prob.kill.at.age.at.time[i,], col = i)
