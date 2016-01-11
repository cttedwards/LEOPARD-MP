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
# and distribution of wating times
n     <- 2000
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
hist(sim, prob = TRUE, main = 'increasing prob.')

# concave catchability
#catchability2 <- 1/(10000 - 9995/99 * (1:100 - 1))
catchability2 <- seq(0.01,2, length.out = 100) ^ 2

plot(catchability)
lines(catchability2)

logitpp2 <- log(catchability2) + theta * log(n)

pp2 <- ilogit(logitpp2)

plot(1:length(pp2), pp2, ylim = c(0,1), xlab = 'sequential days', ylab = 'prob. encounter', type = 'l')

plot(catchability2, pp2, ylim = c(0,1), xlab = 'catchability', ylab = 'prob. encounter', type = 'l')

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp2)
hist(sim, prob = TRUE)

###########
# 3. Simulate waiting time data

# probablity of kill per day since hunting began
prob.encounter <- function(d, beta = c(0.01, 0.02, 2), n = 2000, theta = 0.1) {
    
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
pp <- prob.encounter(1:100, beta = c(0.001,0.004,1))[[2]]

plot(pp, type = 'l', ylim = c(0, 1), xlab = 'sequential days', ylab = 'prob. encounter')

sim <- numeric(1000)
for (i in 1:length(sim)) sim[i] <- rdays(pp)
hist(sim, prob = TRUE, main = 'days')

# CHECK LIBRARY
library(leopard)

age.structure <- c(0.35,0.2,0.06,0.05,0.04,0.04,0.02,0.1,0.04,0.02,0.01,0.01,0.01,0.05)

xx <- newleopard(2000 * age.structure)

days <- numeric(1000)
for (i in 1:1000) {
    xx <- newleopard(rbeta(1,  20, 2) * 2000 * age.structure)
    days[i] <- kill(xx,beta = c(0.0001, 0.0004, 1))$waiting.time
}

# geometric distribution
hist(days)





