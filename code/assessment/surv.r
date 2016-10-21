
source("utils/reader.r")

# get data
sabi.surv <- reader("sabi.surv")

dat <- data.frame(time1 = sabi.surv$AGE.START.MONTHS, time2 = sabi.surv$AGE.DEATH.MONTHS, event = sabi.surv$CENSUS.STATUS)

# convert time to ages
dat$age1 <- dat$time1 / 12
dat$age2 <- dat$time2 / 12

dat$age3 <- cut(dat$age2, breaks = 0:20)

#dat <- subset(dat, age2 > 0.5)

# exponential survival function
ll <- function(lambda, d, t) {
    logS <- t * lambda
    logh <- log(lambda)
    return(- d * logh - logS)
}

obj <- function(par) {
    sum(apply(with(dat, data.frame(d = event, t = age2)), 1, function(x) ll(par, x['d'], x['t'])))
}

# check profile
plot(seq(0.0001, 10, length.out = 100), vapply(seq(0.0001, 10, length.out = 100), obj, numeric(1)), type = 'l')

# estimate
lambda.hat <- optimize(obj, interval = c(0.0001, 2))$min
abline(v = lambda.hat, col = 2)

# fit to data
hist(dat$age2, prob = TRUE)
curve(dexp(x, lambda.hat), from = 0, to = 20, add = TRUE, col = 2)

# weibull survival function
ll <- function(lambda, p, d, t) {
    logS <- - (lambda * t) ^ p
    logh <- p * log(lambda) + log(p) + (t - 1) * log(t)
    return(- d * logh - logS)
}

obj <- function(par) {
    sum(apply(with(dat, data.frame(d = event, t = age2)), 1, function(x) ll(par[1], par[2], x['d'], x['t'])))
}

# estimate
par.hat <- optim(par = c(0.5, 1), obj, lower = c(0.001, 0.001), upper = c(10, 10), method = "L-BFGS-B")$par

# check profiles
plot(seq(0.0001, 2, length.out = 100), vapply(seq(0.0001, 2, length.out = 100), function(x) obj(c(x, par.hat[2])), numeric(1)), type = 'l', main = "lambda")
abline(v = par.hat[1], col = 2)
plot(seq(0.0001, 2, length.out = 100), vapply(seq(0.0001, 2, length.out = 100), function(x) obj(c(par.hat[1], x)), numeric(1)), type = 'l', main = "p")
abline(v = par.hat[2], col = 2)

# fit to data
hist(dat$age2, prob = TRUE, breaks = 40)
curve(dweibull(x, par.hat[1], par.hat[2]), from = 1, to = 20, add = TRUE, col = 2)

# plot hazard function
h <- function(x, lambda, p) lambda ^ p * p * x ^ (p - 1)
plot(1:20, h(1:20, par.hat[1], par.hat[2]), type = 'l', ylim = c(0,1))

# plot survival probability
S <- function(x, lambda, p) exp(- (lambda * x) ^ p)
plot(1:20, S(1:20, par.hat[1], par.hat[2]), type = 'l', ylim = c(0,1))

# compare to KM estimate
sabi.surv.fit <- survfit(Surv(time=age1,time2=age2,event=event,type='counting')~1, data = dat)
plot(sabi.surv.fit, xlab="AGE (YEARS)", ylab="SURVIVAL PROBABILITY", conf.int = TRUE)
lines(0:20, S(0:20, par.hat[1], par.hat[2]), col = 4)


