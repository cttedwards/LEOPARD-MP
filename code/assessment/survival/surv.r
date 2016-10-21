
# get data
sabi.surv <- read.csv("C:/RESEARCH/LEOPARD-MP/data/sabi.surv.csv")

dat <- data.frame(time1 = sabi.surv$AGE.START.MONTHS, time2 = sabi.surv$AGE.DEATH.MONTHS, event = sabi.surv$CENSUS.STATUS)

# convert time to ages
dat$age1 <- dat$time1 / 12
dat$age2 <- dat$time2 / 12

dat$age3 <- cut(dat$age2, breaks = 0:20)

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
plot(seq(0.001, 2, length.out = 100), vapply(seq(0.001, 2, length.out = 100), obj, numeric(1)), type = 'l')

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
curve(dweibull(x, par.hat[1], par.hat[2]), from = 0, to = 20, add = TRUE, col = 2)

# plot hazard function
h <- function(x, lambda, p) lambda ^ p * p * x ^ (p - 1)
plot(0:20, h(0:20, par.hat[1], par.hat[2]), type = 'l', ylim = c(0,1))

# plot survival probability
S <- function(x, lambda, p) exp(- (lambda * x) ^ p)
plot(1:20, S(1:20, par.hat[1], par.hat[2]), type = 'l', ylim = c(0,1))

# compare to KM estimate
sabi.surv.fit <- survival::survfit(survival::Surv(time=age1,time2=age2,event=event,type='counting')~1, data = dat)
plot(sabi.surv.fit, xlab="AGE (YEARS)", ylab="SURVIVAL PROBABILITY", conf.int = TRUE)
lines(0:20, S(0:20, par.hat[1], par.hat[2]), col = 4)

# check linearity
dfr <- with(sabi.surv.fit, data.frame(surv = surv, llsurv = log(-log(surv)), time = time, ltime = log(time)))
dfr <- subset(dfr, surv > 0 & surv < 1)
plot(llsurv ~ time, dfr)
abline(lm(llsurv ~ time, subset(dfr, time >= 1)))
abline(lm(llsurv ~ time, subset(dfr, time < 1)))

# mixture weibull model
S     <- expression(pi + (1 - pi) * exp(-(lambda * t)^p))
logS  <- expression(log(pi + (1 - pi) * exp(-(lambda * t)^p)))
nlogS <- expression(-log(pi + (1 - pi) * exp(-(lambda * t)^p)))
h     <- D(nlogS, 't')

ll <- function(lambda, d, t, p = 1, pi = 0) {
    x1 <- eval(logS)
    x2 <- log(eval(h))
    return(-d * x2 - x1)
}

# check with exponential distribution fit
obj <- function(par) {
    sum(apply(with(dat, data.frame(d = event, t = age2)), 1, function(x) ll(par, x['d'], x['t'])))
}

# check profile estimate
plot(seq(0.001, 2, length.out = 100), vapply(seq(0.001, 2, length.out = 100), obj, numeric(1)), type = 'l')
lambda.hat <- optimize(obj, interval = c(0.0001, 2))$min
abline(v = lambda.hat, col = 2)

# now estimates p and pi
obj <- function(par) {
    sum(apply(with(dat, data.frame(d = event, t = age2)), 1, function(x) ll(par[1], x['d'], x['t'], par[2], par[3])))
}

# estimate
par.hat <- optim(par = c(0.5, 1, 0.01), obj, lower = c(0.001, 0.001, 0.001), upper = c(3, 3, 0.999), method = "L-BFGS-B")$par

# check marginal profiles
plot(seq(0.001, 3, length.out = 100), vapply(seq(0.001, 3, length.out = 100), function(x) obj(c(x, par.hat[2], par.hat[3])), numeric(1)), type = 'l', main = "lambda")
abline(v = par.hat[1], col = 2)
plot(seq(0.001, 3, length.out = 100), vapply(seq(0.001, 3, length.out = 100), function(x) obj(c(par.hat[1], x, par.hat[3])), numeric(1)), type = 'l', main = "p")
abline(v = par.hat[2], col = 2)
plot(seq(0.001, 0.999, length.out = 100), vapply(seq(0.001, 0.999, length.out = 100), function(x) obj(c(par.hat[1], par.hat[2], x)), numeric(1)), type = 'l', main = "pi")
abline(v = par.hat[3], col = 2)


# fit to data
hist(dat$age2, prob = TRUE, breaks = 40)
curve(dweibull(x, par.hat[1], par.hat[2]), from = 0, to = 20, add = TRUE, col = 2)

# plot survival and hazard function
windows(width = 14); par(mfrow = c(1, 2))
tt <- seq(0, 20, length = 1000)
hh <- vapply(tt, function(x) { t <- x; lambda <- par.hat[1]; p <- par.hat[2]; pi <- par.hat[3]; eval(h) }, numeric(1))
plot(tt, hh, type = 'l', main = "Hazard")
SS <- vapply(tt, function(x) { t <- x; lambda <- par.hat[1]; p <- par.hat[2]; pi <- par.hat[3]; eval(S) }, numeric(1))
plot(tt, SS, type = 'l', main = "Survival", ylim = c(0, 1))

# compare to KM estimate
sabi.surv.fit <- survival::survfit(survival::Surv(time=age1,time2=age2,event=event,type='counting')~1, data = dat)
plot(sabi.surv.fit, xlab="AGE (YEARS)", ylab="SURVIVAL PROBABILITY", conf.int = TRUE)
lines(tt, SS, col = 4)

