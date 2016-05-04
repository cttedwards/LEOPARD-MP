


library(plyr)
library(rstan)

load('../../results/survival_estimates.Rdata')

S_estimates  <- c(S.1   = surv.diff.c[2] * (1/surv.diff.c[1]),       # cub
                 S.2   = surv.diff.j[2] * (1/surv.diff.j[1]),       # juvenile
                 S.3   = surv.diff.f[4] * (1/surv.diff.f[3]),       # sub-adult female
                 S.4   = surv.diff.f[5] * (1/surv.diff.f[4]),       # adult female (36-48)
                 S.5   = surv.diff.f[6] * (1/surv.diff.f[5]),       # adult female (48-60)
                 S.6   = surv.diff.f[7] * (1/surv.diff.f[6]),       # adult female (60-72)
                 S.7   = surv.diff.f[8] * (1/surv.diff.f[7]),       # adult female (72-84)
                 S.8   = surv.diff.f[9] * (1/surv.diff.f[8]),       # adult female (>84)
                 S.9   = surv.diff.m[4] * (1/surv.diff.m[3]),       # sub-adult male
                 S.10  = surv.diff.m[5] * (1/surv.diff.m[4]),       # adult male (36-48)
                 S.11  = surv.diff.m[6] * (1/surv.diff.m[5]),       # adult male (48-60)
                 S.12  = surv.diff.m[7] * (1/surv.diff.m[6]),       # adult male (60-72)
                 S.13  = surv.diff.m[8] * (1/surv.diff.m[7]),       # adult male (72-84)
                 S.14  = surv.diff.m[9] * (1/surv.diff.m[8]))       # adult male (>84)

N_estimates <- c('2008' = NA, '2009' = NA, '2010' = NA, '2011' = NA, '2012' = NA, '2013' = 2008, '2014' = NA, '2015' = NA) 
N_estimates[is.na(N_estimates)] <- -1

classes <- c('cub', 'juv', 'saf', 'adf', 'sam', 'adm')

# geometric means of rate data
geomean <- function(x) prod(x) ^ (1/length(x))

S_estimates <- c('S.cub' = S_estimates[1], 'S.juv' = S_estimates[2], 'S.saf' = S_estimates[3], 'S.adf' = geomean(S_estimates[4:8]), 'S.sam' = S_estimates[9], 'S.adm' = geomean(S_estimates[10:14]))

# kill data
kills <- data.frame(year = 2008:2015, 
                    trophy = c(35, 35, 35, 35, 35, 28, 41, 36), 
                    problem_animal = c(97, 54, 85, 40, 58, 35, 60, 60))
kills$total <- kills$trophy + kills$problem

# proportions data
age_comp <- read.csv('age_class_composition.csv', row.names = 1)
age_comp <- as.matrix(t(age_comp))

# abundance
abundance <- read.csv('density_estimates.csv')
abundance$density[is.na(abundance$density)] <- -1

# compile
mdl <- stan_model(file = 'leopard.stan')

dat <- list(T = 8, A = 6, S = as.numeric(S_estimates), k = 2, proportions = as.matrix(age_comp), kills = kills$total, density = as.numeric(abundance[,2]), numbers = as.numeric(N_estimates))

par_init <- function() list(N0 = c(600, 400, 100, 400, 100, 400), H = runif(1, 0, 0.5), logq = runif(1, -8, -3), h = runif(1, 1, 5), selectivity = c(0.1, 0.1, 0.9, 0.9, 0.9, 0.9))

mdl_fit <- sampling(mdl, data = dat, init = par_init, chains = 10, iter = 1e4, thin = 10)

#Warning messages:
#1: There were 4217 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. 
#2: There were 2131 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. 
#3: Examine the pairs() plot to diagnose sampling problems

# STILL NOT FITTING TO DENSITY DATA.... 

traceplot(mdl_fit, pars = 'h', inc_warmup = TRUE)
traceplot(mdl_fit, pars = 'H', inc_warmup = TRUE)
traceplot(mdl_fit, pars = 'logq', inc_warmup = TRUE)
traceplot(mdl_fit, pars = 'N0')
traceplot(mdl_fit, pars = 'selectivity')

mdl_res <- extract(mdl_fit)

hist(mdl_res[['H']],, main = 'Harvest rate')
hist(mdl_res[['logq']], main = 'log(catchability)')
hist(apply(mdl_res[['N0']], 1, sum), main = 'Total numbers (initial)')
hist(mdl_res[['h']], main = 'Harem size')

kills$totalHat <- apply(mdl_res[['killsHat']], 2, median)
plot(total ~ year, kills)
lines(totalHat ~ year, kills)

numbers <- data.frame(year = 2008:2015, N = N_estimates, Nhat = apply(mdl_res[['numbersHat']], 2, median))
numbers[numbers < 0] <- NA
plot(N ~ year, numbers)
lines(Nhat ~ year, numbers)

abundance$densityHat <- apply(mdl_res[['densityHat']], 2, median)
abundance[abundance < 0] <- NA
plot(density ~ year, abundance)
lines(densityHat ~ year, abundance)

thetaHat <- apply(mdl_res[['thetaHat']], 2:3, median)

dfr <- data.frame(age = colnames(age_comp), theta = age_comp[6,], thetaHat = thetaHat[6, ] / sum(thetaHat[6,]))
plot(theta ~ age, dfr)
points(thetaHat ~ age, dfr, col = 2)

dfr <- data.frame(age = colnames(age_comp), theta = age_comp[7,], thetaHat = thetaHat[7, ] / sum(thetaHat[7,]))
plot(theta ~ age, dfr)
points(thetaHat ~ age, dfr, col = 2)

dfr <- data.frame(age = colnames(age_comp), theta = age_comp[8,], thetaHat = thetaHat[8, ] / sum(thetaHat[8,]))
plot(theta ~ age, dfr)
points(thetaHat ~ age, dfr, col = 2)

