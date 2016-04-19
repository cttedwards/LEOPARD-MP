



library(rstan)

load('../../results/survival_estimates.Rdata')

Sestimates  <- c(S.1   = surv.diff.c[2] * (1/surv.diff.c[1]),       # cub
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

kestimates <- c(br.36 = 1.857142857,                               # birth rate (mother 36-48)     
                br.48 = 1.842105263,                               # birth rate (mother 48-60)              
                br.60 = 1.857142857,                               # birth rate (mother 60-72)               
                br.72 = 1.923076923,                               # birth rate (mother 72-84)               
                br.84 = 1.845070423)                               # birth rate (mother >84) )

Nestimates <- c(nc = 694,
                nj = 410,
                saf = 117,
                f36 = 99,
                f48 = 80,
                f60 = 73,
                f72 = 45,
                f84 = 194,
                sam = 84,
                m36 = 41,
                m48 = 24,
                m60 = 24,
                m72 = 28,
                m84 = 95)

classes <- c('cub', 'juv', 'saf', 'adf', 'sam', 'adm')

# geometric means of rate data
geomean <- function(x) prod(x) ^ (1/length(x))

S.cub <- Sestimates[1]
S.juv <- Sestimates[2]
S.saf <- Sestimates[3]
S.adf <- geomean(Sestimates[4:8])


kills <- data.frame(year = 2008:2015, trophy = c(35, 35, 35, 35, 35, 28, 41, 36), problem_animal = c(97, 54, 85, 40, 58, 35, 60, 60))
kills$total <- kills$trophy + kills$problem

mdl <- stan_model(file = 'leopard.stan')

dat <- list(T = 8, A = 14, S = Sestimates, k = kestimates, h = 1.5, kills = kills$total)

par_init <- function() list(N0 = Nestimates, H = 0.05)

mdl_fit <- sampling(mdl, data = dat, init = par_init, chains = 1, iter = 1)

M <- extract(mdl_fit, pars = 'M')[[1]]
M[1,,]







