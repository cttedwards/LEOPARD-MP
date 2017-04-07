

library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("brooks.rda")

y[is.na(y)] <- -1

dat <- list(T = T, T1 = T1, T2 = T2, f = f, y = y, m = m)


mdl <- stan_model(stanc_ret = stanc(file = "brooks.stan", model_name = "lapwing"))


ini <- function() return(list(sigy = 1, 
            Na = c(1000.,1000,1092.23,1100.01,1234.32,1460.85,1570.38,1819.79,
                       1391.27,1507.60,1541.44,1631.21, 1628.60,1609.33,1801.68,1809.08,1754.74,
                       1779.48,1699.13,1681.39,1610.46,1918.45,1717.07,1415.69, 1229.02,1082.02,
                       1096.61,1045.84,1137.03,981.1,647.67,992.65,968.62,926.83,952.96,865.64),
      N1 = c(400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,
             400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400), 
      alpha1 = 1, alphaa = 2, alphar = -2, alphal = -4, beta1 = -2, betaa = 0.1, betar = -0.7,betal = -0.3))

ini <- function() list(sigy = 1, 
            Na = vapply(rnorm(36, 1000, 10), function(x) max(x, 1), numeric(1)),
            N1 = vapply(rnorm(36, 400, 10), function(x) max(x, 1), numeric(1)),
            alpha1 = rnorm(1,1,1), 
            alphaa = rnorm(1,2,1), 
            alphar = rnorm(1,-2,1), 
            alphal = rnorm(1,-4,1), 
            beta1 = rnorm(1,-2,1), 
            betaa = rnorm(1,0.1,1), 
            betar = rnorm(1,-0.7,1),
            betal = rnorm(1,-0.3,1))


res <- sampling(mdl, init = ini, data = dat, chains = 10)

traceplot(res, pars = "Na", inc_warmup = TRUE)
