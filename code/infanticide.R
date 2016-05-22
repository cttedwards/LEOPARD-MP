

library(ggplot2)
library(plyr)

source('utils/pdfr.r')

# numbers (mature)
# from Swanepoel et al. 2014
nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

nfemales
nmales
hsize <- 1.14
nmales * hsize
nfemales / (nmales * hsize)
pfbreed <- (nmales * hsize) / nfemales
sratio <- nfemales / nmales
# => at the current sex ratio only half the 
# => females are in harems assuming that all 
# => males breed

# the probability that a female breeds follows a binomial distribution
# with the number of trials equal to the number of females and the probability
# of success calculated so that nfemales * psuccess -> nmales * haremsize
# currently we have 491 * psuccess = 212 * 1.14, giving psuccess = 0.492

# we can write the relationship between psuccess and the sex ratio as a binomial
# regression and fit to the known values to give the psuccess and harem size

pfbreed.ff <- function(mu, beta = 1, sexratio = sratio) (exp(-beta * sexratio - mu)) / (exp(-beta * sexratio - mu) + 1)

obj <- function(mu) {
    pfbreed - pfbreed.ff(mu)
}

muhat <- uniroot(obj, interval = c(-10, 10))$root

# checks
pfbreed.ff(mu = muhat) # == pfbreed = 0.492
sratio * pfbreed.ff(mu = muhat) # == hsize = 1.14

# when equal number of males and females
pfbreed.ff(mu = muhat, sexratio = 1)
1 * pfbreed.ff(mu = muhat, sexratio = 1)

pfbreed.ff(mu = muhat, sexratio = 2)
2 * pfbreed.ff(mu = muhat, sexratio = 2)

pfbreed.ff(mu = muhat, sexratio = 3)
3 * pfbreed.ff(mu = muhat, sexratio = 3)

# as the sex ratio changes we expect psuccess to 
# also change and we can use this to back calculate
# the harem size at different sex ration values

sr <- seq(0.5, 3, length = 101)
plot(sr, pfbreed.ff(muhat, sexratio = sr), main = 'Prob. that a female breeds \n(i.e. is in a harem)', ylab = '', xlab = 'sex ratio (N females / N males)')
abline(v = sratio); abline(h = pfbreed)

plot(sr, sr * pfbreed.ff(muhat, sexratio = sr), main = 'Harem size', ylab = '', xlab = 'sex ratio (N females / N males)')
abline(v = sratio); abline(h = hsize)

# create orthogonal data.frame
#nmales   <- seq(0, 600, by = 100)#nmales * seq(0, 2, length = 11)
nfemales <- seq(0, 1000, length = length(sr))

dfr <- expand.grid(nf = nfemales, sr = sr)
dfr$nm   <- dfr$nf / dfr$sr

dfr$pfb <- apply(dfr, 1, function(x) pfbreed.ff(mu = muhat, sexratio = x['sr']))
dfr$hs  <- apply(dfr, 1, function(x) x['sr'] * pfbreed.ff(mu = muhat, sexratio = x['sr']))

dfr$nfb <- dfr$nf * dfr$pfb
dfr$nmb <- dfr$nm

ggplot(dfr) + geom_line(aes(sr, pfb)) + geom_line(aes(sr, hs)) #, col = as.factor(nm))) + facet_wrap(~nf) + ggtitle('number_encounters()')

ggplot(dfr) + geom_line(aes(nm, nfb, col = sr))

# finally we can estimate the rate of increase in number of female
# breeders per male and plot this as a function of the sex ratio
dfr.rho <- ddply(dfr, .(sr), summarize, rho = coef(lm(nfb ~ nm))[2])

plot(rho ~ sr, dfr.rho, main = 'Marginal increase in number \nof females breeding per male', ylab = '', xlab = 'sex ratio (N females / N males)')
abline(v = sratio)

# this could be used to calculate exile probability for males
# => if marginal increase in number of females drops below 1 then males will start to leave
plot(rho ~ I(1/sr), dfr.rho, main = 'Marginal increase in number \nof females breeding per male', ylab = '', xlab = 'sex ratio (N males / N females)')
abline(v = 1/sratio); abline(h = 1, lty = 2)



# create orthogonal data.frame
nmales   <- seq(0, 600, by = 100)#nmales * seq(0, 2, length = 11)
nfemales <- seq(0, 600, by = 100)

dfr <- expand.grid(nm = nmales, nf = nfemales, sm = seq(0, 1, length = 101))
dfr$mk <- (1 - dfr$sm) * dfr$nm
dfr$fk <- 0
#dfr$sm <- dfr$nm / max(dfr$nm)
rm(nmales, nfemales)

# number of encounters between males and females
# follows the harmonic mean function
number_encounters <- function(nmales, nfemales) {
    
    # number of encounters between males and females
    nencounters <- 2 * nmales * nfemales / (nmales + nfemales)
    
    # return
    return(nencounters)
}

dfr$ne <- apply(dfr, 1, function(x) number_encounters(x[1], x[2]))

ggplot(dfr) + geom_line(aes(H <- 1 - sm, ne, col = as.factor(nm))) + facet_wrap(~nf) + ggtitle('number_encounters()')

# harem size
harem_size <- function(nmales, nfemales) {
    
    # females per male
    lambda <- nfemales / nmales
    
    prob_harem <- 1 - ppois(0, lambda)
    #prob_harem <- ifelse(nmales > 0, 1 - pnbinom(0, size = 100, mu = lambda), 0)
    #prob_harem <- pnbinom(lambda, size = nfemales, mu = lambda)
    
    number_harems <- nmales * prob_harem
    
    harem_size <- ifelse(nfemales > 0, nfemales / number_harems, 1)
    
    return(list(hs = harem_size, ehs = lambda))
}

dfr$hs  <- apply(dfr, 1, function(x) harem_size(x[1], x[2])[['hs']])
dfr$ehs <- apply(dfr, 1, function(x) harem_size(x[1], x[2])[['ehs']])

ggplot(dfr) + geom_line(aes(nf / nm, hs)) + 
    geom_line(aes(nf / nm, ehs), col = 'red') + 
    geom_point(aes(x = 491 / 212, y = 1.14), size = 3) +
    ggtitle('harem_size()')


# number of exta encounters between males and females
# as a result of hunting
number_extra_encounters <- function(nmales, nfemales, mkills = 0, fkills = 0) {
    
    nmalesk   <- nmales - mkills
    nfemalesk <- nfemales - fkills
    
    # number of extra encounters between males and females
    #nencounters <- 2 * (nmales * nfemales / (nmales + nfemales) - nmalesk * nfemalesk / (nmalesk + nfemalesk))
    nencounters <- 2 * nmalesk * nfemalesk / (nmalesk + nfemalesk)
    
    # return
    return(nencounters)
}

dfr$nee <- apply(dfr, 1, function(x) number_extra_encounters(x[1], x[2], x[4], x[5]))

ggplot(dfr) + geom_line(aes(H <- 1 - sm, nee, col = as.factor(nm))) + facet_wrap(~nf) + ggtitle('number_encounters()')



# probability of at least one encounter per female
# follows a Poisson probability
probability_encounter <- function(nmales, nfemales) {
    
    nencounters <- number_encounters(nmales, nfemales)

    nencounters_perfemale <- nencounters / nfemales
    
    pencounter <- 1 - exp(-nencounters_perfemale)
    
    return(pencounter)
}

dfr$pe <- apply(dfr, 1, function(x) probability_encounter(x[1], x[2]))

ggplot(dfr) + geom_line(aes(nm, pe, col = as.factor(nf))) + ggtitle('probability_encounter()')

#gg <- ggplot(dfr) + geom_line(aes(nm, pe, col = as.factor(nf))) + ggtitle('Probability of male occupancy per harem\n') + labs(x = 'Number of mature males', y = '', col = 'Number of mature \nfemales')
#pdfr(gg, width = 10, name = 'pocc')

# cub survivorship is a function of male survival
# and the probability of encounter
cub_survivorship <- function(nmales, nfemales, smales) {
    
    # equilibrium cub survial
    scub_star <- 0.3270833
    
    # probability of mother encountering a male
    pencounter <- probability_encounter(nmales, nfemales)
    
    # probability of infanticide
    pinfanticide <- pencounter * (1 - smales)
    
    scub <- scub_star * (1 - pinfanticide)
    
    return(scub)
    
}

dfr$sc <- apply(dfr, 1, function(x) cub_survivorship(x[1], x[2], x[3]))

ggplot(dfr) + geom_line(aes(sm, sc, col = as.factor(nf))) + ggtitle('Survivorship')

gg <- ggplot(dfr) + geom_line(aes(nm, sc, col = as.factor(nf))) + ggtitle('Cub survivorship') + labs(x = 'Number males', y = '', col = 'Number\nfemales')

pdfr(gg, width = 10, name = 'cub_survivorship')

birth_rate <- function(nfemales, nmales, smales, harem_size_equilibrium = 2, clutch_size = 2) {
    
    pencounter <- probability_encounter(nmales, nfemales)
    
    pinfanticide <- pencounter * (1 - smales)
    
    ninfanticide <- nfemales * pinfanticide
    
    harem_size  <- (harem_size_equilibrium * nfemales + ninfanticide) / nfemales
    
    nharems <- nfemales / harem_size
    
    cubs <- harem_size * clutch_size * 2 * nmales * nharems / (nmales + nharems)
    
    return(list(h = harem_size, b = cubs))
}

dfr$hs <- apply(dfr, 1, function(x) birth_rate(x[1], x[2], x[3])[[1]])
dfr$br <- apply(dfr, 1, function(x) birth_rate(x[1], x[2], x[3])[[2]])

ggplot(dfr) + geom_point(aes(hs, br, col = nm)) + facet_wrap(~nf)

ggplot(dfr) + geom_line(aes(hs, sc, col = as.factor(nf))) + ggtitle('Harem size')

ggplot(dfr) + geom_line(aes(br, sc, col = as.factor(nf))) + ggtitle('Birth rate')

ggplot(dfr) + geom_line(aes(sm, hs, col = as.factor(nf))) + ggtitle('Harem size')

ggplot(dfr) + geom_line(aes(sm, br, col = nm), size = 2) + facet_wrap(~nf) + ggtitle('Birth rate')

ggplot(dfr) + geom_line(aes(nm, hs, col = as.factor(nf))) + ggtitle('Harem size')

ggplot(dfr) + geom_line(aes(nm, br, col = as.factor(nf))) + ggtitle('Birth rate')


gg <- ggplot(dfr) + geom_line(aes(sm, ss)) + facet_wrap(~age) + ggtitle('Survivorship') + labs(x = 'Survivorship of mature males', y = '')

pdfr(gg, width = 10, name = 'sinf')


