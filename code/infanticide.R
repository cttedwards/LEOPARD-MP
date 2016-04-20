

library(ggplot2)

source('utils/pdfr.r')

# numbers
# from Swanepoel et al. 2014
nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

# create orthogonal data.frame
nmales   <- nmales * seq(0, 2, length = 101)
nfemales <- seq(100, 600, by = 100)

dfr <- expand.grid(nm = nmales, nf = nfemales)

# number of encounters between males and females
# follows the harmonic mean function
number_encounters <- function(nmales, nfemales) {
    
    # number of encounters between males and females
    nencounters <- 2 * nmales * nfemales / (nmales + nfemales)
    
    # return
    return(nencounters)
}

dfr$ne <- apply(dfr, 1, function(x) number_encounters(x[1], x[2]))

ggplot(dfr) + geom_line(aes(nm, ne, col = as.factor(nf))) + ggtitle('number_encounters()')

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

dfr$sm <- seq(0, 1, length = 101)

#dfr <- data.frame(nm = nmales, nf = nfemales, sm = smales)
#dfr$nm <- dfr$nm * dfr$sm

dfr$sc <- apply(dfr, 1, function(x) cub_survivorship(x[1], x[2], x[3]))

ggplot(dfr) + geom_line(aes(sm, sc)) + ggtitle('Survivorship')

gg <- ggplot(dfr) + geom_line(aes(sm, ss)) + facet_wrap(~age) + ggtitle('Survivorship') + labs(x = 'Survivorship of mature males', y = '')

pdfr(gg, width = 10, name = 'sinf')


birth_rate <- function(nf, nm, sm) {
    
    nmales   <- nm
    nfemales <- nf
    smales   <- sm
    
    pencounter <- probability_encounter(nmales, nfemales)
    
    pinfanticide <- pencounter * (1 - smales)
    
    ni <- nf * pinfanticide
    
    harem_size  <- (3 * nf + ni) / nf
    clutch_size <- 2
    
    nh <- nf / harem_size
    
    cubs <- harem_size * clutch_size * 2 * nm * nh / (nm + nh)
    
    return(list(h = harem_size, b = cubs))
}

nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

# create orthogonal data.frame
smales   <- seq(0, 1, length = 101)

dfr <- data.frame(nm = nmales, nf = nfemales, sm = smales)

dfr$sc <- apply(dfr, 1, function(x) cub_survivorship(x[1], x[2], x[3]))
dfr$hs <- apply(dfr, 1, function(x) birth_rate(x[1], x[2], x[3])[[1]])
dfr$br <- apply(dfr, 1, function(x) birth_rate(x[1], x[2], x[3])[[2]])

ggplot(dfr) + geom_line(aes(sc, hs)) + ggtitle('Harem size')

ggplot(dfr) + geom_line(aes(sc, br)) + ggtitle('Birth rate')

ggplot(dfr) + geom_line(aes(sm, hs)) + ggtitle('Harem size')

ggplot(dfr) + geom_line(aes(sm, br)) + ggtitle('Birth rate')


gg <- ggplot(dfr) + geom_line(aes(sm, ss)) + facet_wrap(~age) + ggtitle('Survivorship') + labs(x = 'Survivorship of mature males', y = '')

pdfr(gg, width = 10, name = 'sinf')


