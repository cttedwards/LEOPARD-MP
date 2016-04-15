

library(ggplot2)

# numbers
# from Swanepoel et al. 2014
nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

# harem size
harem_size <- 5

# create orthogonal data.frame
nseq <- seq(0, 2, length = 6)

nmales   <- nmales * nseq
nfemales <- nfemales * nseq

dfr <- expand.grid(nm = nmales, nf = nfemales)
dfr$nh <- dfr$nf / harem_size

number_encounters <- function(nm, nf) {
    
    nmales   <- nm
    nfemales <- nf
    
    # number of harems
    nharems     <- nfemales / harem_size
    
    # number of encounters between males and harems
    nencounters <- 2 * nmales * nharems / (nmales + nharems)
    
    # maximum number of encounters per harem assuming
    # females mate 3 times per year
    nencounters_max <- nharems * 3
    
    # NB this is obsolete because nencounters approaches nharems * 2 as nmales tends to infinity
    
    # truncated at <= nencounters_max
    nencounters <- vapply(nencounters, function(x) min(x, nencounters_max), numeric(1))
    
    # return
    return(nencounters)
}

dfr$ne <- apply(dfr, 1, function(x) number_encounters(x[1], x[2]))

ggplot(dfr) + geom_line(aes(nm, ne, col = as.factor(nf))) + ggtitle('number_encounters()')

probability_encounter <- function(nm, nf) {
    
    nmales   <- nm
    nfemales <- nf
    
    nencounters <- number_encounters(nmales, nfemales)
    
    nharems     <- nfemales / harem_size
    
    nencounters_perharem <- nencounters / nharems
    
    pencounter <- 1 - exp(-nencounters_perharem)
    
    return(pencounter)
}

dfr$pe <- apply(dfr, 1, function(x) probability_encounter(x[1], x[2]))

ggplot(dfr) + geom_line(aes(nm, pe, col = as.factor(nf))) + ggtitle('probability_encounter()')


# probability of encounter should have an origin at zero and asymptote at 1 as nencounters approaces its
# asymptote at nharems * 2

logistic <- function(x) 1 / (1 + exp(-x))

dfr$value2 <- logistic(dfr$value)


# number of harems over the number of encounters
probability_encounter <- function(nm, nf) {
    
    nmales   <- nm
    nfemales <- nf
    
    # number of harems
    nharems     <- nfemales / harem_size
    
    # number of encounters (truncated)
    nencounters <- number_encounters(nmales, nfemales)
    
    # probability of encounter
    pencounter    <- nharems / nencounters
    
    # truncated at <= 1
    pencounter    <- vapply(pencounter, function(x) min(x, 1), numeric(1))
    
    # return
    return(pencounter)
}

dfr$value <- apply(dfr, 1, function(x) probability_encounter(x[1], x[2]))

ggplot(dfr) + geom_line(aes(nm, value, col = as.factor(nf))) + ggtitle('probability_encounter()')

ggplot(dfr) + geom_line(aes(nf, value, col = as.factor(nm))) + ggtitle('probability_encounter()')


# probability of infanticide
probability_infanticide <- function() {
    
    
    
    
    
}


# probability of occupancy
# i.e. probability of a male being present
# (number of males over the number of male/female
# encounters)
pocc <- function(nm, nf, h = 5) {
    nh <- nf / h
    nh / nocc(nm, nf)
}

# probability of infanticide
# (probaility of male being killed multiplied by 
# probability of occupancy)
pinf <- function(sm, nm, nf) {
    sm * pocc(nm, nf)
}

xx <- pocc(nm * nseq, nf)

plot(nm * nseq, xx)



