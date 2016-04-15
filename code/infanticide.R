

library(ggplot2)

source('utils/pdfr.r')

# numbers
# from Swanepoel et al. 2014
nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

# harem size
harem_size <- 5

# create orthogonal data.frame
nmales   <- nmales * seq(0, 2, length = 101)
nfemales <- seq(100, 600, by = 100)

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

gg <- ggplot(dfr) + geom_line(aes(nm, pe, col = as.factor(nf))) + ggtitle('Probability of male occupancy per harem\n') + labs(x = 'Number of mature males', y = '', col = 'Number of mature \nfemales')

pdfr(gg, width = 10, name = 'pocc')

cub_survivorship <- function(nm, nf, sm) {
    
    nmales   <- nm
    nfemales <- nf
    smales   <- sm
    
    scub_star <- 0.3270833
    sjuv_star <- 0.7197452
    
    pencounter <- probability_encounter(nmales, nfemales)
    
    scub <- scub_star * (1 - pencounter * (1 - smales))
    sjuv <- sjuv_star * (1 - pencounter * (1 - smales))
    
    return(list(scub = scub, sjuv = sjuv))
    
}

nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

# create orthogonal data.frame
smales   <- seq(0, 1, length = 101)

dfr <- data.frame(nm = nmales, nf = nfemales, sm = smales)
dfr$nm <- dfr$nm * dfr$sm

dfr <- rbind(data.frame(dfr, age = 'Cub',      ss = apply(dfr, 1, function(x) cub_survivorship(x[1], x[2], x[3])[[1]])),
             data.frame(dfr, age = 'Juvenile', ss = apply(dfr, 1, function(x) cub_survivorship(x[1], x[2], x[3])[[2]])))

ggplot(dfr) + geom_line(aes(sm, ss)) + facet_wrap(~age) + ggtitle('Survivorship')

gg <- ggplot(dfr) + geom_line(aes(sm, ss)) + facet_wrap(~age) + ggtitle('Survivorship') + labs(x = 'Survivorship of mature males', y = '')

pdfr(gg, width = 10, name = 'sinf')






