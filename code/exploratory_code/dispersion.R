# explore dispersal of males

dfr <- expand.grid(nm = 500, sr = sratio, sm = seq(0.1, 2.0, length = 100))

# number of females at current sratio
dfr$nf   <- dfr$nm / dfr$sr

# adjust number of males
dfr$nm <- dfr$nm * dfr$sm

# recalculate sex ratio
dfr$sr   <- dfr$nm / dfr$nf

# remove survivorship
dfr$sm <- NULL

#nmales   <- seq(0, 1000, by = 10)
#nfemales <- seq(0, 1000, by = 10)
#
#dfr <- expand.grid(nm = nmales, sr = seq(0.1, 1.2, length = 3))
#
#dfr$nf <- dfr$nm / dfr$sr
#
#dfr <- subset(dfr, nm > 0 & nf < 1000)

dfr <- dfr[order(dfr$nm),]

#dfr <- subset(dfr, sr < 2)

pfb <- function(sr, beta) 1 - exp(-beta * sr)
hsz <- function(sr, beta) (1 / sr) * pfb(sr, beta)

# number of births per male
bpm <- function(nm, nf, beta, k = 2) {
    
    hs <- hsz(nm / nf, beta)
    
    nh <- nf / hs
    
    nb <- k * 2 * nm * nf / (nm + nh)
    
    nb / nm
}

dfr$bm <- bpm(dfr$nm, dfr$nf, betahat)

ggplot(dfr) + 
    geom_line(aes(nm, bm, col = sr), size = 2) + 
    geom_hline(yintercept = bpm(nm = 556, nf = 846, beta = betahat), col = "red") + geom_vline(xintercept = 556, col = "red") +
    theme_bw() + 
    labs(x = 'Number of adult males', y = '', col = 'Sex\nratio') +
    ggtitle('Births per male')


# number of dispersing adult males
ndp <- function(nm, nf, beta, bcrit = 4) {
    
    bm <- bpm(nm, nf, beta)
    
    nb <- bm * nm
    
    ff <- function(x, y, z) {a <- 0; if (y <= bcrit) a <- x - z / bcrit; a}
    
    if (length(bm > 1)) { apply(data.frame(nm, bm, nb), 1, function(x) ff(x['nm'], x['bm'], x['nb']))
    } else ff(nm, bm, nb)
}

dfr$nd <- ndp(dfr$nm, dfr$nf, betahat)

ggplot(dfr) + 
    geom_line(aes(nm, nd, col = sr), size = 2) + 
    geom_hline(yintercept = ndp(nm = 556, nf = 846, beta = betahat), col = "red") + geom_vline(xintercept = 556, col = "red") +
    theme_bw() + 
    labs(x = 'Number of adult males', y = '', col = 'Sex\nratio') +
    ggtitle('Number of dispersers')

dfr$bm <- bpm(dfr$nm - dfr$nd, dfr$nf, betahat)

ggplot(dfr) + 
    geom_line(aes(nm, bm, col = sr), size = 2) + 
    geom_hline(yintercept = bpm(nm = 556, nf = 846, beta = betahat), col = "red") + geom_vline(xintercept = 556, col = "red") +
    theme_bw() + 
    labs(x = 'Number of adult males', y = '', col = 'Sex\nratio') +
    ggtitle('Births per male (adjusted)')

dfr$dm <- dfr$bm / dfr$nm

ggplot(dfr) + 
    geom_line(aes(bm, dm, col = sr), size = 2)




