


library(ggplot2)
library(plyr)

source('utils/pdfr.r')

# numbers (mature)
# from Swanepoel et al. 2014
nfemales <- sum(c(92,67,63,62,93,469))
nmales   <- sum(c(39,62,47,78,84,246))


nfemales
nmales
hsize <- 1.14
nmales * hsize
nfemales / (nmales * hsize)
pfbreed <- (nmales * hsize) / nfemales
sratio <- nfemales / nmales
sprop  <- nfemales / (nmales + nfemales)
# => at the current sex ratio only half the 
# => females are in harems assuming that all 
# => males breed

# the probability that a female breeds follows a binomial distribution
# with the number of trials equal to the number of females and the probability
# of success calculated so that nfemales * psuccess -> nmales * haremsize
# currently we have 846 * psuccess = 556 * 1.14, giving psuccess = 0.75

# we can write the relationship between psuccess and the sex proportion as a binomial
# regression and fit to the known value to give the psuccess and harem size

pfbreed.ff <- function(sp, beta = 10) {
    
    ff  <- function(mu) (exp(-mu - beta * sp)) / (exp(-mu - beta * sp) + 1)
    obj <- function(mu) pfbreed - ff(mu)
    
    muhat <- uniroot(obj, interval = c(-100, 100))$root
    
    return(list(pfb = ff(muhat), muhat = muhat))
}

# checks
pfbreed.ff(sprop)[['pfb']]          # == pfbreed = 0.749
sratio * pfbreed.ff(sprop)[['pfb']] # == hsize = 1.14

# plot
beta.seq <- seq(0, 10, by = 1)
muhat.seq <- numeric(length(beta.seq))
for (i in 1:length(muhat.seq)) muhat.seq[i] <- pfbreed.ff(sprop, beta.seq[i])[['muhat']]

sp.seq <- seq(0, 1, length = 101)
pfb  <- function(mu, sp, beta) (exp(-mu - beta * sp)) / (exp(-mu - beta * sp) + 1)
plot(sp.seq, sp.seq, type = 'n', ylab = 'Probability of female breeding', xlab = 'Proportion females')
for (i in 1:length(muhat.seq))
    lines(sp.seq, pfb(muhat.seq[i], sp.seq, beta.seq[i]), lty = 2)
abline(v = sprop, col = 2); abline(h = pfbreed, col = 2)

# expression for harem size
#rho <- sprop
#gam <- nfemales / nmales
#x <- expression(gam * (exp(-muhat - 10 * sprop)) / (exp(-muhat - 10 * rho) + 1))
#eval(x)

# create orthogonal data.frame
nmales   <- seq(0, 1000, by = 10)
nfemales <- seq(0, 1000, by = 10)

dfr <- expand.grid(nm = nmales, nf = nfemales)
dfr <- subset(dfr, nm > 0 | nf > 0)
dfr$sp   <- dfr$nf / (dfr$nm + dfr$nf)
dfr$sr   <- dfr$nf / dfr$nm

birth <- function(nm, nf, hs) 2 * nm / (nm + (nf/hs))

encounter <- function(nm, nf, hs) 2 * nf * 1 / (nm + (nf / hs))

res <- data.frame(nm = numeric(), nf = numeric(), sp = numeric(), sr = numeric(), mu = numeric(), beta = numeric(), fb = numeric(), hs = numeric())

for (i in 1:length(muhat.seq)) {
    
    dfr$mu   <- muhat.seq[i]
    dfr$beta <- beta.seq[i]
    
    dfr$fb <- apply(dfr, 1, function(x) pfb(x['mu'], x['sp'], x['beta']))
    dfr$hs <- apply(dfr, 1, function(x) min(3, x['sr'] * pfb(x['mu'], x['sp'], x['beta'])))
    
    dfr$br <- apply(dfr, 1, function(x) birth(x['nm'], x['nf'], x['hs']))
    
    dfr$ne <- apply(dfr, 1, function(x) encounter(x['nm'], x['nf'], x['hs']))
    
    dfr$ii <- apply(dfr, 1, function(x) {
        
        nf <- x['nf']
        nm <- x['nm']
        hs <- x['hs']
        sr <- x['sr']
        
        eval(D(expression(2 * sr * nm / ((nf / sr) + (sr * nm / hs))), "sr"))
    
    })
    
    res <- rbind(res, dfr)
    
}

ggplot(res) + geom_line(aes(sp, hs, col = sr)) + labs(x = 'Proportion females', y = 'Harem size') + geom_hline(yintercept = hsize) + geom_vline(xintercept = sprop)

ggplot(res) + geom_line(aes(sp, hs, col = as.factor(beta))) + labs(x = 'Proportion females', y = 'Harem size') + geom_hline(yintercept = hsize) + geom_vline(xintercept = sprop)

ggplot(res) + geom_line(aes(sp, br, col = as.factor(beta))) + labs(x = 'Proportion females', y = 'Litters per female') + geom_vline(xintercept = sprop)

ggplot(res) + geom_line(aes(sp, ne, col = as.factor(beta))) + labs(x = 'Proportion females', y = 'Number female encounters per male') + geom_vline(xintercept = sprop)

ggplot(res) + geom_line(aes(sp, 1 - ii, col = as.factor(beta))) + labs(x = 'Proportion females', y = 'Rate of change in number female encounters per male') + geom_vline(xintercept = sprop)

# beta should be >1 because proportion of females breeding should approach zero as
# proportion of females approaches one. However we still expect the harem size to 
# increase as the proportion of females increases.
# Therefore choose maximum beta so that harem size increases monotonically with proportion females

# probability of infanticide is:
# new_nmales * new_hs / (old_nmales * old_hs) = new_nfemales_occupied / old_nfemales_occupied = 1 + prob_infanticide if > 1

#### END #####

# as the sex ratio changes we expect psuccess to 
# also change and we can use this to back calculate
# the harem size at different sex ratio values

plot(pfb ~ sp, dfr, main = 'Prob. that a female breeds \n(i.e. is in a harem)', ylab = '', xlab = 'sex ratio (N females / N males)')
abline(v = sprop); abline(h = pfbreed)

plot(hs ~ sp, dfr, main = 'Harem size', ylab = '', xlab = 'sex ratio (N females / N males)')
abline(v = sprop); abline(h = hsize)

# create orthogonal data.frame
nmales   <- seq(0, 1000, by = 100)#nmales * seq(0, 2, length = 11)
nfemales <- seq(0, 1000, by = 100)

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

