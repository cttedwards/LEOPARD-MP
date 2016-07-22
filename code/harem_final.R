


library(ggplot2)
library(plyr)

# numbers (mature)
nfemales <- sum(c(2,2,2,3,17))
nmales   <- sum(c(2,2,3,3,9))

mean.smales <- sum(c(0.9642857, 1.0000000, 1.0000000, 0.9000000, 0.2857143) * c(2,2,3,3,9)) / sum(c(2,2,3,3,9))


nfemales
nmales
hsize <- 1.14
nmales * hsize
nfemales / (nmales * hsize)
pfbreed <- (nmales * hsize) / nfemales
sratio <- nmales / nfemales
sprop  <- nmales / (nmales + nfemales)
# => at the current sex ratio only half the 
# => females are in harems assuming that all 
# => males breed

# the probability that a female breeds follows a binomial distribution
# with the number of trials equal to the number of females and the probability
# of success calculated so that nfemales * psuccess -> nmales * haremsize
# currently we have 846 * psuccess = 556 * 1.14, giving psuccess = 0.75

# we can write the relationship between psuccess and the sex proportion as a binomial
# regression and fit to the known value to give the psuccess and harem size

# probability of female breeding is
# 0.75 at current sex ratio and zero
# when there are no males
pfbreed.ff <- function(sratio) {
    
    ff <- function(beta, sr) 1 - exp(-beta * sr)
    
    obj3 <- function(par) {
        
        beta <- par
        
        (pfbreed - ff(beta, sr = sratio))^2 + (0 - ff(beta, sr = 0))^2 
        
    }
    
    hat <- optimize(obj3, interval = c(-100, 100))$minimum
    
    return(list(pfb = ff(hat, sratio), beta = hat))
}

# checks
pfbreed.ff(sratio)[['pfb']]          # == pfbreed = 0.749
(1/sratio) * pfbreed.ff(sratio)[['pfb']] # == hsize = 1.14

# functions
pfb <- function(sr, beta) 1 - exp(-beta * sr)
hsz <- function(sr, beta) (1 / sr) * pfb(sr, beta)

# hat values
betahat <- pfbreed.ff(sratio)[['beta']]

# plot
sr.seq <- seq(0, 2, length = 101)

windows(width = 14)
par(mfrow = c(1, 2))
plot(range(sr.seq), c(0, 1.1), type = 'n', main = 'Probability of female breeding', xlab = 'Sex ratio (M:F)', ylab = '')
lines(sr.seq, pfb(sr.seq, betahat), lty = 1)
abline(v = sratio, col = 2); abline(h = pfbreed, col = 2, lty = 2)

plot(range(sr.seq), c(0, 3), type = 'n', main = 'Harem size', xlab = 'Sex ratio (M:F)', ylab = '')
lines(sr.seq, hsz(sr.seq, betahat), lty = 1)
abline(v = sratio, col = 2); abline(h = hsize, col = 2)
savePlot(filename = '../report/harem_size.pdf', type = "pdf")
dev.off()

# probability of cub survivorship is
# 0.33 at current sex ratio and zero
# when there are no females

# encounter prob. per harem between males and harems

nmales   <- seq(0, 100, by = 1)
nfemales <- seq(0, 100, by = 1)

dfr <- expand.grid(nm = nmales, nf = nfemales)
dfr <- subset(dfr, nm > 0 & nf > 0)

dfr$sr   <- dfr$nm / dfr$nf

dfr <- dfr[order(dfr$sr),]

dfr <- subset(dfr, sr < 2)

dfr$sm <- mean.smales

# encounter rate per harem
ecr <- function(sm, nm, nf, beta) {
    
    nm <- nm * sm
    
    hs <- hsz(nm / nf, beta)
    
    nh <- nf / hs
    
    2 * nm / (nm + nh)
}


dfr$ec <- ecr(dfr$sm, dfr$nm, dfr$nf, betahat)

# p infanticide = p encounter * p male dieing
pin <- function(sm, nm, nf, beta) {
    
    # probability of encounter following
    # male mortality
    ec <- ecr(sm, nm, nf, beta)
    
    # probability of mortality
    # and encounter
    (1 - sm) * ec
}

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, betahat)

windows(width = 14)
par(mfrow = c(1, 2))
plot(ec~sr, dfr, type = 'l', main = "Encounter probability", xlab = 'Sex ratio (M:F)', ylab = '', ylim = c(0,1))
abline(v = sratio, col = 2); abline(h = ecr(mean.smales, 19, 26, betahat), col = 2, lty = 2)

plot(pi~sr, dfr, type = 'l', main = "Probability of infanticide", xlab = 'Sex ratio (M:F)', ylab = '', ylim = c(0,1))
abline(v = sratio, col = 2); abline(h = pin(mean.smales, 19, 26, betahat), col = 2, lty = 2)
savePlot(filename = '../report/prob_inf.pdf', type = "pdf")
dev.off()

# base level of cub survivorship
# cub can die either from natural causes or infanticide
# 0.33 = 1 - ((1 - csurv) + (1 - pinf) - intersect)

# csurv (if intersect = 0)
csurv <- 0.3270833 / (1 - pin(sm = mean.smales, nm = 19, nf = 26, beta = betahat))

# check
csv <- function(sm, nm, nf, beta) csurv * (1 - pin(sm, nm, nf, beta))

csv(sm = mean.smales, nm = 19, nf = 26, beta = betahat)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, betahat)

plot(cs~sr, dfr, type = 'l', main = "Cub survival", xlab = 'Sex ratio (M:F)', ylab = '')
abline(v = sratio, col = 2); abline(h = csv(sm = mean.smales, nm = 19, nf = 26, beta = betahat), col = 2)

# explore impact of male survivorship
nmales   <- seq(0, 100, by = 1)
nfemales <- seq(0, 100, by = 1)

dfr <- expand.grid(nm = nmales, nf = nfemales, sm = seq(0.1, 0.9, length = 5))
dfr <- subset(dfr, nm > 0 & nf > 0)

dfr$nm <- dfr$nm * dfr$sm

dfr$sr   <- dfr$nm / dfr$nf

dfr <- dfr[order(dfr$sr),]

dfr <- subset(dfr, sr < 2)

dfr$ec <- ecr(dfr$sm, dfr$nm, dfr$nf, betahat)

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, betahat)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, betahat)

windows(width = 14)
#ggplot(dfr) + geom_line(aes(sr, pi, col = as.factor(sm)), size = 2) + 
#    theme_bw() + 
#    labs(x = 'Sex ratio (M:F)', y = '', col = 'Male\nsurvivorship') +
#    ggtitle('Probability of infanticide')
ggplot(dfr) + geom_line(aes(sr, cs, col = as.factor(sm)), size = 2) + 
    theme_bw() + 
    labs(x = 'Sex ratio (M:F)', y = '', col = 'Male\nsurvivorship') +
    ggtitle('Cub survivorship')
savePlot(filename = '../report/cub_surv.pdf', type = "pdf")
dev.off()


# re-explore impact of male survivorship
dfr <- expand.grid(nm = seq(0, 100, by = 1), sr = seq(0.1, 1.2, length = 3), sm = seq(0.1, 0.9, length = 101))

dfr$nm <- dfr$nm * dfr$sm

dfr$nf   <- dfr$nm / dfr$sr

dfr <- subset(dfr, nm > 0 & nf < 100)

dfr <- dfr[order(dfr$sm),]

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, betahat)

ggplot(dfr) + geom_line(aes(sm, pi)) + facet_wrap(~sr)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, betahat)

ggplot(dfr) + geom_line(aes(sm, cs)) + facet_wrap(~sr)

# re-explore impact of male survivorship
# at current sex ratio
#dfr <- expand.grid(nm = seq(1, 100, by = 1), sr = sratio, sm = seq(0.1, 1.0, length = 100))
dfr <- expand.grid(nm = 19, sr = sratio, sm = seq(0.1, 1.0, length = 100))

# number of males adjusted for survivorship
dfr$nm <- dfr$nm * dfr$sm

# number of females at current sratio
dfr$nf   <- dfr$nm / dfr$sr

# adjust sex ratio
dfr$sr   <- dfr$nm / dfr$nf

# order by sm
dfr <- dfr[order(dfr$sm),]

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, betahat)

#ggplot(dfr) + geom_line(aes(sm, pi))

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, betahat)

ggplot(dfr) + geom_line(aes(sm, cs)) + geom_hline(yintercept = 0.3270833) + geom_vline(xintercept = mean.smales)

# re-explore impact of male survivorship
# at adjusted sex ratio
dfr <- expand.grid(nm = 19, sr = sratio, sm = seq(0, 1.0, length = 100))

# number of females at current sratio
dfr$nf   <- dfr$nm / dfr$sr

# adjust sex ratio
dfr$sr   <- dfr$nm * dfr$sm / dfr$nf

# order by sm
dfr <- dfr[order(dfr$sm),]

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, betahat)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, betahat)

windows(width = 14)
ggplot(dfr) + 
    geom_line(aes(sm, cs, col = sr), size = 2) + 
    #geom_hline(yintercept = 0.3270833, col = "red") + geom_vline(xintercept = mean.smales, col = "red") +
    geom_hline(yintercept = csv(mean.smales, 19, 26, betahat), col = "red") + geom_vline(xintercept = mean.smales, col = "red") +
    theme_bw() + 
    labs(x = 'Male survivorhsip', y = '', col = 'Sex\nratio') +
    ggtitle('Cub survivorship')
savePlot(filename = '../report/cub_surv2.pdf', type = "pdf")
dev.off()









