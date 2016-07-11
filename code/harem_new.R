


library(ggplot2)
library(plyr)

source('utils/pdfr.r')

# numbers (mature)
nfemales <- sum(c(92,67,63,62,93,469))
nmales   <- sum(c(39,62,47,78,84,246))


nfemales
nmales
hsize <- 1.14
nmales * hsize
nfemales / (nmales * hsize)
pfbreed <- (nmales * hsize) / nfemales
sratio <- nmales / nfemales
#sprop  <- nmales / (nmales + nfemales)
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
pfbreed.ff <- function(sr) {
    
    #ff  <- function(mu, beta) (exp(-mu - beta * sr)) / (exp(-mu - beta * sr) + 1)
    
    ff <- function(mu, beta) 1 - exp(-mu - beta * sr)
    
    #obj <- function(pars) {
    #    
    #    mu   <- pars[1]
    #    beta <- pars[2]
    #    
    #    (pfbreed - ff(mu, beta))^2 + (0 - (exp(-mu - beta * 0)) / (exp(-mu - beta * 0) + 1))^2 + ((0 - mu)^2) / 1e4
    #    
    #}
    #
    #hat <- optim(c(1, -1), obj, control = list(reltol = 1e-16))$par
    
    #obj2 <- function(pars) {
    #    
    #    mu   <- 0
    #    beta <- pars[1]
    #    
    #    (pfbreed - ff(mu, beta))^2 + (0 - (exp(-mu - beta * 0)) / (exp(-mu - beta * 0) + 1))^2 
    #    
    #}
    #
    #hat <- c(0, optimize(obj2, interval = c(-100, 100))$minimum)
    
    obj3 <- function(pars) {
        
        mu   <- 0
        beta <- pars[1]
        
        (pfbreed - ff(mu, beta))^2 + (0 - (1 - exp(-mu - beta * 0)))^2 
        
    }
    
    hat <- c(0, optimize(obj3, interval = c(-100, 100))$minimum)
    
    return(list(pfb = ff(hat[1], hat[2]), mu = hat[1], beta = hat[2]))
}

#pfbreed.ff <- function(sratio) 0.7492199

# checks
pfbreed.ff(sratio)[['pfb']]          # == pfbreed = 0.749
(1/sratio) * pfbreed.ff(sratio)[['pfb']] # == hsize = 1.14

# functions
#pfb <- function(sr, mu, beta) rep(0.7492199, length(sr))
#pfb <- function(sr, mu, beta) (exp(-mu - beta * sr)) / (exp(-mu - beta * sr) + 1)
pfb <- function(sr, mu, beta) 1 - exp(-mu - beta * sr)
hsz <- function(sr, mu, beta) (1 / sr) * pfb(sr, mu, beta)

hsz.se <- function(nm, nf, mu, beta) {
    p <- pfb(sr, mu, beta)
    sqrt(nf * p * (1 - p) / (nm ^ 2))
}

# encounter rate per harem
ecr <- function(nm, nf, mu, beta) {
    
    hs <- hsz(nm / nf, mu, beta)
    
    nh <- nf / hs
    
    2 * nh * nm / (nm + nh) * (1 / nh)
}


# hat values
muhat   <- pfbreed.ff(sratio)[['mu']]
betahat <- pfbreed.ff(sratio)[['beta']]

# plot
sr.seq <- seq(0, 2, length = 101)

plot(range(sr.seq), c(0, 1.1), type = 'n', main = 'Probability of female breeding', xlab = 'Sex ratio (M:F)', ylab = '')
lines(sr.seq, pfb(sr.seq, muhat, betahat), lty = 1)
abline(v = sratio, col = 2); abline(h = pfbreed, col = 2)

plot(range(sr.seq), c(0, 3), type = 'n', main = 'Harem size', xlab = 'Sex ratio (M:F)', ylab = '')
lines(sr.seq, hsz(sr.seq, muhat, betahat), lty = 1)
abline(v = sratio, col = 2); abline(h = hsize, col = 2)


# probability of cub survivorship is
# 0.33 at current sex ratio and zero
# when there are no females

# encounter prob. per harem between males and harems

nmales   <- seq(0, 1000, by = 10)
nfemales <- seq(0, 1000, by = 10)

dfr <- expand.grid(nm = nmales, nf = nfemales)
dfr <- subset(dfr, nm > 0 & nf > 0)

dfr$sr   <- dfr$nm / dfr$nf

dfr <- dfr[order(dfr$sr),]

dfr <- subset(dfr, sr < 2)

dfr$ec <- ecr(dfr$nm, dfr$nf, muhat, betahat)

plot(ec~sr, dfr, type = 'l', main = "Encounter probability", xlab = 'Sex ratio (M:F)', ylab = '')
abline(v = sratio, col = 2); abline(h = 1, col = 2, lty = 2)

# p infanticide = p encounter * p male dieing
pin <- function(sm, nm, nf, mu, beta) {
    
    # probability of encounter following
    # male mortality
    ec <- ecr(sm * nm, nf, mu, beta)
    
    # probability of mortality
    # and encounter
    (1 - sm) * ec
}

# base level of cub survivorship
# cub can die either from natural causes or infanticide
# 0.33 = 1 - ((1 - csurv) + (1 - pinf) - intersect)

msurv <- 0.9
nfemales  <- sum(c(92,67,63,62,93,469))
nmales    <- sum(c(39,62,47,78,84,246))

# csurv (if intersect = 0)
csurv <- 0.3270833 / (1 - pin(msurv, nmales, nfemales, muhat, betahat))

# check
csv <- function(sm, nm, nf, mu, beta) csurv * (1 - pin(sm, nm, nf, mu, beta))

csv(msurv, nmales, nfemales, muhat, betahat)

dfr$sm <- 0.9

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

plot(pi~sr, dfr, type = 'l', main = "Probability of infanticide", xlab = 'Sex ratio (M:F)', ylab = '')
abline(v = sratio, col = 2); abline(h = pin(msurv, nmales, nfemales, muhat, betahat), col = 2, lty = 2)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

plot(cs~sr, dfr, type = 'l', main = "Cub survival", xlab = 'Sex ratio (M:F)', ylab = '')
abline(v = sratio, col = 2); abline(h = csv(msurv, nmales, nfemales, muhat, betahat), col = 2, lty = 2)

# explore impact of male survivorship
nmales   <- seq(0, 1000, by = 10)
nfemales <- seq(0, 1000, by = 10)

dfr <- expand.grid(nm = nmales, nf = nfemales, sm = seq(0.1, 0.9, length = 11))
dfr <- subset(dfr, nm > 0 & nf > 0)

dfr$nm <- dfr$nm * dfr$sm

dfr$sr   <- dfr$nm / dfr$nf

dfr <- dfr[order(dfr$sr),]

dfr <- subset(dfr, sr < 2)

dfr$ec <- ecr(dfr$nm, dfr$nf, muhat, betahat)

ggplot(dfr) + geom_line(aes(sr, ec, col = as.factor(sm)))

plot(ec~sr, dfr, type = 'l', main = "Encounter probability", xlab = 'Sex ratio (M:F)', ylab = '')
abline(v = sratio, col = 2); abline(h = 1, col = 2, lty = 2)

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

ggplot(dfr) + geom_line(aes(sr, pi, col = as.factor(sm)))

boxplot(pi~sm, dfr, main = "Probability of infanticide", xlab = 'Male survivorship', ylab = '', outline = FALSE)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

boxplot(cs~sm, dfr, main = "Cub survivorship", xlab = 'Male survivorship', ylab = '', outline = FALSE)


# re-explore impact of male survivorship
dfr <- expand.grid(nm = seq(0, 1000, by = 10), sr = seq(0.1, 1.2, length = 11), sm = seq(0.1, 0.9, length = 101))

dfr$nm <- dfr$nm * dfr$sm

dfr$nf   <- dfr$nm / dfr$sr

dfr <- subset(dfr, nm > 0 & nf < 1000)

dfr <- dfr[order(dfr$sm),]

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

ggplot(dfr) + geom_line(aes(sm, pi)) + facet_wrap(~sr)

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

ggplot(dfr) + geom_line(aes(sm, cs)) + facet_wrap(~sr)

# re-explore impact of male survivorship
# at current sex ratio
dfr <- expand.grid(nm = seq(100, 1000, by = 1), sr = sratio, sm = seq(0.1, 1.0, length = 100))

# number of males adjusted for survivorship
dfr$nm <- dfr$nm * dfr$sm

# number of females at current sratio
dfr$nf   <- dfr$nm / dfr$sr

# adjust sex ratio
dfr$sr   <- dfr$nm / dfr$nf

# order by sm
dfr <- dfr[order(dfr$sm),]

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

#ggplot(dfr) + geom_line(aes(sm, pi))

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

ggplot(dfr) + geom_line(aes(sm, cs)) + geom_hline(yintercept = 0.3270833) + geom_vline(xintercept = 0.9)

# re-explore impact of male survivorship
# at adjusted sex ratio
dfr <- expand.grid(nm = 600, sr = sratio, sm = seq(0.1, 1.0, length = 100))

# number of females at current sratio
dfr$nf   <- dfr$nm / dfr$sr

# number of males adjusted for survivorship
dfr$nm <- dfr$nm * dfr$sm

# adjust sex ratio
dfr$sr   <- dfr$nm / dfr$nf

# order by sm
dfr <- dfr[order(dfr$sm),]

dfr$pi <- pin(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

#ggplot(dfr) + geom_line(aes(sm, pi))

dfr$cs <- csv(dfr$sm, dfr$nm, dfr$nf, muhat, betahat)

ggplot(dfr) + geom_line(aes(sm, cs)) + geom_hline(yintercept = 0.3270833) + geom_vline(xintercept = 0.9)

csv(0.9, nmales, nfemales, muhat, betahat)
nmales / nfemales



















# expression for harem size
#rho <- sprop
#gam <- nfemales / nmales
#x <- expression(sr * (exp(-mu - beta * sp)) / (exp(-mu - beta * sp) + 1))
#eval(x)


birth <- function(nm, nf, hs) 2 * nm / (nm + (nf/hs))
encounter <- function(nm, nf, hs) 2 * nf * 1 / (nm + (nf / hs))
pfb  <- function(sr, mu, beta) (exp(-mu - beta * sr)) / (exp(-mu - beta * sr) + 1)

# create orthogonal data.frame
nmales   <- seq(0, 1000, by = 10)
nfemales <- seq(0, 1000, by = 10)

dfr <- expand.grid(nm = nmales, nf = nfemales)
dfr <- subset(dfr, nm > 0 & nf > 0)

dfr$sr   <- dfr$nm / dfr$nf

dfr <- subset(dfr, sr < 2)

dfr$mu   <- pfbreed.ff(sratio)[['mu']]
dfr$beta <- pfbreed.ff(sratio)[['beta']]

dfr$fb <- apply(dfr, 1, function(x) pfb(x['sr'], x['mu'], x['beta']))
dfr$hs <- apply(dfr, 1, function(x) (1 / x['sr']) * pfb(x['sr'], x['mu'], x['beta']))
    
dfr$br <- apply(dfr, 1, function(x) birth(x['nm'], x['nf'], x['hs']))
    
dfr$ne <- apply(dfr, 1, function(x) encounter(x['nm'], x['nf'], x['hs']))
    
dfr$ii <- apply(dfr, 1, function(x) {
        
        nf <- x['nf']
        nm <- x['nm']
        hs <- x['hs']
        sr <- x['sr']
        
        eval(D(expression(2 * sr * nm / ((nf / sr) + (sr * nm / hs))), "sr"))
    
    })
    

ggplot(dfr) + geom_line(aes(sr, fb)) + labs(x = 'Sex ratio (M:F)', y = 'Prob. female breeding') + geom_hline(yintercept = pfbreed) + geom_vline(xintercept = sratio)

ggplot(res) + geom_line(aes(sr, hs)) + labs(x = 'Sex ratio (M:F)', y = 'Harem size') + geom_hline(yintercept = hsize) + geom_vline(xintercept = sratio)

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

