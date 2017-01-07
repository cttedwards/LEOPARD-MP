## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(leopard)
library(ggplot2)
library(reshape2)

## ------------------------------------------------------------------------
# survival rates per demographic category
ss <- c(0.3270833, 0.7197452, 0.9615385, 0.8882435, 0.9729730, 0.9382353, 0.9230769, 0.7219583, 0.9124424, 0.9642857, 1.0000000, 1.0000000, 0.9000000, 0.2857143)

# initial population numbers
nn <- c(nc  = 14,
        nj  = 7,
        saf = 3,
        f36 = 2,
        f48 = 2,
        f60 = 2,
        f72 = 3,
        f84 = 17,
        sam = 1,
        m36 = 2,
        m48 = 2,
        m60 = 3,
        m72 = 3,
        m84 = 9)

## ------------------------------------------------------------------------
# numbers (mature)
nfemales <- sum(nn[4:8])
nmales   <- sum(nn[10:14])
# current known harem size and sex ratio
hsize <- 1.14
sratio <- nmales / nfemales
# proportion of females breeding
pfbreed <- sratio * hsize

## ------------------------------------------------------------------------
xx <- leopard(nn, survival.rates = ss, litter.sizes = 2, deterministic = TRUE, harem.size = hsize)

## ------------------------------------------------------------------------
harem(xx)@harem.size

## ---- fig.width = 7, echo = FALSE----------------------------------------
pfb <- function(nf, nm) harem(c(nf, nm), beta = xx@harem.size.par) * (nm / nf)
hsz <- function(nf, nm) harem(c(nf, nm), beta = xx@harem.size.par)

nf <- 10
nm.seq <- seq(1, 20, length = 101)
sr.seq <- nm.seq/ nf

pfb.seq <- vapply(nm.seq, function(x) pfb(nf, x), numeric(1))
hsz.seq <- vapply(nm.seq, function(x) hsz(nf, x), numeric(1))

par(mfrow = c(1, 2))
plot(range(sr.seq), c(0, 1.1), type = 'n', main = 'Probability of female breeding', xlab = 'Sex ratio (M:F)', ylab = '')
lines(sr.seq, pfb.seq, lty = 1)
abline(v = sratio, col = 2); abline(h = pfbreed, col = 2, lty = 2)

plot(range(sr.seq), c(0, 3), type = 'n', main = 'Harem size', xlab = 'Sex ratio (M:F)', ylab = '')
lines(sr.seq, hsz.seq, lty = 1)
abline(v = sratio, col = 2); abline(h = hsize, col = 2)

## ------------------------------------------------------------------------
# encounter rate per harem
ecr <- function(nm, nf) {
    
    hs <- hsz(nf, nm)
    
    nh <- nf / hs
    
    harmonic_mean(nm, nh) * (1 / nh)
}

## ------------------------------------------------------------------------
# p infanticide = p encounter * p male dieing
pin <- function(sm, nm, sf, nf) {
    
    # probability of encounter following
    # male mortality
    ec <- ecr(sm * nm, sf * nf)
    
    # probability of mortality
    # and encounter
    (1 - sm) * ec
}

## ---- echo = FALSE-------------------------------------------------------
# set up data frame with a range of 
# realistic sex ratios
dfr    <- expand.grid(nm = seq(0, 100, by = 1), nf = seq(0, 100, by = 1))
dfr    <- subset(dfr, nm > 0 & nf > 0)
dfr$sr <- dfr$nm / dfr$nf
dfr    <- dfr[order(dfr$sr),]
dfr    <- subset(dfr, sr < 2)

# average survival
smales   <- sum(ss[10:14] * nn[10:14]) / sum(nn[10:14])
dfr$sm   <- smales
sfemales <- sum(ss[4:8] * nn[4:8]) / sum(nn[4:8])
dfr$sf   <- sfemales

# encounter rate
dfr$ec <- apply(dfr, 1, function(x) ecr(x['sm'] * x['nm'], x['sf'] * x['nf']))

# infanticide probability
dfr$pi <- apply(dfr, 1, function(x) pin(x['sm'], x['nm'], x['sf'], x['nf']))

## ---- echo = FALSE, fig.width = 7----------------------------------------
par(mfrow = c(1, 2))
plot(ec~sr, dfr, type = 'l', main = "Encounter probability", xlab = 'Sex ratio (M:F)', ylab = '', ylim = c(0,1))
abline(v = sratio, col = 2); abline(h = ecr(smales * 19, sfemales * 26), col = 2, lty = 2)

plot(pi~sr, dfr, type = 'l', main = "Probability of infanticide", xlab = 'Sex ratio (M:F)', ylab = '', ylim = c(0,1))
abline(v = sratio, col = 2); abline(h = pin(smales, 19, sfemales, 26), col = 2, lty = 2)

## ------------------------------------------------------------------------
infanticide(survival(xx))@realised.survival.rate[1]

## ---- echo = FALSE, fig.width=6------------------------------------------

simulated.survival <- numeric(1e1)
for (i in 1:length(simulated.survival)) {
    xx.tmp <- leopard(nn, survival.rates = ss, litter.sizes = 2, deterministic = FALSE, harem.size = hsize)
    simulated.survival[i] <- infanticide(survival(xx.tmp))@realised.survival.rate[1]
}

hist(simulated.survival, prob = TRUE)
abline(v = ss[1], col = 2)

## ------------------------------------------------------------------------
xx.init <- xx

## ---- echo = FALSE, fig.width = 7----------------------------------------
nyr <- 20

# results vectors
hsize.proj <- csurv.proj <- sr.proj <-  numeric(nyr)

# record initial state
xx <- xx.init

xx <- survival(xx)
xx <- infanticide(xx)

sr.proj[1]    <- sex_ratio(xx)
hsize.proj[1] <- xx@harem.size
csurv.proj[1] <- xx@realised.survival.rate[1]

xx <- transition(xx)

# project forward to adjust numbers
for (j in 2:nyr) {
    
    xx <- survival(xx)
    xx <- infanticide(xx)
    
    # record harem size
    # and cub survival
    sr.proj[j]    <- sex_ratio(xx)
    hsize.proj[j] <- xx@harem.size
    csurv.proj[j] <- xx@realised.survival.rate[1]
        
    xx <- transition(xx)
    
}

dfr <- rbind(
    data.frame(year = 1:length(hsize.proj), label = "Harem size", value = hsize.proj, model = "Non-equ.\nstart"),
    data.frame(year = 1:length(csurv.proj), label = "Cub surivival", value = csurv.proj, model = "Non-equ.\nstart"),
    data.frame(year = 1:length(sr.proj),    label = "Sex ratio", value = sr.proj, model = "Non-equ.\nstart"))


ggplot(dfr) + geom_line(aes(year, value, col = model)) + facet_grid(label~., scales = "free_y") + labs(y = '', x = 'Time', col = '')


## ------------------------------------------------------------------------
xx.equ <- equilibrium(xx)

## ---- echo = FALSE, fig.width = 7----------------------------------------
#xx <- equilibrium(xx, target.harem.size = 1.14, target.cub.survival = 0.327, verbose = TRUE)
hsize.proj[] <- harem(xx.equ)@harem.size
csurv.proj[] <- infanticide(survival(xx.equ))@realised.survival.rate[1]

# plot median probability of infanticide
dfr <- rbind(dfr, data.frame(year = 1:nyr, label = "Harem size", value = hsize.proj, model = "Equ.\nstart"), data.frame(year = 1:nyr, label = "Cub surivival", value = csurv.proj, model = "Equ.\nstart"))
ggplot(dfr) + geom_line(aes(year, value, col = model)) + facet_grid(label~., scales = "free_y") + labs(y = '', x = 'Time', col = '')


## ---- warning=FALSE------------------------------------------------------
harvest.rate <- seq(0, 0.4, length = 5)
selectivity  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)

n.births <- array(0, dim = c(length(harvest.rate), nyr))
cub.surv <- array(0, dim = c(length(harvest.rate), nyr))

# loop over harvest rates
for (i in 1:length(harvest.rate)) {
    
    # re-initialise
    yy <- xx.equ
    
    yy <- survival(yy)
    yy <- infanticide(yy)
    
    n.births[i, 1] <- yy[1]
    cub.surv[i, 1] <- yy@realised.survival.rate[1]
    
    yy <- transition(yy)
    
    # project forward to adjust numbers
    for (j in 2:nyr) {
        
        # kill individuals
        #yy[10:14] <- yy[10:14] - yy[10:14] * harvest.rate[i]
        
        yy <- survival(x = yy, y = yy * selectivity * harvest.rate[i])
        yy <- infanticide(yy)
        
        # record probability of infanticide
        # and cub survival
        n.births[i, j] <- yy[1]
        cub.surv[i, j] <- yy@realised.survival.rate[1]
            
        yy <- transition(yy)
        
    }
}

## ---- fig.width=6, echo = FALSE------------------------------------------
# plot median probability of infanticide
dimnames(n.births) <- list(H = harvest.rate, Year = 1:nyr)
ggplot(melt(n.births)) + geom_line(aes(Year, value, col = as.factor(H))) + 
    ggtitle('Number of births per year') + labs(y = '', col = 'Harvest\nRate')

# plot median cub survivorship
dimnames(cub.surv) <- list(H = harvest.rate, Year = 1:nyr)
ggplot(melt(cub.surv)) + geom_line(aes(Year, value, col = as.factor(H))) + 
    ggtitle('Expected cub survival') + labs(y = '', col = 'Harvest\nRate')

## ---- echo = FALSE, fig.width=6, warning=FALSE---------------------------
harvest.rate <- -seq(0.4, 0.0, length = 5)
selectivity  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)

n.births <- array(0, dim = c(length(harvest.rate), nyr))
cub.surv <- array(0, dim = c(length(harvest.rate), nyr))

# loop over harvest rates
for (i in 1:length(harvest.rate)) {
    
    # re-initialise
    yy <- xx.equ
    
    yy <- survival(yy)
    yy <- infanticide(yy)
    
    n.births[i, 1] <- yy[1]
    cub.surv[i, 1] <- yy@realised.survival.rate[1]
    
    yy <- transition(yy)
    
    # project forward to adjust numbers
    for (j in 2:nyr) {
        
        yy <- birth(yy)
        
        # kill individuals
        
        yy <- survival(x = yy, y = yy * selectivity * harvest.rate[i])
        yy <- infanticide(yy)
        
        # record probability of infanticide
        # and cub survival
        n.births[i, j] <- yy[1]
        cub.surv[i, j] <- yy@realised.survival.rate[1]
            
        yy <- transition(yy)
        
    }
}

# plot median probability of infanticide
dimnames(n.births) <- list(H = harvest.rate, Year = 1:nyr)
ggplot(melt(n.births)) + geom_line(aes(Year, value, col = as.factor(H))) + 
    ggtitle('Number of births per year') + labs(y = '', col = 'Immigration\nRate')

# plot median cub survivorship
dimnames(cub.surv) <- list(H = harvest.rate, Year = 1:nyr)
ggplot(melt(cub.surv)) + geom_line(aes(Year, value, col = as.factor(H))) + 
    ggtitle('Expected cub survival') + labs(y = '', col = 'Harvest\nRate')

