

# female birth rates per age category
k <- c(1.857143, 1.842105, 1.857143, 1.923077, 1.845070)


# harem size
h <- 1.5

# numbers
nf <- sum(c(99,80,73,45,194))
nm <- sum(c(41,24,24,28,95))
nh <- nf / h


# male encounters per harem
encounter.h <- function(nm, nh) {
  
  2 * nm / (nm + nh)
  
}

# male encounters
encounter.n <- function(nm, nh) {
  
  (2 * nm / (nm + nh)) * nh

}


scaler <- seq(0, 2, length = 101)

x1 <- x2 <- p <- numeric(101)

for (i in 1:101) {
  
  males   <- nm * scaler[i]
  harem <- nh * (max(scaler) - scaler[i])
  
  p[i] <- males
  x1[i] <- encounter.h(males, harem) 
  x2[i] <- encounter.n(males, harem) 
  
}

plot(x1 ~ p, type = 'l', xlab = 'number males', ylab = 'male encounters per harem')

plot(x2 ~ p, type = 'l', xlab = 'number males', ylab = 'number of encounters')

######

library(ggplot2)

source('utils/pdfr.r')

# numbers
# from Swanepoel et al. 2014
nfemales <- sum(c(99,80,73,45,194))
nmales   <- sum(c(41,24,24,28,95))

# harem size
harem_size <- 1.5

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
  # females mate 2 times per year
  nencounters_max <- nharems * 2
  
  # truncated at <= nencounters_max
  nencounters <- vapply(nencounters, function(x) min(x, nencounters_max), numeric(1))
  
  # return
  return(nencounters)
}

dfr$ne <- apply(dfr, 1, function(x) number_encounters(x[1], x[2]))

ggplot(dfr) + geom_line(aes(nm, ne, col = as.factor(nf))) + ggtitle('number of encounters')


######




