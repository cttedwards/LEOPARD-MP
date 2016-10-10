
rm(list=ls())

h <- 9         # average number of breeding females per male
k <- 1.8       # average birth rate
nf <- 100      # number of females
nm <- 1       # number of males

(k * nm) / (nm + (nf * h^-1))  # harmonic-mean birth function

###########################################################################################

# per capita female fecundity (with increasing males, females held constant)

f.func <- function(h = 5, k = 1.8, nf = 100, nm = 1) {
  
  #k <- rpois(1, 1.8)
    
  num <- k * nm
  
  den <- (nm + (nf * h^-1))
  
  Ff <- num / den
  
  return(Ff)
  
}

x <- numeric(200)
males <- seq(0, 200, by = 1)

for(i in 1:200) {
  out   <- f.func(nm = males[i])
  x[i] <- out  
}

plot(x, xlab='number of males : female', ylab='female birth rate', type = 'l')

###########################################################################################

# per capita male fecundity (with increasing females, males held constant)

m.func <- function(h = 5, k = 1.8, nf = 1, nm = 100) {
  
  num <- k * nf
  
  den <- (nm + (nf * h^-1))  
  
  Mf <- num / den
  
  return(Mf)
  
}

x <- numeric(200)
females <- seq(0, 200, by = 1)

for(i in 1:200) {
  out   <- m.func(nf = females[i])
  x[i] <- out  
}

plot(x, xlab='number of females : male', ylab='male birth rate', type = 'l')

###########################################################################################
###########################################################################################

# total births for increasing males with constant females

tb.m.func <- function(h = 5, k = 1.8, nf = 100, nm = 1) {
  
  num <- k * nm
  
  den <- (nm + (nf * h^-1))
  
  Ff <- num / den
  
  tb <- Ff * nf
  
  return(tb)
  
}

x <- numeric(100)
males <- seq(0, 100, by = 1)

for(i in 1:100) {
  out   <- tb.m.func(nm = males[i])
  x[i] <- out  
}

plot(x, xlab='number of males : female', ylab='total births', type = 'l')


###########################################################################################

# total births for increasing females with constant males

tb.f.func <- function(h = 5, k = 1.8, nf = 1, nm = 100) {
  
  num <- k * nf
  
  den <- (nm + (nf * h^-1))
  
  Mf <- num / den
  
  tb <- Mf * nf
  
  return(tb)
  
}

x <- numeric(100)
females <- seq(0, 100, by = 1)

for(i in 1:100) {
  out   <- tb.f.func(nf = females[i])
  x[i] <- out  
}

plot(x, xlab='number of females : male', ylab='total births', type = 'l')

###########################################################################################





