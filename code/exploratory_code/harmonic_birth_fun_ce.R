
rm(list=ls())

#h <- 5         # average number of breeding females per male
#k <- 1.8       # average birth rate
#nf <- 100      # number of females
#nm <- 1       # number of males

#k * nm / (nm + (nf * h^-1))  # harmonic-mean birth function

###########################################################################################

# female birth rates per age category
k <- c(1.857143, 1.842105, 1.857143, 1.923077, 1.845070)

# numbers
nf <- c(99,80,73,45,194)
nm <- c(41,24,24,28,95)

# harem size
h <- 5

# per capita female fecundity

f.func <- function(nf, nm) {
    
    k * sum(nm) / (sum(nm) + sum(nf/h))
}

# per capita male fecundity

m.func <- function(nf, nm) {
  
    sum(k * nf) / (sum(nm) + sum(nf/h))
}

# total births

b.func1 <- function(nf, nm) {
  
    Ff <- f.func(nf, nm)
    Fm <- m.func(nf, nm)
    
    births <- sum(nf * Ff) + sum(nm) * Fm
    
    return(births) 
  
}

b.func2 <- function(nf, nm) {
    
    births <- sum(nf * k) * 2 * sum(nm) / (sum(nm) + sum(nf/h))
    
    return(births) 
    
}

b.func3 <- function(nf, nm) {
    
    births <- nf * k * 2 * sum(nm) / (sum(nm) + sum(nf/h))
    
    return(births) 
    
}

scaler <- seq(0, 2, length = 101)

x1 <- x2 <- x3 <- p <- numeric(101)

for (i in 1:101) {
    
    males   <- nm * scaler[i]
    females <- nf * (max(scaler) - scaler[i])
    
    p[i] <- sum(males)/sum(males + females)
    x1[i] <- b.func1(females, males) 
    x2[i] <- b.func2(females, males) 
    x3[i] <- sum(b.func3(females, males)) 
}

plot(x1 ~ p, type = 'l', xlab = 'proportion males', ylab = 'total births')
lines(x2 ~ p, col = 2)
lines(x3 ~ p, col = 4)

for (i in 1:101) {
    
    males   <- nm * scaler[i]
    females <- nf
    
    p[i] <- sum(males)/sum(females)
    x1[i] <- b.func1(females, males)  
}

plot(x1 ~ p, xlab = 'number of males per female', ylab = 'total births', type = 'l')

for (i in 1:101) {
    
    males   <- nm
    females <- nf * scaler[i]
    
    p[i] <- sum(females)/sum(males)
    x1[i] <- b.func1(females, males)  
}

plot(x1 ~ p, xlab = 'number of females per male', ylab = 'total births', type = 'l')
