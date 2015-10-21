



###########################################################################################
# LOAD PACKAGES
###########################################################################################

#library(msm)

rm(list=ls())
source('utils/saver.r')
source('utils/loader.r')
source('utils/reader.r')
source('S4objects/S4classes.r')
source('S3objects/S3methods.r')

loader('survival_estimates')
loader('hunt.stats')

###########################################################################################
# SET WORKING DIRECTORY/SOURCE & LOAD PARAMETERS 
###########################################################################################

# parameter vector of vital rates
param <- c(S.1   = surv.diff.c[2] * (1/surv.diff.c[1]),       # cub
           S.2   = surv.diff.j[2] * (1/surv.diff.j[1]),       # juvenile
           S.3   = surv.diff.f[4] * (1/surv.diff.f[3]),       # sub-adult female
           S.4   = surv.diff.f[5] * (1/surv.diff.f[4]),       # adult female (36-48)
           S.5   = surv.diff.f[6] * (1/surv.diff.f[5]),       # adult female (48-60)
           S.6   = surv.diff.f[7] * (1/surv.diff.f[6]),       # adult female (60-72)
           S.7   = surv.diff.f[8] * (1/surv.diff.f[7]),       # adult female (72-84)
           S.8   = surv.diff.f[9] * (1/surv.diff.f[8]),       # adult female (>84)
           S.9   = surv.diff.m[4] * (1/surv.diff.m[3]),       # sub-adult male
           S.10  = surv.diff.m[5] * (1/surv.diff.m[4]),       # adult male (36-48)
           S.11  = surv.diff.m[6] * (1/surv.diff.m[5]),       # adult male (48-60)
           S.12  = surv.diff.m[7] * (1/surv.diff.m[6]),       # adult male (60-72)
           S.13  = surv.diff.m[8] * (1/surv.diff.m[7]),       # adult male (72-84)
           S.14  = surv.diff.m[9] * (1/surv.diff.m[8]),       # adult male (>84)
           br.36 = 1.857142857,                               # birth rate (mother 36-48)     
           br.48 = 1.842105263,                               # birth rate (mother 48-60)              
           br.60 = 1.857142857,                               # birth rate (mother 60-72)               
           br.72 = 1.923076923,                               # birth rate (mother 72-84)               
           br.84 = 1.845070423,                               # birth rate (mother >84)                
           p     = 0.5)                                       # sex ratio at birth

# multiplicative maternal effects on cub/juvenile survival
maternal.effects <- c(S.1.1 = (surv.diff.c1[2] * (1/surv.diff.c1[1])) / (surv.diff.c[2] * (1/surv.diff.c[1])),
                      S.1.2 = (surv.diff.c2[2] * (1/surv.diff.c2[1])) / (surv.diff.c[2] * (1/surv.diff.c[1])),
                      S.1.3 = (surv.diff.c3[2] * (1/surv.diff.c3[1])) / (surv.diff.c[2] * (1/surv.diff.c[1])),
                      S.1.4 = (surv.diff.c4[2] * (1/surv.diff.c4[1])) / (surv.diff.c[2] * (1/surv.diff.c[1])),
                      S.1.5 = (surv.diff.c5[2] * (1/surv.diff.c5[1])) / (surv.diff.c[2] * (1/surv.diff.c[1])),
                      S.2.1 = (surv.diff.j1[2] * (1/surv.diff.j1[1])) / (surv.diff.j[2] * (1/surv.diff.j[1])),
                      S.2.2 = (surv.diff.j2[2] * (1/surv.diff.j2[1])) / (surv.diff.j[2] * (1/surv.diff.j[1])),
                      S.2.3 = (surv.diff.j3[2] * (1/surv.diff.j3[1])) / (surv.diff.j[2] * (1/surv.diff.j[1])),
                      S.2.4 = (surv.diff.j4[2] * (1/surv.diff.j4[1])) / (surv.diff.j[2] * (1/surv.diff.j[1])),
                      S.2.5 = (surv.diff.j5[2] * (1/surv.diff.j5[1])) / (surv.diff.j[2] * (1/surv.diff.j[1]))
)

# initial population size (Swanepoel et al. 2014; ±2000 leopard)
x.initial <- c(nc=694,
               nj=410,
               saf=117,
               f36=99,
               f48=80,
               f60=73,
               f72=45,
               f84=194,
               sam=84,
               m36=41,
               m48=24,
               m60=24,
               m72=28,
               m84=95)

# create new object to hold leopard numbers
# and vital rates
xx <- leopard(x.initial, param[1:14], param[15:19])

# check assignments
all(xx@.Data == x.initial)

all(xx@expected.survival.rate == param[1:14])
all(xx@female.birth.rate == param[15:19])

# assign multiplicative maternal effects
xx@maternal.effect[] <- matrix(maternal.effects, nrow=2, ncol=5, byrow=T)

# initial assignment of cubs at random to maternal age class
prob.maternal.age <- xx[4:8]/sum(xx[4:8])
xx@realised.birth[1,] <- rmultinom(1, xx[1], prob.maternal.age)
xx@realised.birth[2,] <- rmultinom(1, xx[2], prob.maternal.age)

# check assignments
all(apply(xx@realised.birth, 1, sum) == xx@.Data[1:2])

# calculate stochastic survival
## 1. calculate juvenile survival
## using juveniles in each maternal age
## category from previous time step [stored 
## in xx@realised.birth]. Maternal age dependent
## survival is calculated using the
## xx@maternal.effect multiplier
## 2. update juvenile numbers in 
## xx@realised.birth using cubs
## in each maternal age category and
## maternal age dependent survival
## 3. calculate overall cub survival from
## this transition step
## 4. update stochastic survival rates for 
## remaining age classes
xx <- survival(xx)

# check that realised survival rates return integer values
# [this works within the function but not here for some reason]
all((xx@realised.survival.rate * xx@.Data) %% 1 == 0)

# check cub transition to juveniles
sum(xx@realised.birth[2, ]) == xx@realised.survival.rate[1] * xx[1]

# check stochasticity
ss <- c()
for (i in 1:1000) ss <- c(ss, survival(xx)@realised.survival.rate)
ss <- matrix(ss, ncol = 14, byrow = TRUE)

for (i in 1:14) {
    plot(ss[,i], type = 'l', main = names(xx)[i]) 
    abline(h = xx@expected.survival.rate[i], col = 2)
    abline(h = mean(ss[,i]), col = 4)
}

# construct transition matrix
# for survival only
M <- tmatrix(xx)

# checks
all(apply(M, 2, sum) == xx@realised.survival.rate)

sum(M %*% xx) == sum(xx * xx@realised.survival.rate)


# calculate stochastic birth
## execute stochastic birth process
## and record actual number of births
## per maternal age category in 
## xx@realised.birth['cub',]
xx <- birth(xx)

# check
all(xx@realised.birth['cub',] == xx[4:8] * xx@realised.birth.rate)

# check stochasticity
bb <- c()
for (i in 1:1000) bb <- c(bb, birth(xx)@realised.birth.rate)
bb <- matrix(bb, ncol = 5, byrow = TRUE)

for (i in 1:5) {
    plot(bb[,i], type = 'l', main = names(xx@expected.birth.rate)[i]) 
    abline(h = xx@expected.birth.rate[i], col = 2)
    abline(h = mean(bb[,i]), col = 4)
}


# construct transition matrix
M <- tmatrix(xx)

sum(M[1,]) == sum(xx@realised.birth.rate)

# step forward
## 1. construct transistion matrix using
## xx@realised.survival.rate and
## xx@realised.birth.rate and
## matrix multiply to update the numbers
## vector
## 2. reset realised rates to zero
## [NB the numbers in xx[1] and xx[2]
## should now match those in xx@realised.birth. You
## can check this using apply(xx@realised.birth,1,sum)]

xx[] <- as.integer(M %*% xx)

# check
xx[1] == sum(xx@realised.birth[1,])
xx[2] == sum(xx@realised.birth[2,])

