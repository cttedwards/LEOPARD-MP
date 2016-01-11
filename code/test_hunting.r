
###########################################################################################
# LOAD PACKAGES
###########################################################################################

#library(msm)

rm(list=ls())
#setwd("/Users/RossTyzackPitman/Documents/OneDrive/PhD/Data/R_Database/R_PROJECTS/PhD_Chapter3/MSE_Paper/code")
source('utils/saver.r')
source('utils/loader.r')
source('utils/reader.r')
source('S4objects/S4classes.r')
source('S3objects/S3methods.r')

loader('survival_estimates')
loader('hunt.stats')

s1 <- Sys.time()

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

# parameter vector of standard error values
param.se <- c(S.1   = surv.diff.c.se[2],
              S.2   = surv.diff.j.se[2],
              S.3   = surv.diff.f.se[4],
              S.4   = surv.diff.f.se[5],
              S.5   = surv.diff.f.se[6],
              S.6   = surv.diff.f.se[7],
              S.7   = surv.diff.f.se[8],
              S.8   = surv.diff.f.se[9],
              S.9   = surv.diff.m.se[4],
              S.10  = surv.diff.m.se[5],
              S.11  = surv.diff.m.se[6],
              S.12  = surv.diff.m.se[7],
              S.13  = surv.diff.m.se[8],
              S.14  = surv.diff.m.se[9],
              br.36 = 0.08,                    
              br.48 = 0.07,                   
              br.60 = 0.12,                    
              br.72 = 0.18,                    
              br.84 = 0.06,                    
              p     = 0)


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

# dimensions

## number of monte carlo samples
nreps <- 1000
## number of projection years
nyr.proj <- 50

## matrix of paramter values
params <- matrix(param,nrow=length(param),ncol=nreps)  
rownames(params) <- names(param)                       
colnames(params) <- 1:nreps

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

# population projection array
x <- array(x.initial,dim=c(length(x.initial),nreps,nyr.proj))
dimnames(x) <- list(age.class = names(x.initial),rep = 1:nreps, year = 1:nyr.proj)
x[,,2:nyr.proj] <- NA

# projection
for(i in 1:nreps) {
    
    # sample parameters
    param.sample <- params[,i]
    
    # create new object to hold leopard numbers
    # and vital rates
    xx <- newleopard(x.initial, param.sample[1:14], param.sample[15:19])
    
    # assign multiplicative maternal effects
    xx@maternal.effect[] <- matrix(maternal.effects, nrow=2, ncol=5, byrow=T)
    
    # assign cubs at random to maternal age class
    prob.maternal.age <- xx[4:8]/sum(xx[4:8])
    xx@realised.birth[1,] <- rmultinom(1, xx[1], prob.maternal.age)
    xx@realised.birth[2,] <- rmultinom(1, xx[2], prob.maternal.age)
    
    # data.frame to record kills for one rep
    kill.summary <- data.frame(year = numeric(), age = character(), waiting.time = numeric())
    
    # loop forward over years
    for(y in 2:nyr.proj) {
        
        # correlated deviation in survival: 
        # log-normal with cv = 0.2
        # truncated at 0 and 1
        sdev <- rnorm(1)
        sigma <- sqrt(log(1+0.20^2))
        
        # (check this: if you take 1000 values of sdev the right hand side should have a mean value
        # approximately equal to param.sample)
        param.sample[1:14] <- exp(log(param.sample[1:14]) +  sigma * sdev - sigma^2/2)
        param.sample[1:14] <- vapply(vapply(param.sample[1:14],function(x) max(x,0),numeric(1)),function(x) min(x,1),numeric(1))
        
        # quota setting function
        # (to be updated)
        quota <- 1
        
        # how will tenure rate affect cub survival??????? Infanticide
        
        ## kill
        kills <- vector('list', quota)
        for(k in 1:quota) {
            
            # stochastic kill with associated
            # waiting time
            kills[[k]] <- kill(xx)
            
            # remove individual from population
            # numbers vector
            xx <- xx - kills[[k]]$kill.at.age
            
            # if you kill a cub or juvenile then this will
            # have an effect on the average cub/juvenile survival
            # rate because of the maternal age effect. Because the kill
            # function does not descriminate the maternal age when
            # hunting cubs/juveniles (why should it?) we just remove
            # them at random [NB over lots of iterations this is
            # the same as doing nothing because the average survival
            # will be unchanged]
            if(kills[[k]]$kill.at.age[1]) {
                kill.prob <- xx@realised.birth[1,]/sum(xx@realised.birth[1,])
                xx@realised.birth[1,] <- xx@realised.birth[1,] - rmultinom(1, 1, kill.prob)
            }
            if(kills[[k]]$kill.at.age[2]) {
                kill.prob <- xx@realised.birth[2,]/sum(xx@realised.birth[2,])
                xx@realised.birth[2,] <- xx@realised.birth[2,] - rmultinom(1, 1, kill.prob)
            }
            
            # update summary data.frame
            kill.summary <- rbind(kill.summary, data.frame(year = y, age = xx@names[which(kills[[k]]$kill.at.age)], waiting.time = kills[[k]]$waiting.time))
        }
        
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
        ## [NB this won't work for first time step
        ## because initial cub numbers are not assigned
        ## to a maternal age class]
        xx <- survival(xx)
        
        # calculate stochastic birth
        ## execute stochastic birth process
        ## and record actual number of births
        ## per maternal age category in 
        ## xx@realised.birth['cub',]
        xx <- birth(xx)
        
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
        xx <- transition(xx)
                
        # record numbers
        x[,i,y] <- xx
                
    }
}

###########################################################################################
# PLOTS
###########################################################################################

x.tot <- apply(x, 2:3, sum)
boxplot(x.tot,ylab="N",xaxt="n",xlab="Step",outline=FALSE,ylim=c(0,5000))
axis(side=1,at=1:nyr.proj)

###########################################################################################
# COMPUTATION TIME
###########################################################################################

Sys.time()-s1

###########################################################################################
# END
###########################################################################################
