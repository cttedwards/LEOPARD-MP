
#{{{ kill function
# generic function
"kill" <- function(x, ...) UseMethod("kill")

# method
"kill.leopard" <- function(x, catchability = 0.001, theta = 0.4) {
    
    ilogit <- function(x) exp(x) / (exp(x) + 1)
    
    numbers.at.age <- x@.Data
    
    # hunter selectivity (Braczkowski et al 2015)
    # the probability of shooting each leopard age class if encountered
    hunter.preference <- c(0.256,0.244,0.712,0.462,0.462,0.462,0.359,0.726,0.712,0.897,0.974,0.949,1.00,0.99)
    
    # assuming independent encounters for each age class
    # the probability of encounter at age is a function of
    # the numbers at age and the catchability
    # (this is a multinomial regression)
    prob.encounter.at.age <- ilogit(log(catchability) + theta * log(numbers.at.age))
    
    # the probability of a kill is the probability of encouter multiplied
    # by the conditional probability of killing (i.e. the probability that 
    # a hunter will kill an individual of that age given that it has been 
    # encountered)
    prob.kill.at.age <- hunter.preference * prob.encounter.at.age
    
    # probability of not killing
    prob.not.kill.at.age <- 1 - prob.kill.at.age
    
    # using the probability of kill at age select an
    # individual to kill by sampling from a multinomial
    # distribution
    kill.at.age <- as.logical(rmultinom(1, 1, prob.kill.at.age))
    
    # the binomial probability of killing that one individual and 
    # not killing any of the others is:
    prob.kill <- prob.kill.at.age[kill.at.age] * prod(prob.not.kill.at.age[!kill.at.age]) * choose(14,1)
    
    # record waiting time
    waiting.time <- rgeom(1, prob.kill) + 1
    
    # return
    return(list(kill.at.age = kill.at.age, waiting.time = waiting.time))
    
}
#}}}

#{{{ survival function
# generic function
"survival" <- function(x, ...) UseMethod("survival")

# method
"survival.leopard" <- function(x) {
              
              survival.rates  <- x@expected.survival.rate
              
              maternal.effect <- x@maternal.effect
              
              # if no females in that age category
              # then cubs and juveniles do not survive
              number.females <- x@.Data[4:8]
              if (any(number.females == 0))
                  maternal.effect[which(number.females == 0),] <- 0
              
              # calculate weighted average juvenile survival using 
              # the number of juveniles in each maternal age class
              number.juvenile   <- x@realised.birth['juvenile',]
              juvenile.survival <- survival.rates[2] * maternal.effect['juvenile',]
              if (sum(number.juvenile) > 0) {
                  survival.rates[2] <- sum(juvenile.survival * number.juvenile) / sum(number.juvenile)
              } else {
                  stop('no juveniles\n')
              }
              
              # update juvenile survival assuming stochastic sampling process
              survival.rates[2] <- ifelse(x[2]  > 0, rbinom(1, x[2],  survival.rates[2])  / x[2] , 0)
              
              # next calculate survival for cubs
              # in each maternal age class
              cub.survival  <- survival.rates[1] * maternal.effect['cub',]
              
              # sample cubs and allocate to juvenile age category
              number.cubs   <- as.integer(x@realised.birth['cub',])
              x@realised.birth['juvenile',] <- rbinom(5, number.cubs, cub.survival)
              
              # now calculate the average survival rate associated with 
              # the stochastic process just implemented
              if (sum(number.cubs) > 0) {
                survival.rates[1] <- sum(x@realised.birth['juvenile',]) / sum(x@realised.birth['cub',])
              } else {
                stop('no cubs\n')
              }
              
              # implement stochastic survival for remaining age classes
              for (i in 3:14) {
                survival.rates[i] <- ifelse(x[i]  > 0, rbinom(1, x[i],  survival.rates[i])  / x[i] , 0)
                #cat(i,': ',as.integer(survival.rates[i] * x[i]) == survival.rates[i] * x[i], '\n')
              }    
              
              # finally record realised rates
              # [NB these should return integer values when applied
              # to the current numbers in each age class]
              x@realised.survival.rate <- survival.rates
              
              return(x)
              
          }
#}}}

#{{{ birth function
# generic function
"birth" <- function(x, ...) UseMethod("birth")

# method
# using harmonic-mean birth function (k = birth rates from Sabi Sands; h = 5 average number of females a male can monopolise)
"birth.leopard" <- function(x) {
            
            # number of females 
            # per maternal age
            nf <- x[4:8]
            
            # number of males
            nm <- sum(x[10:14])
    
            # number of harems
            h <- 5
            nh <- sum(nf) / h
            
            # expected births per female age category
            expected.birth <- 2 * nm * nf * x@litter.size / (nm + nh)
    
            # realised births per maternal age category
            for (a in 1:length(nf)) {
                x@realised.birth['cub', a] <- ifelse(nm > 0 & nh > 0, rpois(1, expected.birth[a]), 0)
            }
            
            # maternal birth rates
            # [NB realised rates should return integer values when applied
            # to the current numbers in each age class]
            x@expected.birth.rate <- expected.birth / nf
            x@realised.birth.rate <- x@realised.birth['cub', ] / nf
              
            return(x)
              
          }
#}}}

#{{{ transition matrix function
# generic function
"tmatrix" <- function(x, ...) UseMethod("tmatrix")

# method
"tmatrix.leopard" <- function(x, y) {
              
              S <- x@realised.survival.rate
              B <- x@realised.birth.rate
              
              # the sex of each cub is
              # a binomial process
              p <- ifelse(x['nj'] > 0, rbinom(1, as.integer(x['nj'] * S[2]), x@sex.ratio) / (x['nj'] * S[2]), 0)
              
              # empty transition matrix
              M <- array(0, dim = c(length(x), length(x)))
              
              # populate matrix with survival rates
              
              # cub
              M[2,1] <- S[1]
              
              # sub-adult
              M[3,2] <- (1 - p) * S[2]  # juvenile to sub-adult female
              M[9,2] <- p * S[2]      # juvenile to sub-adult male
              
              # adult
              M[4,3]   <- S[3]        # sub-adult female
              M[5,4]   <- S[4]        # adult female (36-48)
              M[6,5]   <- S[5]        # adult female (48-60)
              M[7,6]   <- S[6]        # adult female (60-72)
              M[8,7]   <- S[7]        # adult female (72-84)
              M[8,8]   <- S[8]        # adult female (>84)
              M[10,9]  <- S[9]        # sub-adult male
              M[11,10] <- S[10]       # adult male (36-48)
              M[12,11] <- S[11]       # adult male (48-60)
              M[13,12] <- S[12]       # adult male (60-72)
              M[14,13] <- S[13]       # adult male (72-84)
              M[14,14] <- S[14]       # adult male (>84)
              
              # populate matrix with birth rates
              M[1,4] <- B[1]   # birth rate of female aged 36-48 months
              M[1,5] <- B[2]   # birth rate of female aged 48-60 months
              M[1,6] <- B[3]   # birth rate of female aged 60-72 months
              M[1,7] <- B[4]   # birth rate of female aged 72-84 months
              M[1,8] <- B[5]   # birth rate of female aged >84 months
              
              # return
              return(M)
          }
#}}}

#{{{ transition/dynamics function
# generic function
"transition" <- function(x, ...) UseMethod("transition")

# method
"transition.leopard" <- function(x, y) {
    
    M <- tmatrix(x)
    
    # project forward
    x[] <- as.integer(M %*% x)
    
    # re-set vital rates
    x@realised.survival.rate[] <- 0
    x@expected.birth.rate[]    <- 0
    x@realised.birth.rate[]    <- 0
    
    # return
    return(x)
}
#}}}


