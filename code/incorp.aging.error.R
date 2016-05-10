
source('aging_error.r')

#####################

#vec <- c(5, 0, 1, 2, 3, 4, 5, 2, 5, 8, 5, 3, 10, 5)
vec <- removals$trophy@kills
vec

vec.sa.male <- vec
vec.sa.male

vec.male36_48 <- vec
vec.male36_48

vec.male48_60 <- vec
vec.male48_60

vec.male60_72 <- vec
vec.male60_72

vec.male72_84 <- vec
vec.male72_84

vec.male84 <- vec
vec.male84

vec.female <- vec
vec.female

#####################

#----sub-adult male
if(vec[9] > 0){
  
  subadult_male.vec.final <- rep(0, 12)
  for(i in 1:vec[9]){
    subadult_male.vec.single <- c(0,0,0,0,0,0,1,0,0,0,0,0)
    subadult_male.vec.single.random <- sample(subadult_male.vec.single, length(subadult_male.vec.single), prob = c(aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], 1-aging.error[4], aging.error[4], aging.error[4], aging.error[4], aging.error[4], aging.error[4]))
    subadult_male.vec.final <- subadult_male.vec.final + subadult_male.vec.single.random
  }
  vec.sa.male[3:14] <- subadult_male.vec.final
} else {
  vec <- vec
}

#----adult male (36-48)
if(vec[10] > 0){
  
  male36_48.vec.final <- rep(0, 12)
  for(i in 1:vec[10]){
    male36_48.vec.single <- c(0,0,0,0,0,0,0,1,0,0,0,0)
    male36_48.vec.single.random <- sample(male36_48.vec.single, length(male36_48.vec.single), prob = c(aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[5], 1-aging.error[5], aging.error[5], aging.error[5], aging.error[5], aging.error[5]))
    male36_48.vec.final <- male36_48.vec.final + male36_48.vec.single.random
  }
  vec.male36_48[3:14] <- male36_48.vec.final
} else {
  vec <- vec
}

#----adult male (48-60)
if(vec[11] > 0){
  
  male48_60.vec.final <- rep(0, 12)
  for(i in 1:vec[11]){
    male48_60.vec.single <- c(0,0,0,0,0,0,0,0,1,0,0,0)
    male48_60.vec.single.random <- sample(male48_60.vec.single, length(male48_60.vec.single), prob = c(aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[6], aging.error[6], 1-aging.error[6], aging.error[6], aging.error[6], aging.error[6]))
    male48_60.vec.final <- male48_60.vec.final + male48_60.vec.single.random
  }
  vec.male48_60[3:14] <- male48_60.vec.final
} else {
  vec <- vec
}

#----adult male (60-72)
if(vec[12] > 0){
  
  male60_72.vec.final <- rep(0, 12)
  for(i in 1:vec[12]){
    male60_72.vec.single <- c(0,0,0,0,0,0,0,0,0,1,0,0)
    male60_72.vec.single.random <- sample(male60_72.vec.single, length(male60_72.vec.single), prob = c(aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[7], aging.error[7], aging.error[7], 1-aging.error[7], aging.error[7], aging.error[7]))
    male60_72.vec.final <- male60_72.vec.final + male60_72.vec.single.random
  }
  vec.male60_72[3:14] <- male60_72.vec.final
} else {
  vec <- vec
}

#----adult male (72-84)
if(vec[13] > 0){
  
  male72_84.vec.final <- rep(0, 12)
  for(i in 1:vec[13]){
    male72_84.vec.single <- c(0,0,0,0,0,0,0,0,0,0,1,0)
    male72_84.vec.single.random <- sample(male72_84.vec.single, length(male72_84.vec.single), prob = c(aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[8], aging.error[8], aging.error[8], aging.error[8], 1-aging.error[8], aging.error[8]))
    male72_84.vec.final <- male72_84.vec.final + male72_84.vec.single.random
  }
  vec.male72_84[3:14] <- male72_84.vec.final
} else {
  vec <- vec
}

#----adult male (>84)
if(vec[14] > 0){
  
male84.vec.final <- rep(0, 12)
for(i in 1:vec[14]){
  male84.vec.single <- c(0,0,0,0,0,0,0,0,0,0,0,1)
  male84.vec.single.random <- sample(male84.vec.single, length(male84.vec.single), prob = c(aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[2], aging.error[9], aging.error[9], aging.error[9], aging.error[9], aging.error[9], 1-aging.error[9]))
  male84.vec.final <- male84.vec.final + male84.vec.single.random
  }
  vec.male84[3:14] <- male84.vec.final
} else {
  vec <- vec
}


final.male.vec <- (if(vec[9]>0){ 
  vec.sa.male[3:14] 
  } else {
  rep(0, 12)
  }) +
  (if(vec[10]>0){ 
    vec.male36_48[3:14] 
  } else {
    rep(0, 12)
  }) +
  (if(vec[11]>0){ 
    vec.male48_60[3:14] 
  } else {
    rep(0, 12)
  }) +
  (if(vec[12]>0){ 
    vec.male60_72[3:14] 
  } else {
    rep(0, 12)
  }) +
  (if(vec[13]>0){ 
    vec.male72_84[3:14] 
  } else {
    rep(0, 12)
  }) +
  (if(vec[14]>0){ 
    vec.male84[3:14] 
  } else {
    rep(0, 12)
  })

final.male.vec
final.male.vec <- c(0,0,final.male.vec)

#####################

#----females
if(sum(vec[3:8]) > 0){
  
  female.vec.final <- rep(0, 12)
  for(i in 1:sum(vec[3:8])){
    female.vec.single <- c(1,0,0,0,0,0,0,0,0,0,0,0)
    female.vec.single.random <- sample(female.vec.single, length(female.vec.single), prob = c(1-aging.error[1], 1-aging.error[1], 1-aging.error[1], 1-aging.error[1], 1-aging.error[1], 1-aging.error[1], aging.error[1], aging.error[1], aging.error[1], aging.error[1], aging.error[1], aging.error[1]))
    female.vec.final <- female.vec.final + female.vec.single.random
  }
  vec.female[3:14] <- female.vec.final
} else {
  vec <- vec
}

final.female.vec <- if(sum(vec[3:8]>0)){
  vec.female[3:14]
} else {
  rep(0, 12)
}

final.female.vec
final.female.vec <- c(vec[1:2], final.female.vec)

#####################

final.vec.all.leopard <- final.female.vec + final.male.vec
final.vec.all.leopard
sum(final.vec.all.leopard)
sum(vec)
removals$trophy@kills <- final.vec.all.leopard
