
#vec <- c(rep(0, 13), 5)
vec <- removals$trophy@kills
vec

vec.male84 <- vec
vec.male84

#----adult male (>84)
if(vec[14] > 0){
  
  male84.vec.final <- rep(0, 12)
  
  for(t in 1:vec[14]){
    male84.vec.single.random <- as.vector(rmultinom(1, size = 1, prob=c(rep(0, 6), rep(1-0.91, 5), 0.91)))
    male84.vec.final <- male84.vec.final + male84.vec.single.random
  }
  vec.male84[3:14] <- male84.vec.final
} else {
  vec <- vec
}


final.male.vec <- (if(vec[14]>0){ 
    vec.male84[3:14] 
  } else {
    rep(0, 12)
  })

final.male.vec <- c(0,0,final.male.vec)
final.male.vec

removals$trophy@kills <- final.male.vec
