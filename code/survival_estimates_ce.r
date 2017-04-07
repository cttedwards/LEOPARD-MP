
###########################################################################################
# SURVIVAL ANALYSES
###########################################################################################

library(survival)


source('utils/reader.r')
source('utils/saver.r')
source('utils/writer.r')
source('utils/loader.r')

sabi.surv.cj <- reader("sabi.surv.cub.juv")
sabi.surv    <- reader("sabi.surv")

###########################################################################################
# ALL DATA (SEPARATED BY SEX)
###########################################################################################

sabi.surv.f <- sabi.surv[ which( ! sabi.surv$SEX %in% "MALE") , ]     # remove male records
sabi.surv.m <- sabi.surv[ which( ! sabi.surv$SEX %in% "FEMALE") , ]   # remove female records

sabi.surv.f.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.DEATH.MONTHS,event=CENSUS.STATUS,type='counting')~1, data = sabi.surv.f)
print(sabi.surv.f.fit,digits=6)
plot(sabi.surv.f.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="FEMALE SURVIVAL PROBABILITY", conf.int = TRUE)

sabi.surv.m.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.DEATH.MONTHS,event=CENSUS.STATUS,type='counting')~1, data = sabi.surv.m)
print(sabi.surv.m.fit,digits=6)
plot(sabi.surv.m.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="MALE SURVIVAL PROBABILITY", conf.int = TRUE)

###########################################################################################
# DEPENDENT DATA (SURVIVAL RATES DEPENDENT ON MATERALE AGE)
###########################################################################################

# cj survival (from all mothers)
sabi.surv.cj.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, data = sabi.surv.cj)
plot(sabi.surv.cj.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="All DEPENDENT SURVIVAL PROBABILITY", conf.int = TRUE)


###########################################################################################
# ESTIMATE SURVIVAL PROBABILITIES
###########################################################################################

# dependent young

# cub and juveniles (from all mothers)
surv.est.cj <- stepfun(sabi.surv.cj.fit$time, c(1, sabi.surv.cj.fit$surv))

plot(surv.est.cj)

age.month <- c(0,12,24)

abline(v = age.month, lty = 2)
abline(h = surv.est.cj(age.month), lty = 2)

surv.diff.cj <- surv.est.cj(age.month)

surv.cj <- surv.diff.cj[-1] / surv.diff.cj[1:(length(surv.diff.cj) - 1)] 

names(surv.cj) <- age.month[1:(length(age.month) - 1)]

# sub-adult and adults
age.month <- seq(24,230,12)

# females
surv.est.f <- stepfun(sabi.surv.f.fit$time, c(1, sabi.surv.f.fit$surv))

plot(surv.est.f)

abline(v = age.month, lty = 2)
abline(h = surv.est.f(age.month), lty = 2)

# proportion surviving to age x
surv.est.f <- surv.est.f(age.month)
names(surv.est.f) <- age.month
      
# proportion surviving age x
surv.f <- surv.est.f[-1] / surv.est.f[1:(length(surv.est.f) - 1)] 
names(surv.f) <- age.month[1:(length(age.month) - 1)]

# males
surv.est.m <- stepfun(c(sabi.surv.m.fit$time, 218), c(1, sabi.surv.m.fit$surv,0))

plot(surv.est.m)

abline(v = age.month, lty = 2)
abline(h = surv.est.m(age.month), lty = 2)

# proportion surviving to age x
surv.est.m <- surv.est.m(age.month)
names(surv.est.m) <- age.month

# proportion surviving age x
surv.m <- surv.est.m[-1] / surv.est.m[1:(length(surv.est.m) - 1)] 
names(surv.m) <- age.month[1:(length(age.month) - 1)]

###########################################################################################
# FINAL SURVIVAL ESTIMATES
###########################################################################################

surv <- c('sc0' = as.numeric(surv.cj['0']),
          'sc1' = as.numeric(surv.cj['12']), 
          'sf0' = as.numeric(surv.f['24']), 
          'sf+' = as.numeric(sum(surv.est.f[3:length(surv.est.f)]) / sum(surv.est.f[2:length(surv.est.f)])), 
          'sm0' = as.numeric(surv.m['24']), 
          'sm1' = as.numeric(surv.m['36']),
          'sm2' = as.numeric(surv.m['48']), 
          'sm+' = as.numeric(sum(surv.est.m[5:length(surv.est.m)]) / sum(surv.est.m[4:length(surv.est.m)])))

saver(surv, name = 'survival_estimates_ce')

###########################################################################################
# END
###########################################################################################
