
###########################################################################################
# SURVIVAL ANALYSES
###########################################################################################

library(OIsurv)

rm(list=ls())
setwd("/Users/RossTyzackPitman/Documents/OneDrive/PhD/Data/R_Database/R_PROJECTS/PhD_Chapter3/MSE_Paper/LEOPARD-MP/code")

source('utils/reader.r')
source('utils/saver.r')
source('utils/writer.r')
source('utils/loader.r')

kzn.surv.cub = reader("kzn.surv.cub")
kzn.surv.juv = reader("kzn.surv.juv")
kzn.surv.sa = reader("kzn.surv.sa")
kzn.surv.a = reader("kzn.surv.a")
kzn.surv.cj = reader("kzn.surv.cub.juv")
kzn.surv = reader("kzn.surv")

###########################################################################################
# ALL DATA (SEPARATED BY SEX)
###########################################################################################

kzn.surv.f <- kzn.surv[ which( ! kzn.surv$SEX %in% "MALE") , ]     # remove male records
kzn.surv.m <- kzn.surv[ which( ! kzn.surv$SEX %in% "FEMALE") , ]   # remove female records

kzn.surv.f.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.DEATH.MONTHS,event=CENSUS.STATUS,type='counting')~1, data = kzn.surv.f)
print(kzn.surv.f.fit,digits=6)
plot(kzn.surv.f.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="FEMALE SURVIVAL PROBABILITY", conf.int = TRUE)

kzn.surv.m.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.DEATH.MONTHS,event=CENSUS.STATUS,type='counting')~1, data = kzn.surv.m)
print(kzn.surv.m.fit,digits=6)
plot(kzn.surv.m.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="MALE SURVIVAL PROBABILITY", conf.int = TRUE)

###########################################################################################
# DEPENDENT DATA OF CUBS AND JUVENILES | we cannot estimate cub and juvenile surv based
# on mothers age since we do not have the data from KZN (as we do with Sabi Sands)
###########################################################################################

kzn.surv.c.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS,event=AGE.MET,type='counting')~1, 
                           data = kzn.surv.cub)                                                              # cj survival (from all mothers)
plot(kzn.surv.c.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="CUB SURVIVAL PROBABILITY", conf.int = TRUE)


kzn.surv.j.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS,event=AGE.MET,type='counting')~1, 
                          data = kzn.surv.juv)                                                              # cj survival (from all mothers)
plot(kzn.surv.j.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="JUVENILE SURVIVAL PROBABILITY", conf.int = TRUE)

###########################################################################################
# ESTIMATE SURVIVAL PROBABILITIES
###########################################################################################

### dependent young

# cub (from all mothers)
surv.est.c <- stepfun(kzn.surv.c.fit$time, c(1, kzn.surv.c.fit$surv))
surv.diff.c <- surv.est.c(c(0,12))
surv.est.c.se <- stepfun(kzn.surv.c.fit$time, c(1, kzn.surv.c.fit$std.err))
surv.diff.c.se <- surv.est.c.se(c(0,12))

# juvenile (from all mothers)
surv.est.j <- stepfun(kzn.surv.j.fit$time, c(1, kzn.surv.j.fit$surv))
surv.diff.j <- surv.est.j(c(12,24))
surv.est.j.se <- stepfun(kzn.surv.j.fit$time, c(1, kzn.surv.j.fit$std.err))
surv.diff.j.se <- surv.est.j.se(c(12,24))

### sub-adult and adults

# females
surv.est.f <- stepfun(kzn.surv.f.fit$time, c(1, kzn.surv.f.fit$surv))
surv.diff.f <- surv.est.f(c(0,12,24,36,48,60,72,84,150))
surv.est.f.se <- stepfun(kzn.surv.f.fit$time, c(1, kzn.surv.f.fit$std.err))
surv.diff.f.se <- surv.est.f.se(c(0,12,24,36,48,60,72,84,150))

# males
surv.est.m <- stepfun(kzn.surv.m.fit$time, c(1, kzn.surv.m.fit$surv))
surv.diff.m <- surv.est.m(c(0,12,24,36,48,60,72,84,150))
surv.est.m.se <- stepfun(kzn.surv.m.fit$time, c(1, kzn.surv.m.fit$std.err))
surv.diff.m.se <- surv.est.m.se(c(0,12,24,36,48,60,72,84,150))

###########################################################################################
# FINAL SURVIVAL ESTIMATES
###########################################################################################

S.1   <- surv.diff.c[2]  * (1/surv.diff.c[1])

S.2   <- surv.diff.j[2]  * (1/surv.diff.j[1])

S.3 <- surv.diff.f[4] * (1/surv.diff.f[3])
S.4 <- surv.diff.f[5] * (1/surv.diff.f[4])
S.5 <- surv.diff.f[6] * (1/surv.diff.f[5])
S.6 <- surv.diff.f[7] * (1/surv.diff.f[6])
S.7 <- surv.diff.f[8] * (1/surv.diff.f[7])
S.8 <- surv.diff.f[9] * (1/surv.diff.f[8])

S.9 <- surv.diff.m[4] * (1/surv.diff.m[3])
S.10 <- surv.diff.m[5] * (1/surv.diff.m[4])
S.11 <- surv.diff.m[6] * (1/surv.diff.m[5])
S.12 <- surv.diff.m[7] * (1/surv.diff.m[6])
S.13 <- surv.diff.m[8] * (1/surv.diff.m[7])
S.14 <- surv.diff.m[9] * (1/surv.diff.m[8])

# Standard error of survival estimates

S.1.se <- surv.diff.c.se[2]

S.2.se <- surv.diff.j.se[2]

S.3.se <- surv.diff.f.se[4]
S.4.se <- surv.diff.f.se[5]
S.5.se <- surv.diff.f.se[6]
S.6.se <- surv.diff.f.se[7]
S.7.se <- surv.diff.f.se[8]
S.8.se <- surv.diff.f.se[9]

S.9.se <- surv.diff.m.se[4]
S.10.se <- surv.diff.m.se[5]
S.11.se <- surv.diff.m.se[6]
S.12.se <- surv.diff.m.se[7]
S.13.se <- surv.diff.m.se[8]
S.14.se <- surv.diff.m.se[9]

###########################################################################################
# df
###########################################################################################

age.class <- c("S.1",
               "S.2",
               "S.3", "S.4", "S.5", "S.6", "S.7", "S.8", "S.9", "S.10",
               "S.11", "S.12", "S.13", "S.14")
#sex <- c("mix", "mix", "f", "f", "f", "f", "f", "f", "m", "m", "m", "m", "m", "m")
#order <- c(1,2,1,2,3,4,5,6,1,2,3,4,5,6)
surv <- c(S.1,
          S.2,
          S.3, S.4, S.5, S.6, S.7, S.8, S.9, S.10,
          S.11, S.12, S.13, S.14)
surv.se <- c(S.1.se,
             S.2.se,
             S.3.se, S.4.se, S.5.se, S.6.se, S.7.se, S.8.se, S.9.se, S.10.se,
             S.11.se, S.12.se, S.13.se, S.14.se)

kzn.surv.df <- data.frame(age.class, surv, surv.se)

###########################################################################################
# SAVE
###########################################################################################

saver(surv.diff.c,
      surv.diff.j,
      surv.diff.f, surv.diff.m,
      surv.diff.c.se,
      surv.diff.j.se,
      surv.diff.f.se, surv.diff.m.se,
      kzn.surv.df,
      name = 'survival_estimates_KZN')

###########################################################################################
# END
###########################################################################################
