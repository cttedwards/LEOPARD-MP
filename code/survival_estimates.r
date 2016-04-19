
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

#sabi.surv.cub = reader("sabi.surv.cub")
#sabi.surv.juv = reader("sabi.surv.juv")
#sabi.surv.sa = reader("sabi.surv.sa")
#sabi.surv.a = reader("sabi.surv.a")
sabi.surv.cj = reader("sabi.surv.cub.juv")
sabi.surv = reader("sabi.surv")

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

sabi.surv.cj.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, 
                             data = sabi.surv.cj)                                                              # cj survival (from all mothers)
plot(sabi.surv.cj.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="All DEPENDENT SURVIVAL PROBABILITY", conf.int = TRUE)


sabi.surv.cj1.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, 
                              data = subset(sabi.surv.cj,MOTHER.AGE.CLASS.LITTER=="ADULT_36_48"))              # cj survival (from mother 36-48 months)
plot(sabi.surv.cj1.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="DEPENDENT1 SURVIVAL PROBABILITY", conf.int = TRUE)

sabi.surv.cj2.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, 
                              data = subset(sabi.surv.cj,MOTHER.AGE.CLASS.LITTER=="ADULT_48_60"))              # cj survival (from mother 48-60 months)
plot(sabi.surv.cj2.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="DEPENDENT2 SURVIVAL PROBABILITY", conf.int = TRUE)

sabi.surv.cj3.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, 
                              data = subset(sabi.surv.cj,MOTHER.AGE.CLASS.LITTER=="ADULT_60_72"))              # cj survival (from mother 60-72 months)
plot(sabi.surv.cj3.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="DEPENDENT3 SURVIVAL PROBABILITY", conf.int = TRUE)

sabi.surv.cj4.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, 
                              data = subset(sabi.surv.cj,MOTHER.AGE.CLASS.LITTER=="ADULT_72_84"))              # cj survival (from mother 72-84 months)
plot(sabi.surv.cj4.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="DEPENDENT4 SURVIVAL PROBABILITY", conf.int = TRUE)

sabi.surv.cj5.fit <- survfit(Surv(time=AGE.START.MONTHS,time2=AGE.MONTHS.24,event=AGE.MET.24,type='counting')~1, 
                              data = subset(sabi.surv.cj,MOTHER.AGE.CLASS.LITTER=="ADULT_84"))                 # cj survival (from mother >84 months)
plot(sabi.surv.cj5.fit, col=c(1:2), xlab="AGE (MONTHS)", ylab="DEPENDENT5 SURVIVAL PROBABILITY", conf.int = TRUE)

###########################################################################################
# ESTIMATE SURVIVAL PROBABILITIES
###########################################################################################

# dependent young

# cub (from all mothers)
surv.est.c <- stepfun(sabi.surv.cj.fit$time, c(1, sabi.surv.cj.fit$surv))
surv.diff.c <- surv.est.c(c(0,12))
surv.est.c.se <- stepfun(sabi.surv.cj.fit$time, c(1, sabi.surv.cj.fit$std.err))
surv.diff.c.se <- surv.est.c.se(c(0,12))

# cub from mother aged 36-48 months
surv.est.c1 <- stepfun(sabi.surv.cj1.fit$time, c(1, sabi.surv.cj1.fit$surv))
surv.diff.c1 <- surv.est.c1(c(0,12))
surv.est.c1.se <- stepfun(sabi.surv.cj1.fit$time, c(1, sabi.surv.cj1.fit$std.err))
surv.diff.c1.se <- surv.est.c1.se(c(0,12))

# cub from mother aged 48-60 months
surv.est.c2 <- stepfun(sabi.surv.cj2.fit$time, c(1, sabi.surv.cj2.fit$surv))
surv.diff.c2 <- surv.est.c2(c(0,12))
surv.est.c2.se <- stepfun(sabi.surv.cj2.fit$time, c(1, sabi.surv.cj2.fit$std.err))
surv.diff.c2.se <- surv.est.c2.se(c(0,12))

# cub from mother aged 60-72 months
surv.est.c3 <- stepfun(sabi.surv.cj3.fit$time, c(1, sabi.surv.cj3.fit$surv))
surv.diff.c3 <- surv.est.c3(c(0,12))
surv.est.c3.se <- stepfun(sabi.surv.cj3.fit$time, c(1, sabi.surv.cj3.fit$std.err))
surv.diff.c3.se <- surv.est.c3.se(c(0,12))

# cub from mother aged 72-84 months
surv.est.c4 <- stepfun(sabi.surv.cj4.fit$time, c(1, sabi.surv.cj4.fit$surv))
surv.diff.c4 <- surv.est.c4(c(0,12))
surv.est.c4.se <- stepfun(sabi.surv.cj4.fit$time, c(1, sabi.surv.cj4.fit$std.err))
surv.diff.c4.se <- surv.est.c4.se(c(0,12))

# cub from mother aged 84 months
surv.est.c5 <- stepfun(sabi.surv.cj5.fit$time, c(1, sabi.surv.cj5.fit$surv))
surv.diff.c5 <- surv.est.c5(c(0,12))
surv.est.c5.se <- stepfun(sabi.surv.cj5.fit$time, c(1, sabi.surv.cj5.fit$std.err))
surv.diff.c5.se <- surv.est.c5.se(c(0,12))

# juvenile (from all mothers)
surv.est.j <- stepfun(sabi.surv.cj.fit$time, c(1, sabi.surv.cj.fit$surv))
surv.diff.j <- surv.est.j(c(12,24))
surv.est.j.se <- stepfun(sabi.surv.cj.fit$time, c(1, sabi.surv.cj.fit$std.err))
surv.diff.j.se <- surv.est.j.se(c(12,24))

# juvenile from mother aged 36-48 months
surv.est.j1 <- stepfun(sabi.surv.cj1.fit$time, c(1, sabi.surv.cj1.fit$surv))
surv.diff.j1 <- surv.est.j1(c(12,24))
surv.est.j1.se <- stepfun(sabi.surv.cj1.fit$time, c(1, sabi.surv.cj1.fit$std.err))
surv.diff.j1.se <- surv.est.j1.se(c(12,24))

# juvenile from mother aged 48-60 months
surv.est.j2 <- stepfun(sabi.surv.cj2.fit$time, c(1, sabi.surv.cj2.fit$surv))
surv.diff.j2 <- surv.est.j2(c(12,24))
surv.est.j2.se <- stepfun(sabi.surv.cj2.fit$time, c(1, sabi.surv.cj2.fit$std.err))
surv.diff.j2.se <- surv.est.j2.se(c(12,24))

# juvenile from mother aged 60-72 months
surv.est.j3 <- stepfun(sabi.surv.cj3.fit$time, c(1, sabi.surv.cj3.fit$surv))
surv.diff.j3 <- surv.est.j3(c(12,24))
surv.est.j3.se <- stepfun(sabi.surv.cj3.fit$time, c(1, sabi.surv.cj3.fit$std.err))
surv.diff.j3.se <- surv.est.j3.se(c(12,24))

# juvenile from mother aged 72-84 months
surv.est.j4 <- stepfun(sabi.surv.cj4.fit$time, c(1, sabi.surv.cj4.fit$surv))
surv.diff.j4 <- surv.est.j4(c(12,24))
surv.est.j4.se <- stepfun(sabi.surv.cj4.fit$time, c(1, sabi.surv.cj4.fit$std.err))
surv.diff.j4.se <- surv.est.j4.se(c(12,24))

# juvenile from mother aged 84 months
surv.est.j5 <- stepfun(sabi.surv.cj5.fit$time, c(1, sabi.surv.cj5.fit$surv))
surv.diff.j5 <- surv.est.j5(c(12,24))
surv.est.j5.se <- stepfun(sabi.surv.cj5.fit$time, c(1, sabi.surv.cj5.fit$std.err))
surv.diff.j5.se <- surv.est.j5.se(c(12,24))

# sub-adult and adults

# females
surv.est.f <- stepfun(sabi.surv.f.fit$time, c(1, sabi.surv.f.fit$surv))
surv.diff.f <- surv.est.f(c(0,12,24,36,48,60,72,84,150))
surv.est.f.se <- stepfun(sabi.surv.f.fit$time, c(1, sabi.surv.f.fit$std.err))
surv.diff.f.se <- surv.est.f.se(c(0,12,24,36,48,60,72,84,150))

# males
surv.est.m <- stepfun(sabi.surv.m.fit$time, c(1, sabi.surv.m.fit$surv))
surv.diff.m <- surv.est.m(c(0,12,24,36,48,60,72,84,150))
surv.est.m.se <- stepfun(sabi.surv.m.fit$time, c(1, sabi.surv.m.fit$std.err))
surv.diff.m.se <- surv.est.m.se(c(0,12,24,36,48,60,72,84,150))

###########################################################################################
# FINAL SURVIVAL ESTIMATES
###########################################################################################

S.1   <- surv.diff.c[2]  * (1/surv.diff.c[1])
S.1.1 <- surv.diff.c1[2] * (1/surv.diff.c1[1])
S.1.2 <- surv.diff.c2[2] * (1/surv.diff.c2[1])
S.1.3 <- surv.diff.c3[2] * (1/surv.diff.c3[1])
S.1.4 <- surv.diff.c4[2] * (1/surv.diff.c4[1])
S.1.5 <- surv.diff.c5[2] * (1/surv.diff.c5[1])

S.2   <- surv.diff.j[2]  * (1/surv.diff.j[1])
S.2.1 <- surv.diff.j1[2] * (1/surv.diff.j1[1])
S.2.2 <- surv.diff.j2[2] * (1/surv.diff.j2[1])
S.2.3 <- surv.diff.j3[2] * (1/surv.diff.j3[1])
S.2.4 <- surv.diff.j4[2] * (1/surv.diff.j4[1])
S.2.5 <- surv.diff.j5[2] * (1/surv.diff.j5[1])

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
S.1.1.se <- surv.diff.c1.se[2]
S.1.2.se <- surv.diff.c2.se[2]
S.1.3.se <- surv.diff.c3.se[2]
S.1.4.se <- surv.diff.c4.se[2]
S.1.5.se <- surv.diff.c5.se[2]

S.2.se <- surv.diff.j.se[2]
S.2.1.se <- surv.diff.j1.se[2]
S.2.2.se <- surv.diff.j2.se[2]
S.2.3.se <- surv.diff.j3.se[2]
S.2.4.se <- surv.diff.j4.se[2]
S.2.5.se <- surv.diff.j5.se[2]

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

age.class <- c("S.1", "S.1.1", "S.1.2", "S.1.3", "S.1.4", "S.1.5",
               "S.2", "S.2.1", "S.2.2", "S.2.3", "S.2.4", "S.2.5",
               "S.3", "S.4", "S.5", "S.6", "S.7", "S.8", "S.9", "S.10",
               "S.11", "S.12", "S.13", "S.14")
#sex <- c("mix", "mix", "f", "f", "f", "f", "f", "f", "m", "m", "m", "m", "m", "m")
#order <- c(1,2,1,2,3,4,5,6,1,2,3,4,5,6)
surv <- c(S.1, S.1.1, S.1.2, S.1.3, S.1.4, S.1.5,
          S.2, S.2.1, S.2.2, S.2.3, S.2.4, S.2.5,
          S.3, S.4, S.5, S.6, S.7, S.8, S.9, S.10,
          S.11, S.12, S.13, S.14)
surv.se <- c(S.1.se, S.1.1.se, S.1.2.se, S.1.3.se, S.1.4.se, S.1.5.se,
             S.2.se, S.2.1.se, S.2.2.se, S.2.3.se, S.2.4.se, S.2.5.se,
             S.3.se, S.4.se, S.5.se, S.6.se, S.7.se, S.8.se, S.9.se, S.10.se,
             S.11.se, S.12.se, S.13.se, S.14.se)

sabi.sands.surv.df <- data.frame(age.class, surv, surv.se)

###########################################################################################
# SAVE
###########################################################################################

saver(surv.diff.c,surv.diff.c1,surv.diff.c2,surv.diff.c3,surv.diff.c4,surv.diff.c5,
      surv.diff.j,surv.diff.j1,surv.diff.j2,surv.diff.j3,surv.diff.j4,surv.diff.j5,
      surv.diff.f, surv.diff.m,
      surv.diff.c.se,surv.diff.c1.se,surv.diff.c2.se,surv.diff.c3.se,surv.diff.c4.se,surv.diff.c5.se,
      surv.diff.j.se,surv.diff.j1.se,surv.diff.j2.se,surv.diff.j3.se,surv.diff.j4.se,surv.diff.j5.se,
      surv.diff.f.se, surv.diff.m.se,
      sabi.sands.surv.df,
      name = 'survival_estimates')

###########################################################################################
# END
###########################################################################################
