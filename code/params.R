
loader('survival_estimates')
loader('hunt.stats')

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

