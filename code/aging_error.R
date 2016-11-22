
# proportion of age/sex classes misidentified (Balme et al. 2012)
aging.error <- c(id.f    = 0.32,    # proportion of females misidentified as males
                 id.m    = 0.31,    # proportion of males misidentified as females
                 age.m.1 = 0.61,    # male age class misidentified (age.m.1 = males at age 1)
                 age.m.2 = 0.61,    # male age class misidentified
                 age.m.3 = 0.60,    # male age class misidentified
                 age.m.4 = 0.47,    # male age class misidentified
                 age.m.5 = 0.47,    # male age class misidentified
                 age.m.6 = 0.47,    # male age class misidentified
                 age.m.7 = 0.46)    # male age class misidentified

# se for aging.error (Balme et al. 2012)
aging.error.se <- c(id.f    = 0.01,
                    id.m    = 0.01,
                    age.m.1 = 0.01,
                    age.m.2 = 0.01,
                    age.m.3 = 0.01,
                    age.m.4 = 0.01,
                    age.m.5 = 0.01,
                    age.m.6 = 0.01,
                    age.m.7 = 0.01)

saver(aging.error, aging.error.se, name = 'aging_error')


