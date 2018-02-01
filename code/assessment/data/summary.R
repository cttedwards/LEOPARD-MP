
dat <- read.csv("summary.csv", header = TRUE, row.names = 1)

windows(); par(mfrow = c(4,2), mar = c(4,4,1,1))
plot(n_infanticide~n_adult, dat, type = 'b', pch = 19, cex = 2, xlim = c(40, 60))
plot(n_infanticide~adult_sex_ratio, dat, type = 'b', pch = 19, cex = 2, xlim = c(0, 2))

plot(harem_size~n_adult, dat, type = 'b', pch = 19, cex = 2, xlim = c(40, 60))
plot(harem_size~adult_sex_ratio, dat, type = 'b', pch = 19, cex = 2, xlim = c(0, 2))

plot(interbirth_interval~n_adult, dat, type = 'b', pch = 19, cex = 2, xlim = c(40, 60))
plot(interbirth_interval~adult_sex_ratio, dat, type = 'b', pch = 19, cex = 2, xlim = c(0, 2))

plot(n_male_dispersed~n_adult, dat, type = 'b', pch = 19, cex = 2, xlim = c(40, 60))
plot(n_male_dispersed~adult_sex_ratio, dat, type = 'b', pch = 19, cex = 2, xlim = c(0, 2))


