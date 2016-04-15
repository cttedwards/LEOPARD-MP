##################################################################################################
# SET WORKING DIRECTORY AND LOAD UTILITY FUNCTIONS
##################################################################################################

rm(list=ls())

setwd("/Users/RossTyzackPitman/Documents/OneDrive/PhD/Data/R_Database/R_PROJECTS/PhD_Chapter3/MSE_Paper/LEOPARD-MP/code")
source('utils/saver.r')
source('utils/loader.r')
source('utils/reader.r')
source('utils/writer.r')

library(ggplot2)
library(reshape)
library(plyr)

loader('secr.results')

##################################################################################################
# Area
##################################################################################################

area <- c(400,400,400,   # welgevonden
          150,150,150,   # wonderkop
          230,230,230,   # atherstone
          320,320,       # vlnr
          310,310,       # lajuma
          540,540,540,   # timbavati
          230,230)       # makalali

##################################################################################################
# Yearly results
##################################################################################################

secr.results.2013 <- subset(secr.results, Year == "2013")
secr.results.2014 <- subset(secr.results, Year == "2014")
secr.results.2015 <- subset(secr.results, Year == "2015")

##################################################################################################
# Averaged secr results
##################################################################################################

ave.secr <- data.frame(year=c("2013", "2014", "2015"),
                          estimate=c(NA,NA,NA),
                          sd=c(NA,NA,NA))

ave.secr[1,2] <- mean(secr.results.2013$Density, na.rm=T)
ave.secr[1,3] <- sd(secr.results.2013$Density, na.rm=T)

ave.secr[2,2] <- mean(secr.results.2014$Density, na.rm=T)
ave.secr[2,3] <- sd(secr.results.2014$Density, na.rm=T)

ave.secr[3,2] <- mean(secr.results.2015$Density, na.rm=T)
ave.secr[3,3] <- sd(secr.results.2015$Density, na.rm=T)

ave.secr$lambda <- NA
ave.secr$lambda[1]  <- 1
ave.secr$lambda[2]  <- ave.secr$estimate[2] / ave.secr$estimate[1]
ave.secr$lambda[3]  <- ave.secr$estimate[3] / ave.secr$estimate[1]

##################################################################################################
# Plot
##################################################################################################

ave.secr.plot <- ggplot(data=ave.secr, aes(x = as.factor(year), y = estimate)) + 
  geom_line() + 
  geom_errorbar(width=.1, aes(ymin=estimate-sd, ymax=estimate+sd)) +
  geom_point(shape=21,size=3,fill="white") +
  theme_bw() + theme(strip.background = element_rect(fill="white")) +
  xlab("") + ylab("Density\n per 100km2") + ylim(0,10)

ave.secr.plot  

##################################################################################################
# End
##################################################################################################




# do not run

secr.results

# Welgevonden
welg.gr.2013_2014 <- (secr.results$Density[6] - secr.results$Density[5]) / secr.results$Density[5] * 100
welg.gr.2013_2015 <- (secr.results$Density[7] - secr.results$Density[5]) / secr.results$Density[5] * 100 / 2

# Atherstone
ath.gr.2013_2014 <- (secr.results$Density[12] - secr.results$Density[11]) / secr.results$Density[11] * 100
ath.gr.2013_2015 <- (secr.results$Density[13] - secr.results$Density[11]) / secr.results$Density[11] * 100 / 2

# Lajuma
lj.gr.2014_2015 <- (secr.results$Density[23] - secr.results$Density[22]) / secr.results$Density[22] * 100

# Makalali
mak.gr.2014_2015 <- (secr.results$Density[28] - secr.results$Density[27]) / secr.results$Density[27] * 100

# Timbavati
tim.gr.2013_2014 <- (secr.results$Density[25] - secr.results$Density[24]) / secr.results$Density[24] * 100

# Venetia
vlnr.gr.2014_2015 <- (secr.results$Density[15] - secr.results$Density[14]) / secr.results$Density[14] * 100

# Wonderkop
won.gr.2013_2014 <- (secr.results$Density[9] - secr.results$Density[8]) / secr.results$Density[8] * 100
won.gr.2013_2015 <- (secr.results$Density[10] - secr.results$Density[8]) / secr.results$Density[8] * 100 / 2

### growth rate
growth.rate <- data.frame(year=c("2013", "2014", "2015"),
                          lambda=c(1,NA,NA),
                          sd=c(0,NA,NA))

growth.rate[2,2] <- 1 + (mean(welg.gr.2013_2014, ath.gr.2013_2014, tim.gr.2013_2014, won.gr.2013_2014) / 100)
growth.rate[2,3] <- sd(c(-25.4902, -16.98113, -2.545455, -70.58824)) / 100

growth.rate[3,2] <- 1 + (mean(welg.gr.2013_2015, ath.gr.2013_2015, lj.gr.2014_2015, mak.gr.2014_2015, vlnr.gr.2014_2015, won.gr.2013_2015) / 100)
growth.rate[3,3] <- sd(c(-23.52941, -0.9433962, -34.31373, 49.20635, -9.615385, -7.843137)) / 100


secr.ave.Combined <- ggplot(data=ave.secr, aes(x = as.factor(year), y = estimate)) + geom_line() + 
  geom_errorbar(width=.1, aes(ymin=estimate-sd, ymax=estimate+sd)) +
  geom_point(shape=21,size=3,fill="white")

secr.ave.Combined = secr.Combined + theme_bw() + theme(strip.background = element_rect(fill="white"))
secr.ave.Combined = secr.Combined + xlab("") + ylab("Density\n per 100km2")
secr.ave.Combined = secr.Combined + theme(axis.title.x=element_text(size=20))
secr.ave.Combined = secr.Combined + theme(axis.title.y=element_text(size=20))
secr.ave.Combined = secr.Combined + theme(axis.text.x=element_text(size=15))
secr.ave.Combined = secr.Combined + theme(axis.text.y=element_text(size=15)) + ylim(0,15)
secr.ave.Combined
# write to file 

ppi = 150
tiff("secr.ave.combined.plot.tiff",width=8*ppi,height=5*ppi,res=ppi, compression = "lzw")
secr.Combined
dev.off()


