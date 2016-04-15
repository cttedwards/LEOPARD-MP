##################################################################################################
# SET WORKING DIRECTORY AND LOAD UTILITY FUNCTIONS
##################################################################################################

rm(list=ls())

setwd("/Users/RossTyzackPitman/Documents/OneDrive/PhD/Data/R_Database/R_PROJECTS/PhD_Chapter3/MSE_Paper/LEOPARD-MP/code")
source('utils/saver.r')
source('utils/loader.r')
source('utils/reader.r')
source('utils/writer.r')
source('utils/pdfr.r')

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
                          se=c(NA,NA,NA))

ave.secr[1,2] <- mean(secr.results.2013$Density, na.rm=T)
ave.secr[1,3] <- sd(secr.results.2013$Density, na.rm=T) / sqrt(4)

ave.secr[2,2] <- mean(secr.results.2014$Density, na.rm=T)
ave.secr[2,3] <- sd(secr.results.2014$Density, na.rm=T) / sqrt(7)

ave.secr[3,2] <- mean(secr.results.2015$Density, na.rm=T)
ave.secr[3,3] <- sd(secr.results.2015$Density, na.rm=T) / sqrt(6)

ave.secr$lambda <- NA
ave.secr$lambda[1]  <- 1
ave.secr$lambda[2]  <- ave.secr$estimate[2] / ave.secr$estimate[1]
ave.secr$lambda[3]  <- ave.secr$estimate[3] / ave.secr$estimate[1]

##################################################################################################
# Plot
##################################################################################################

ave.secr.plot <- ggplot(data=ave.secr, aes(x = as.factor(year), y = estimate)) + 
  geom_line() + 
  geom_errorbar(width=.1, aes(ymin=estimate-se, ymax=estimate+se)) +
  geom_point(shape=21,size=3,fill="white") +
  theme_bw() + theme(strip.background = element_rect(fill="white")) +
  xlab("year") + ylab("Density\n per 100km2") + ylim(0,10)

ave.secr.plot 
pdfr(ave.secr.plot, width = 5, height = 5, name = 'ave.secr.plot')

lambda.plot <- ggplot(data=ave.secr, aes(x = as.factor(year), y = lambda)) + 
  geom_line() +
  geom_point(shape=21,size=3,fill="white") +
  theme_bw() + theme(strip.background = element_rect(fill="white")) +
  xlab("year") + ylab("lambda") + ylim(0.5,1.25) 

lambda.plot
pdfr(lambda.plot, width = 5, height = 5, name = 'lambda.plot')

##################################################################################################
# End
##################################################################################################
