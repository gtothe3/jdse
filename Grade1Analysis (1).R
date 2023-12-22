#########################################################################################
### R Code for JSDSE: Multiple Regression and Transformations at the Industrial-Scale ###
#########################################################################################

##################################
### Import and Clean Data File ###
##################################

rm(list = ls())

dataDir <- "/Users/ggreivel/Desktop/2023-12-19-Paper Revision Files/R-Code and Other Files/"

setwd(dataDir)

forceData <- read.csv('Gr01Data.csv', header = TRUE)

forceData <- na.omit(forceData)

forceData <- forceData[forceData$force > 0 & forceData$temperature > 0,]


#########################################################
### Variable Definitions with Conversions to SI Units ###
#########################################################

forceData$force <- 8896 * forceData$force
forceData$radius <- forceData$roller.diameter/2
forceData$Temp.K <- ((5/9)*(forceData$temperature - 32)) + 273.15
forceData$width <- 25.4 * forceData$avg.width
forceData$entry.gauge <- 25.4 * forceData$entry.gauge
forceData$exit.gauge <- 25.4 * forceData$exit.gauge
forceData$strain <- log (forceData$entry.gauge/forceData$exit.gauge)
forceData$Delta.h <- forceData$entry.gauge - forceData$exit.gauge
forceData$c.length <- sqrt(forceData$radius * (forceData$entry.gauge - forceData$exit.gauge))

forceData$Stand <- as.factor(forceData$stand.id)

forceData$Caster <- as.factor(forceData$caster.id)


###########################################################
### Initial Plots of Force vs. Untransformed Predictors ###
###########################################################

plot(forceData$force~forceData$width,forceData,main="Force vs. Width", xlab="Width [mm]",ylab="Force [N]")

plot(forceData$force~forceData$Temp.K,forceData,main="Force vs. Temperature", xlab="Temperature [K]",ylab="Force [N]")

plot(forceData$force~forceData$entry.gauge,forceData,main="Force vs. Entry Gauge", xlab="Entry Gauge [mm]",ylab="Force [N]")

plot(forceData$force~forceData$exit.gauge,forceData,main="Force vs. Exit Gauge", xlab="Exit Gauge [mm]",ylab="Force [N]")

plot(forceData$force~forceData$radius,forceData,main="Force vs. Roller Radius", xlab="Roller Radius [mm]",ylab="Force [N]")

plot(forceData$force~forceData$Stand,forceData,main="Force vs. Roll Stand", xlab="Roll Stand [-]",ylab="Force [N]")

plot(forceData$force~forceData$Caster,forceData,main="Force vs. Caster", xlab="Caster [-]",ylab="Force [N]")

#Pairs Plots...
df1 <- data.frame(forceData$force, forceData$width, forceData$Temp.K, forceData$entry.gauge, forceData$exit.gauge, forceData$radius)
pairs(df1)

##########################################
### Model N1 Construction and Analysis ###
##########################################

naivelm <- lm(forceData$force ~ forceData$width + forceData$Temp.K + forceData$entry.gauge + forceData$exit.gauge + forceData$radius, data = forceData)
summary(naivelm)
plot(naivelm)




prnaive <- resid(naivelm)/(1 - lm.influence(naivelm)$hat)
(pressnaive <- sum(prnaive^2))
sseN1 <- sum((fitted(naivelm) - forceData$force)^2)
ssrN1 <- sum((fitted(naivelm) - mean(forceData$force))^2)
(sstN1 <- ssrN1 + sseN1)
(R2PredNaive <- 1 - (pressnaive/sstN1))


################################################################################
### Log Transformed Variables and Plots of log(F) vs. Transformed Predictors ###
################################################################################

### Log transformed Regression Models Variables ###

forceData$log.force <- log(forceData$force)
forceData$X.T <- 1/(forceData$Temp.K)
forceData$X.epsilon <- log(forceData$strain)
forceData$X.w <- log(forceData$width)
forceData$X.r <- log(forceData$radius)
forceData$X.Delta.h <- log(forceData$Delta.h)
forceData$X.c.length <- log(forceData$c.length)

### Plots

plot(forceData$log.force ~ forceData$X.T,forceData,main="Log(Force) vs. 1/T", xlab="1/T",ylab="Log(Force)")

plot(forceData$log.force ~ forceData$X.epsilon,forceData,main="Log(Force) vs. Log(Strain)", xlab="Log(Strain)",ylab="Log(Force)")

plot(forceData$log.force ~ forceData$X.w,forceData,main="Log(Force) vs. Log(Width)", xlab="Log(Width)",ylab="Log(Force)")

plot(forceData$log.force ~ forceData$X.c.length,forceData,main="Log(Force) vs. Log(Contact Length)", xlab="Log(Contact Length)",ylab="Log(Force)")


##########################################
### Model T1 Construction and Analysis ###
##########################################

translm1 <- lm(forceData$log.force ~ forceData$X.T + forceData$X.epsilon + forceData$X.w + forceData$X.c.length, data = forceData)
summary(translm1)
plot(translm1)




prT1 <- resid(translm1)/(1 - lm.influence(translm1)$hat)
(pressT1 <- sum(prT1^2))
sseT1 <- sum((fitted(translm1) - forceData$log.force)^2)
ssrT1 <- sum((fitted(translm1) - mean(forceData$log.force))^2)
(sstT1 <- ssrT1 + sseT1)
(R2PredT1 <- 1 - (pressT1/sstT1))


########################################################
### Explore Collinearity (via Pairwise correlations) ###
########################################################

round(cor(forceData[c("X.T", "X.epsilon", "X.w", "X.c.length")]),3)
forceData.cor = cor(forceData[c("X.T", "X.epsilon", "X.w", "X.c.length")])

library(corrplot)
corrplot(forceData.cor)

# VIFs...
library(faraway)
x1 <- model.matrix(translm1) [,-1]
vif(x1)

##########################################
### Model T2 Construction and Analysis ###
##########################################

translm2 <- lm(forceData$log.force ~ forceData$X.T + forceData$X.epsilon + forceData$X.w, data = forceData)
summary(translm2)
#plot(translm2)

prT2 <- resid(translm2)/(1 - lm.influence(translm2)$hat)
(pressT2 <- sum(prT2^2))
sseT2 <- sum((fitted(translm2) - forceData$log.force)^2)
ssrT2 <- sum((fitted(translm2) - mean(forceData$log.force))^2)
(sstT2 <- ssrT2 + sseT2)
(R2PredT2 <- 1 - (pressT2/sstT2))

x2 <- model.matrix(translm2) [,-1]
vif(x2)

##########################################################################################
### Plots of 1/Temp vs. Stand & Color Plots of F varaibles vs. Temp variables by Stand ### 
##########################################################################################

plot(forceData$Temp.K~forceData$Stand,forceData,main="Temperature vs. Roll Stand", xlab="Roll Stand [-]",ylab="Temperature [K]")

plot(forceData$force ~ forceData$Temp.K, col=c("red","orange","yellow","green4","blue")[forceData$Stand], xlim = c(1050,1350), ylim = c(4500000,37500000), xaxp = c(1100,1300,4), yaxp = c(5000000,36000000,6), main="Force vs. Temperature (All Stands, R^2 = 0.745)", xlab="Temperature [K]", ylab="Force [N]")
summary(lm(forceData$force ~ forceData$Temp.K))

plot(forceData$log.force ~ forceData$X.T, col=c("red","orange","yellow","green4","blue")[forceData$Stand], xlim = c(0.00075,0.00095), ylim = c(15,18), xaxp = c(0.00075,0.00095,4), yaxp = c(15,18,6), main="Log(Force) vs. 1/Temperature (All Stands, R^2 = 0.789)", xlab="1/Temperature", ylab="Log(Force)")
summary(lm(forceData$log.force ~ forceData$X.T))


###########################################
#### Model T2b Construction and Analysis ###
###########################################
#
#translm2b <- lm(forceData$log.force ~ forceData$X.epsilon + forceData$X.w, data = forceData)
#summary(translm2b)
#plot(translm2b)
#
#prT2b <- resid(translm2b)/(1 - lm.influence(translm2b)$hat)
#(pressT2b <- sum(prT2^2))
#sseT2b <- sum((fitted(translm2b) - forceData$log.force)^2)
#ssrT2b <- sum((fitted(translm2b) - mean(forceData$log.force))^2)
#(sstT2b <- ssrT2b + sseT2b)
#(R2PredT2b <- 1 - (pressT2b/sstT2b))
#
#x2b <- model.matrix(translm2b) [,-1]
#vif(x2b)
#

##########################################
### Model T3 Construction and Analysis ###
##########################################

translm3 <- lm(forceData$log.force ~ forceData$X.T + forceData$X.epsilon + forceData$X.w + forceData$Stand, data = forceData)
summary(translm3)
#plot(translm3)

prT3 <- resid(translm3)/(1 - lm.influence(translm3)$hat)
(pressT3 <- sum(prT3^2))
sseT3 <- sum((fitted(translm3) - forceData$log.force)^2)
ssrT3 <- sum((fitted(translm3) - mean(forceData$log.force))^2)
(sstT3 <- ssrT3 + sseT3)
(R2PredT3 <- 1 - (pressT3/sstT3))

x3 <- model.matrix(translm3) [,-1]
vif(x3)
# Very Large VIFs

#####################################################################
### Model Selection Analysis with C_p & R^2_adj Based on Model T3 ###
#####################################################################

require(leaps)
selection <- regsubsets(forceData$log.force ~ forceData$X.T + forceData$X.epsilon + forceData$X.w + forceData$Stand, data = forceData)
rselection <- summary(selection)
rselection$which

plot(1:7,rselection$adjr2, main = expression(paste("Adjusted R"^"2"*" vs. Number of Variables")), xlab = "Number of Variables", ylab = expression(paste("Adjusted R"^"2")))
which.max(rselection$adjr2)

plot(1:7,rselection$cp, main = expression(paste("Mallows C"[p]*" vs. Number of Variables")), xlab = "Number of Variables", ylab = expression(paste("Mallows C"[p]*" Statistic")))
abline(0,1)

plot(1:7,rselection$bic, main = "BIC vs. Number of Variables", xlab = "Number of Variables", ylab = "BIC")
which.min(rselection$bic)


##########################################
### Model T4 Construction and Analysis ###
##########################################

translm4 <- lm(forceData$log.force ~ forceData$X.epsilon + forceData$X.w + forceData$Stand, data = forceData)
summary(translm4)
plot(translm4)




#library(car)
#qqPlot(translm4$residuals, main = "Normal Q-Q Plot of Residuals",xlab = "Theoretical Quantiles", ylab = "Residuals")

prT4 <- resid(translm4)/(1 - lm.influence(translm4)$hat)
(pressT4 <- sum(prT4^2))
sseT4 <- sum((fitted(translm4) - forceData$log.force)^2)
ssrT4 <- sum((fitted(translm4) - mean(forceData$log.force))^2)
(sstT4 <- ssrT4 + sseT4)
sstT4
(R2PredT4 <- 1 - (pressT4/sstT4))

x4 <- model.matrix(translm4) [,-1]
vif(x4)


### User-Defined Residual Plots for Model T4

qqnorm(residuals(translm4b), main = "Normal Q-Q Plot of Residuals", xlab = "Theoretical Quantiles", ylab = "Residuals")
qqline(residuals(translm4b))

hist(residuals(translm4b), main = "Histogram of Residuals", xlab = "Residuals", ylab = "Frequency")

plot(fitted(translm4b), residuals(translm4b), main = "Residuals vs. Fitted Values", xlab = "Fitted", ylab = "Residuals")
abline(h=0)

