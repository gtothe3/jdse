####### Correct Initial Plots Code in R #########

rm(list = ls())

dataDir <- "/Users/ggreivel/Desktop/ORwE/Nucor/Paper 3 Files/"

setwd(dataDir)
 
forceData <- read.csv('Gr01_Data.csv', header = TRUE)
 
forceData <- na.omit(forceData)
 
forceData <- forceData[forceData$force > 0 & forceData$temperature > 0,]
 
forceData$Stand <- as.factor(forceData$stand.id)

forceData$force <- 8896 * forceData$force
forceData$radius <- forceData$roller.diameter/2
forceData$Temp.K <- ((5/9)*(forceData$temperature - 32)) + 273.15
forceData$width <- 25.4 * forceData$avg.width
forceData$entry.gauge <- 25.4 * forceData$entry.gauge
forceData$exit.gauge <- 25.4 * forceData$exit.gauge
forceData$strain <- log (forceData$entry.gauge/forceData$exit.gauge)
forceData$Delta.h <- forceData$entry.gauge - forceData$exit.gauge
forceData$c.length <- sqrt(forceData$radius * (forceData$entry.gauge - forceData$exit.gauge))
 

#############################################
### Variables for Strain Rate Model Terms ###
#############################################

forceData$velocity <- (1095.96)*(28.295/forceData$exit.gauge)
forceData$Delta.t <- (forceData$c.length)/(forceData$velocity)
forceData$strain.rate <- (forceData$strain)/(forceData$Delta.t)

##################################
### Regression Model Variables ###
##################################

forceData$log.force <- log(forceData$force)
forceData$X.T <- 1/(forceData$Temp.K)
forceData$X.v <- log(forceData$velocity)
forceData$X.epsilon <- log(forceData$strain)
forceData$X.w <- log(forceData$width)
forceData$X.r <- log(forceData$radius)
forceData$X.Delta.h <- log(forceData$Delta.h)
forceData$X.c.length <- log(forceData$c.length)

#########################

plot(forceData$log.force~forceData$X.T,forceData,main="log(Force) vs. 1/T", xlab="1/T",ylab="log(Force)")

plot(forceData$force~forceData$Temp.K,forceData,main="Force vs. Temperature", xlab="Temperature [K]",ylab="Force [N]")

plot(forceData$force~forceData$entry.gauge,forceData,main="Force vs. Entry Gauge", xlab="Entry Gauge [mm]",ylab="Force [N]")

plot(forceData$force~forceData$exit.gauge,forceData,main="Force vs. Exit Gauge", xlab="Exit Gauge [mm]",ylab="Force [N]")

plot(forceData$force~forceData$radius,forceData,main="Force vs. Roller Radius", xlab="Roller Radius [mm]",ylab="Force [N]")

plot(forceData$force~forceData$Stand,forceData,main="Force vs. Roll Stand", xlab="Roll Stand [-]",ylab="Force [N]")

plot(forceData$Temp.K~forceData$Stand,forceData,main="Temperature vs. Roll Stand", xlab="Roll Stand [-]",ylab="Temperature [K]")


##############################

naivelm <- lm(forceData$force ~ forceData$width + forceData$Temp.K + forceData$entry.gauge + forceData$exit.gauge + forceData$radius, data = forceData)
summary(naivelm)
plot(naivelm)

###########################

plot(forceData$log.force ~ forceData$X.T,forceData,main="Log(Force) vs. 1/T", xlab="1/T",ylab="Log(Force)")

plot(forceData$log.force ~ forceData$X.epsilon,forceData,main="Log(Force) vs. Log(Strain)", xlab="Log(Strain)",ylab="Log(Force)")

plot(forceData$log.force ~ forceData$X.w,forceData,main="Log(Force) vs. Log(Width)", xlab="Log(Width)",ylab="Log(Force)")

plot(forceData$log.force ~ forceData$X.c.length,forceData,main="Log(Force) vs. Log(Contact Length)", xlab="Log(Contact Length)",ylab="Log(Force)")

##############################

translm1 <- lm(forceData$log.force ~ forceData$X.T + forceData$X.epsilon + forceData$X.w + forceData$X.c.length, data = forceData)
summary(translm1)
plot(translm1)


########################### Look for Pairwise correlations/collinearities

round(cor(forceData$X.epsilon,forceData$X.c.length),3)
cor(forceData$X.epsilon,forceData$X.c.length)

round(cor(forceData[c("X.T", "X.epsilon", "X.w", "X.c.length")]),3)
cor(forceData[c("X.T", "X.epsilon", "X.w", "X.c.length")])
forceData.cor = cor(forceData[c("X.T", "X.epsilon", "X.w", "X.c.length")])

library(corrplot)
corrplot(forceData.cor)


########################## Model without Contact Length Term...

translm2 <- lm(forceData$log.force ~ forceData$X.T + forceData$X.epsilon + forceData$X.w, data = forceData)
summary(translm2)
plot(translm2)

################Color Plots of F vs. Temp by Stand ### In Paper3ColorPlot file

plot(forceData$force ~ forceData$Temp.K, col=c("red","orange","yellow","green4","blue")[forceData$Stand], xlim = c(1050,1350), ylim = c(4500000,37500000), xaxp = c(1100,1300,4), yaxp = c(5000000,36000000,6), main="Force vs. Temperature (All Stands, R^2 = 0.745)", xlab="Temperature [K]", ylab="Force [N]")
summary(lm(forceData$force ~ forceData$Temp.K))


plot(forceData$log.force ~ forceData$X.T, col=c("red","orange","yellow","green4","blue")[forceData$Stand], xlim = c(0.00075,0.00095), ylim = c(15,18), xaxp = c(0.00075,0.00095,4), yaxp = c(15,18,6), main="Log(Force) vs. 1/Temperature (All Stands, R^2 = 0.789)", xlab="1/Temperature", ylab="Log(Force)")
summary(lm(forceData$log.force ~ forceData$X.T))


############### Recoding roll stand to numerical values

forceDataS <- read.csv('Gr01_DataForStand.csv', header = TRUE)

forceDataS <- na.omit(forceDataS)

forceDataS <- forceDataS[forceDataS$force > 0 & forceDataS$temperature > 0,]

forceDataS$Stand <- as.factor(forceDataS$stand.id)

forceDataS$force <- 8896 * forceDataS$force
forceDataS$radius <- forceDataS$roller.diameter/2
forceDataS$Temp.K <- ((5/9)*(forceDataS$temperature - 32)) + 273.15
forceDataS$width <- 25.4 * forceDataS$avg.width
forceDataS$entry.gauge <- 25.4 * forceDataS$entry.gauge
forceDataS$exit.gauge <- 25.4 * forceDataS$exit.gauge
forceDataS$strain <- log (forceDataS$entry.gauge/forceDataS$exit.gauge)
forceDataS$Delta.h <- forceDataS$entry.gauge - forceDataS$exit.gauge
forceDataS$c.length <- sqrt(forceDataS$radius * (forceDataS$entry.gauge - forceDataS$exit.gauge))


#############################################
### Variables for Strain Rate Model Terms ###
#############################################

forceDataS$velocity <- (1095.96)*(28.295/forceDataS$exit.gauge)
forceDataS$Delta.t <- (forceDataS$c.length)/(forceDataS$velocity)
forceDataS$strain.rate <- (forceDataS$strain)/(forceDataS$Delta.t)

##################################
### Regression Model Variables ###
##################################

forceDataS$log.force <- log(forceDataS$force)
forceDataS$X.T <- 1/(forceDataS$Temp.K)
forceDataS$X.v <- log(forceDataS$velocity)
forceDataS$X.epsilon <- log(forceDataS$strain)
forceDataS$X.w <- log(forceDataS$width)
forceDataS$X.r <- log(forceDataS$radius)
forceDataS$X.Delta.h <- log(forceDataS$Delta.h)
forceDataS$X.c.length <- log(forceDataS$c.length)

##################

plot(forceDataS$log.force ~ forceDataS$X.S,forceData,main="Log(Force) vs. Stand", xlab="Stand",ylab="Log(Force)")



translm2 <- lm(forceDataS$log.force ~ forceDataS$X.T + forceDataS$X.epsilon + forceDataS$X.w, data = forceDataS)
summary(translm2)
plot(translm2)

translm3a <- lm(forceDataS$log.force ~ forceDataS$X.T + forceDataS$X.epsilon + forceDataS$X.w + forceDataS$X.S, data = forceDataS)
summary(translm3a)
plot(translm3a)

translm3b <- lm(forceDataS$log.force ~ forceDataS$X.T + forceDataS$X.epsilon + forceDataS$X.w + forceDataS$Stand, data = forceDataS)
summary(translm3b)
plot(translm3b)

##### Do best subsets analysis with AIC, C_p,  R^2_adj with the variables in model T3b

require(leaps)
selection <- regsubsets(forceDataS$log.force ~ forceDataS$X.T + forceDataS$X.epsilon + forceDataS$X.w + forceDataS$Stand, data = forceDataS)
rselection$which


#plot(1:7,rselection$adjr2, main = "Adjusted R^2 vs. Number of Parameters", xlab = "Parameters", ylab = "Adj. R^2")
plot(1:7,rselection$adjr2, main = expression(paste("Adjusted R"^"2"*" vs. Number of Variables")), xlab = "Number of Variables", ylab = expression(paste("Adjusted R"^"2")))
which.max(rselection$adjr2)

#plot(1:7,rselection$cp, "Mallows C_p vs. Number of Variables", xlab = "Variables", ylab = "C_p Statistic")
plot(1:7,rselection$cp, main = expression(paste("Mallows C"[p]*" vs. Number of Variables")), xlab = "Number of Variables", ylab = expression(paste("Mallows C"[p]*" Statistic")))
abline(0,1)


############

translm4a <- lm(forceDataS$log.force ~ forceDataS$X.epsilon + forceDataS$X.w + forceDataS$X.S, data = forceDataS)
summary(translm4a)
plot(translm4a)

translm4b <- lm(forceDataS$log.force ~ forceDataS$X.epsilon + forceDataS$X.w + forceDataS$Stand, data = forceDataS)
summary(translm4b)
plot(translm4b)

qqnorm(residuals(translm4b), main = "Normal Q-Q Plot of Residuals", xlab = "Theoretical Quantiles", ylab = "Residuals")
qqline(residuals(translm4b))

hist(residuals(translm4b), main = "Histogram of Residuals", xlab = "Residuals", ylab = "Frequency")

plot(fitted(translm4b), residuals(translm4b), main = "Residuals vs. Fitted Values", xlab = "Fitted", ylab = "Residuals")
abline(h=0)

translm5b <- lm(forceDataS$log.force ~ forceDataS$X.epsilon + forceDataS$Stand, data = forceDataS)
summary(translm5b)
plot(translm5b)



########################### Look for Pairwise correlations/collinearities

round(cor(forceData$X.epsilon,forceData$X.c.length),3)
cor(forceData$X.epsilon,forceData$X.c.length)

round(cor(forceDataS[c("X.T", "X.epsilon", "X.w", "X.S")]),3)
cor(forceDataS[c("X.T", "X.epsilon", "X.w", "X.S")])
forceDataS.cor = cor(forceDataS[c("X.T", "X.epsilon", "X.w", "X.S")])

library(corrplot)
corrplot(forceDataS.cor)



############### RMSE Calculations
#install.packages("qpcR")
#install.packages("mvp")
#library(qpcR)
#library(mvp)

#press(translm3c)

prnaive <- resid(naivelm)/(1 - lm.influence(naivelm)$hat)
(pressnaive <- sum(prnaive^2))
log(pressnaive)
sseN1 <- sum((fitted(naivelm) - forceData$force)^2)
ssrN1 <- sum((fitted(naivelm) - mean(forceData$force))^2)
sstN1 <- ssrN1 + sseN1
sstN1

pr1 <- resid(translm1)/(1 - lm.influence(translm1)$hat)
(press1 <- sum(pr1^2))
sseT1 <- sum((fitted(translm1) - forceData$log.force)^2)
ssrT1 <- sum((fitted(translm1) - mean(forceData$log.force))^2)
sstT1 <- ssrT1 + sseT1
sstT1

pr2 <- resid(translm2)/(1 - lm.influence(translm2)$hat)
(press2 <- sum(pr2^2))
sseT2 <- sum((fitted(translm2) - forceData$log.force)^2)
ssrT2 <- sum((fitted(translm2) - mean(forceData$log.force))^2)
sstT2 <- ssrT2 + sseT2
sstT2

pr3a <- resid(translm3a)/(1 - lm.influence(translm3a)$hat)
(press3a <- sum(pr3a^2))

pr3b <- resid(translm3b)/(1 - lm.influence(translm3b)$hat)
(press3b <- sum(pr3b^2))
sseT3 <- sum((fitted(translm3b) - forceDataS$log.force)^2)
ssrT3 <- sum((fitted(translm3b) - mean(forceDataS$log.force))^2)
sstT3 <- ssrT3 + sseT3
sstT3

pr4b <- resid(translm4b)/(1 - lm.influence(translm4b)$hat)
(press4b <- sum(pr4b^2))
sseT4 <- sum((fitted(translm4b) - forceDataS$log.force)^2)
ssrT4 <- sum((fitted(translm4b) - mean(forceDataS$log.force))^2)
sstT4 <- ssrT4 + sseT4
sstT4

