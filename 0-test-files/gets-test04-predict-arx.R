##################################################
## Test file for gets package. First created
## 25 July 2019.
##
## 1 INITIATE
## 2 TESTS OF MEAN PREDICTIONS
## 3 TESTS OF VARIANCE PREDICTIONS
## 4 TESTS OF plot.options ARGUMENTS
## 5 FURTHER TESTS
##
##################################################

##################################################
## 1 INITIATE
##################################################

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd(choose.dir())

##remove everything in workspace:
rm(list=ls())

##load source code:
source("gets-base-source.R")
source("gets-isat-source.R")


##################################################
## 2 TESTS OF MEAN PREDICTIONS
##################################################

##generate some data:
##===================

##small dgp:
set.seed(123)
vY <- rnorm(20)

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##ar(0) model without constant:
##=============================

mymodel <- arx(vY)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat1 <- yhat2 <- yhat3 <- 0
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(0) model:
##============

mymodel <- arx(vY, mc=TRUE)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat1 <- yhat2 <- yhat3 <- coef(mymodel)[1]
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(1) model:
##============

mymodel <- arx(vY, mc=TRUE, ar=1)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat0 <- vY[length(vY)] #actual value at forecast origin
yhat1 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat0
yhat2 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat1
yhat3 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat2
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(1)-x model:
##==============

mX <- rnorm(length(vY))
mymodel <- arx(vY, mc=TRUE, ar=1, mxreg=mX)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3,
  newmxreg=matrix(1,3,1))

##correct predictions:
yhat0 <- vY[length(vY)] #actual value at forecast origin
yhat1 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat0 + coef(mymodel)[3]*1
yhat2 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat1 + coef(mymodel)[3]*1
yhat3 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat2 + coef(mymodel)[3]*1
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )

##EqWMA(2) model:
##===============

mymodel <- arx(vY, ewma=list(length=2))
regressorsMean(vY, ewma=list(length=2))
mean( c(vY[1],vY[2]) ) #obs no. 3
mean( c(vY[2],vY[3]) ) #obs no. 4
mean( c(vY[18],vY[19]) ) #obs no. 20
mean( c(tail(vY,n=2),tail(vY,n=3)) )

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat1 <- coef(mymodel)[1]*mean( c(vY[20],vY[19]) )
yhat2 <- coef(mymodel)[1]*mean(c(yhat1,vY[20]))
yhat3 <- coef(mymodel)[1]*mean(c(yhat2,yhat1))
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )


##################################################
## 3 TESTS OF VARIANCE PREDICTIONS
##################################################

##generate some data:
##===================

##small dgp:
set.seed(123)
vY <- rnorm(20)

##ar(0) model without constant:
##=============================

mymodel <- arx(vY)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(0) model:
##============

mymodel <- arx(vY, mc=TRUE)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(1) model:
##============

mymodel <- arx(vY, mc=TRUE, ar=1)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##arch(0) model with constant:
##============================

mymodel <- arx(vY, vc=TRUE)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )
all( round(functionVals, digits=10) == round(correctVals, digits=10) )

##arch(1) models:
##===============

mymodel <- arx(vY, arch=1)
predict(mymodel, spec="variance")

mymodel <- arx(vY, mc=TRUE, arch=1)
predict(mymodel, spec="variance")

mymodel <- arx(vY, mc=TRUE, ar=1, arch=1)
predict(mymodel, spec="variance")

##ar(1)-x model:
##==============

mX <- rnorm(length(vY))
mymodel <- arx(vY, mc=TRUE, ar=1, mxreg=mX)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##arch(1)-x model:
##================

mymodel <- arx(vY, arch=1, vxreg=mX)
predict(mymodel, spec="variance", n.ahead=3, newvxreg=matrix(1,3,1))


##################################################
## 4 TESTS OF plot.options ARGUMENTS
##################################################

##ar(1) model:
##============

mymodel <- arx(vY, mc=TRUE, ar=1)
predict(mymodel, plot.options=list(keep=1))
predict(mymodel, plot.options=list(fitted=TRUE))
predict(mymodel, plot.options=list(lty=c(3,2)))
predict(mymodel, plot.options=list(lwd=3))
predict(mymodel, plot.options=list(col=c("darkred","green")))
predict(mymodel, plot.options=list(shades.of.grey=c(80,60)))
predict(mymodel, plot.options=list(newmactual=rep(0,6)))
predict(mymodel, plot.options=list(ylim=c(-8,8)))

##arch(1) model:
##==============

mymodel <- arx(vY, vc=TRUE, arch=1)
predict(mymodel, plot.options=list(keep=1))
predict(mymodel, plot.options=list(fitted=TRUE))
predict(mymodel, plot.options=list(lty=c(3,2)))
predict(mymodel, plot.options=list(lwd=3))
predict(mymodel, plot.options=list(col=c("darkred","green")))
predict(mymodel, plot.options=list(shades.of.grey=c(60,40)))
predict(mymodel, plot.options=list(newvactual=rep(1,6)))
predict(mymodel, plot.options=list(ylim=c(-8,8)))


##################################################
## 5 FURTHER TESTS
##################################################

getOption("plot")
#options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##used to yield error:
set.seed(123)
y <- 2+rnorm(20)
mymodel <- arx(y, mc=TRUE)
predict(mymodel, n.ahead=1)

##used to yield error in the plot:
predict(mymodel, ci.levels=c(0.99,0.80, 0.50))

##used to yield error:
set.seed(123)
y <- 25+rnorm(100)
mymodel <- arx(y, mc=TRUE, ar=1:2) ##y has mean 25
pred <- predict(mymodel, n.ahead=20, plot.options=list(keep=50))
pred
##provide some actual values (needn't be same length as n.ahead)
y.actual <- 25+rnorm(12)
preda <- predict(mymodel, n.ahead=20,
  plot.options=list(newmactual=y.actual))
preda
