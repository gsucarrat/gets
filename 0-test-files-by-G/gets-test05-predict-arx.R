##################################################
## Test file for the predict.arx function. First
## created 25 July 2019.
##
## 1 INITIATE
## 2 TESTS OF MEAN PREDICTIONS
## 3 TESTS OF VARIANCE PREDICTIONS
## 4 TESTS OF plot.options ARGUMENTS
## 5 FURTHER TESTS
##
##################################################

##idea: new argument 'shades.of=c("grey","blue","red","green")'
##this would be in addition to the already existing argument 'shades'
##
##some useful code?:
##
##fc <- colorRampPalette(c("darkgreen", "green"))
##plot(rep(1, 255),col = fc(255), pch = 19, cex = 3)
##
##fc <- colorRampPalette(c("darkgreen", "white"))
##plot(rep(1, 255),col = fc(255), pch = 19, cex = 3)
##
##fc <- colorRampPalette(c("green", "white"))
##plot(rep(1, 255),col = fc(255), pch = 19, cex = 3)
##
##fc <- colorRampPalette(c("red", "white"))
##plot(rep(1, 255),col = fc(255), pch = 19, cex = 3)
##
##fc <- colorRampPalette(c("blue", "white"))
##plot(rep(1, 255),col = fc(255), pch = 19, cex = 3)
##
##fc <- colorRampPalette(c("darkblue", "white"))
##plot(rep(1, 255),col = fc(255), pch = 19, cex = 3)
##
##find colours:
##cols <- colors()
##cols <- cols[grep("blue", cols)]
##plot(rep(1, length(cols)),col = cols, pch = 19, cex = 3)


##################################################
##1 INITIATE
##################################################

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd(choose.dir())

##load required packages:
require(parallel)
require(zoo)

##remove everything in workspace (.GlobaleEnv:
rm(list=ls())

##load source:
source("./gets/R/gets-base-source.R")
source("./gets/R/gets-isat-source.R") #needed for sim()


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

##ar(0) model w/constant:
##=======================

mymodel <- arx(vY)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat1 <- yhat2 <- yhat3 <- coef(mymodel)[1]
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(0) model without constant:
##=============================

mymodel <- arx(vY, mc=FALSE)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat1 <- yhat2 <- yhat3 <- 0
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(1) model w/constant:
##=======================

mymodel <- arx(vY, ar=1)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=12)

##correct predictions:
yhat <- rep(NA,13)
yhat[1] <- vY[length(vY)] #actual value at forecast origin
for(i in 2:13){
  yhat[i] <- coef(mymodel)[1] + coef(mymodel)[2]*yhat[i-1]
}
correctVals <- yhat[-1]

##do they correspond?:
all( functionVals == correctVals )

##ar(1) model without constant:
##=============================

mymodel <- arx(vY, mc=FALSE, ar=1)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=12)

##correct predictions:
yhat <- rep(NA,13)
yhat[1] <- vY[length(vY)] #actual value at forecast origin
for(i in 2:13){
  yhat[i] <- coef(mymodel)[1]*yhat[i-1]
}
correctVals <- yhat[-1]

##do they correspond?:
all( functionVals == correctVals )

##ar(1)-x model w/constant:
##=========================

mX <- rnorm(length(vY))
mymodel <- arx(vY, ar=1, mxreg=mX)

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

##EqWMA(2) model without constant:
##================================

mymodel <- arx(vY, mc=FALSE, ewma=list(length=2))
regressorsMean(vY, mc=FALSE, ewma=list(length=2))
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

##yielded NA-predictions in version 0.27:
##=======================================

##some data
set.seed(123)
y <- rnorm(50)

##predictions by predict.arx():
mymodel <- arx(y, mc=FALSE, ar=1)
functionVals <- predict(mymodel, n.ahead=3)
all( !is.na(functionVals) ) #TRUE since version 0.28

##resolution of issue posted by mdwy62 on Github 18/9-2021
##https://github.com/gsucarrat/gets/issues/50
##========================================================

##the problem: T was not recognised as TRUE
##the following was erroneous before version 0.30:

##ar(0) model w/constant:
mymodel <- arx(vY, mc=T)

##predictions of the mean:
functionVals <- predict(mymodel, spec="mean", n.ahead=3)

##correct predictions:
yhat1 <- yhat2 <- yhat3 <- coef(mymodel)[1]
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals ) #TRUE since version 0.30


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

mymodel <- arx(vY, mc=FALSE)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(0) model w/constant:
##=======================

mymodel <- arx(vY)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)
##note: is it correct that the message "Set 'vc = TRUE'..."
##should be returned?

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##ar(1) model w/constant:
##=======================

mymodel <- arx(vY, ar=1)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##arch(0) model w/constant:
##=========================

mymodel <- arx(vY, vc=TRUE)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- exp(coef(mymodel, spec="v")[1])
#OLD:
#sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )
all( round(functionVals, digits=10) == round(correctVals, digits=10) )

##arch(1) models:
##===============

##just to check that things seem to work

mymodel <- arx(vY, mc=FALSE, arch=1)
predict(mymodel, spec="variance")

mymodel <- arx(vY, arch=1)
predict(mymodel, spec="variance")

mymodel <- arx(vY, ar=1, arch=1)
predict(mymodel, spec="variance")

##ar(1)-x model:
##==============

mX <- rnorm(length(vY))
mymodel <- arx(vY, ar=1, mxreg=mX)

##predictions of the variance:
functionVals <- predict(mymodel, spec="variance", n.ahead=3)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )

##arch(1)-x model:
##================

##just to check that things seem to work

mymodel <- arx(vY, arch=1, vxreg=mX)
predict(mymodel, spec="variance", n.ahead=3, newvxreg=matrix(1,3,1))


##################################################
## 4 TESTS OF plot.options ARGUMENTS
##################################################

##ar(1) model w/constant:
##=======================

mymodel <- arx(vY, ar=1)

##for visual inspection:
predict(mymodel) #benchmark
predict(mymodel, plot.options=list(keep=1))
predict(mymodel, plot.options=list(line.at.origin=TRUE))
predict(mymodel, plot.options=list(start.at.origin=FALSE))
predict(mymodel,
  plot.options=list(start.at.origin=FALSE, fitted=TRUE))
predict(mymodel, plot.options=list(dot.at.origin=FALSE))
predict(mymodel, plot.options=list(hlines=c(-2,-1,0,1,2)))
predict(mymodel, plot.options=list(col=c("darkred","green")))
predict(mymodel, plot.options=list(lty=c(3,2)))
##does not work, but should it?:
predict(mymodel, plot.options=list(lty=3))
##only the forecast is lwd=3, should both be?:
predict(mymodel, plot.options=list(lwd=3))
##now both are lwd=3:
predict(mymodel, plot.options=list(lwd=c(3,3)))
##now only actual is lwd=3:
predict(mymodel, plot.options=list(lwd=c(1,3)))
predict(mymodel, plot.options=list(ylim=c(-8,8)))
predict(mymodel, plot.options=list(ylab="G-values"))
predict(mymodel,
  plot.options=list(main="Plot is slightly lower when 'main' is specified"))
predict(mymodel,
  plot.options=list(legend.text=c("Prognose","Faktisk")))
predict(mymodel, plot.options=list(fitted=TRUE))
predict(mymodel, plot.options=list(newmactual=rep(0,6)))
predict(mymodel, plot.options=list(shades=c(95,50)))
predict(mymodel, plot.options=list(shades=c(50,95))) #invert shades
predict(mymodel, plot.options=list(shades=c(95,50)))

##arch(1) model:
##==============

mymodel <- arx(vY, mc=FALSE, vc=TRUE, arch=1)

##for visual inspection:
predict(mymodel) #benchmark
predict(mymodel, plot.options=list(keep=1))
predict(mymodel, plot.options=list(line.at.origin=TRUE))
predict(mymodel, plot.options=list(start.at.origin=FALSE))
predict(mymodel,
  plot.options=list(start.at.origin=FALSE, fitted=TRUE))
predict(mymodel, plot.options=list(dot.at.origin=FALSE))
predict(mymodel, plot.options=list(hlines=0:4))
predict(mymodel, plot.options=list(col=c("darkred","green")))
predict(mymodel, plot.options=list(lty=c(3,2)))
predict(mymodel, plot.options=list(lwd=3))
predict(mymodel, plot.options=list(ylim=c(-6,8)))
predict(mymodel, plot.options=list(ylab="G-values"))
predict(mymodel,
  plot.options=list(main="Plot is slightly lower when 'main' is specified"))
predict(mymodel,
  plot.options=list(legend.text=c("Prognose","Residualene kvadrert")))
predict(mymodel, plot.options=list(fitted=TRUE))
predict(mymodel, plot.options=list(newvactual=rep(1,6)))
predict(mymodel, plot.options=list(shades=c(95,50)))
predict(mymodel, plot.options=list(shades=c(50,95))) #invert shades


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

##used to yield error in the plotting:
set.seed(123); dgp.n <- 50
y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
y[1:10] <- y[1:10] + 4 #step-shift
arxmod <- arx(y, mc=TRUE, ar=1:2,
  mxreg=cbind(mX, tim(y, which.ones=c(7,16)))) 
predict(arxmod, n.ahead=1, newmxreg=matrix(0,1,5),
  plot.options=list(start.at.origin=FALSE,
  line.at.origin=TRUE, fitted=TRUE))

##used to produce graphical error in the plot; since version 0.25
##the following message is returned to the user: "'newindex' not
##entirely out-of-sample, so no plot produced"
set.seed(123)
y <- rnorm(20)
arxmod <- arx(y, mc=TRUE, ar=1)
##generates predictions with user-specified index, but no plot:
predict(arxmod, newindex=19:30)


##################################################
## 6 SIMULATIONS VS. ANALYTICAL FORMULA
##################################################

##TO DO: COMPARE THE SIMULATED QUANTILES AND THE 
##ANALYTICAL ONES IN THE SPECIAL CASE OF 1-STEP
##AHEAD, WHEN THE MODEL IS AN AR(1) AND THE ERROR
##IS N(0,1).

##true ar(1) parameter:
phi1 <- 0.95

##simulate:
set.seed(123)
y <- arima.sim(list(ar=phi1), 1000)

##large sample (T=1000):
##======================

##estimate ar(1):
mymodel <- arx(y, mc=FALSE, ar=1)

##predictions of the mean:
predict(mymodel, plot=TRUE)
predict(mymodel, n.ahead=24, plot=TRUE)
predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE)
predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE,
  innov=rnorm(10000*24))

##conclusion ("large sample"): the difference is small
##in the fans produced by the bootstrap and innov=rnorm.

##small sample (T=20):
##====================

##estimate ar(1):
mymodel <- arx(y[1:20], ar=1)

##predictions of the mean:
predict(mymodel, n.ahead=24, plot=TRUE) -> tmp1
predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE) -> tmp2
predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE,
  innov=rnorm(10000*24)) -> tmp3
cbind(tmp1,tmp2,tmp3)

##conclusion ("small sample"): the difference can be large
##in the fans produced by the bootstrap and innov=rnorm, or?
##are you sure, G-man?
