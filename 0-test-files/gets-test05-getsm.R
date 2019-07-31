##################################################
## Test file for gets package. First created
## 23 September 2014, Oslo.
##
## 1 INITIATE
## 2 TEST ARGUMENTS OF getsm()
## 3 TEST USER DEFINED ESTIMATION
## 4 SIMULATIONS (FOR THE FUTURE)
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
## 2 TEST ARGUMENTS OF getsm()
##################################################

set.seed(123)
y.n <- 100 #60 or 100. If y.n=60, then usually no specific
y <- arima.sim(list(ar=0.3),y.n)
y <- ts(y, frequency=4, end=c(2015,4))
mX <- matrix(rnorm(4*y.n), y.n, 4)
mX <- ts(mX, frequency=4, end=c(2015,4))
y[1] <- NA; y[y.n] <- NA
getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)
#This yields error!!:
#y <- zooreg(y, frequency=4, start=c(1990,2))

##only mean equation:
gum01 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX)
gum01
getsm01 <- getsm(gum01)
getsm(gum01, t.pval=0.01)
getsm(gum01, wald.pval=0.01)
getsm(gum01, vcov.type="o")
getsm(gum01, vcov.type="w")
getsm(gum01, vcov.type="n")
getsm(gum01, do.pet=FALSE)
getsm(gum01, ar.LjungB=list(lag=6,pval=0.20))
getsm(gum01, ar.LjungB=list(lag=NULL,pval=0.2))
getsm(gum01, ar.LjungB=NULL)
getsm(gum01, arch.LjungB=list(lag=6,pval=0.2))
getsm(gum01, arch.LjungB=list(lag=NULL,pval=0.2))
getsm(gum01, arch.LjungB=NULL)
getsm(gum01, normality.JarqueB=0.025) #to do: make sure the JB-test is included in the print!
getsm(gum01, ar.LjungB=NULL, arch.LjungB=NULL,
  normality.JarqueB=0.025)
getsm(gum01, info.method="hq")
getsm(gum01, keep=1:2)
getsm(gum01, include.gum=TRUE)
getsm(gum01, include.empty=TRUE)
getsm(gum01, max.paths=1)
getsm(gum01, max.paths=2)
suppressMessages(getsm(gum01, print.searchinfo=TRUE))
getsm(gum01, print.searchinfo=FALSE)
getsm(gum01, print.searchinfo=FALSE, estimate.specific=FALSE,
  ar.LjungB=NULL, arch.LjungB=NULL, do.pet=TRUE) #should produce warning: "No estimated model...etc."

##extraction functions (only mean equation):
getsm01 <- getsm(gum01)
suppressMessages( getsm01 <- getsm(gum01) )
print(getsm01)
sigma(getsm01)
rsquared(getsm01)
summary(getsm01)
coef(getsm01)
coef(getsm01, spec="m")
coef(getsm01, spec="v") #should be NULL
coef(getsm01, spec="b")
plot(ES(getsm01, level=c(0.99,0.95,0.9))) #expected shortfall
plot(cbind(fitted(getsm01),
  fitted(getsm01, spec="m"),
  fitted(getsm01, spec="v")))
predict(getsm01)
predict(getsm01, spec="mean")
predict(getsm01, spec="variance") #should return "Set 'vc=TRUE' to plot...
predict(getsm01, spec="both")
predict(getsm01, n.ahead=1)
predict(getsm01, newindex=13:24) #issue in the plotting (see plot), not in the prediction
predict(getsm01, return=FALSE)
predict(getsm01, plot=FALSE)
predict(getsm01, return=FALSE)
predict(getsm01, return=FALSE, plot=FALSE)
logLik(getsm01)
plot(getsm01)
plot(cbind(residuals(getsm01),
  residuals(getsm01, std=FALSE),
  residuals(getsm01, std=TRUE)))
paths(getsm01)
paths(mod01) #should return the error-message: object 'mod01' not found
recursive(getsm01)
terminals(getsm01)
terminals(mod01) #should return the error-message: object 'mod01' not found
plot(VaR(getsm01, level=c(0.99,0.95,0.9)), #value-at-risk
  plot.type="single", col=c("blue","red","green4"))
vcov(getsm01)
vcov(getsm01, spec="m")
vcov(getsm01, spec="v") #should be NULL

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
gum01 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX,
  user.diagnostics=list(name="SWtest", pval=0.025))
gum01
getsm01 <- getsm(gum01,
  user.diagnostics=list(name="SWtest", pval=0.025))
getsm01
getsm(gum01,
  user.diagnostics=list(name="SWtest", pval=1e-10))

##both mean and variance equations:
gum02 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX, arch=1:2, asym=1,
  log.ewma=list(length=3), vxreg=log(mX^2))
gum02
getsm(gum02)
getsm(gum02, t.pval=0.01)
getsm(gum02, wald.pval=0.01)
getsm(gum02, vcov.type="o")
getsm(gum02, vcov.type="w")
getsm(gum02, vcov.type="n")
getsm(gum02, do.pet=FALSE)
getsm(gum02, ar.LjungB=list(lag=6,pval=0.20))
getsm(gum02, arch.LjungB=list(lag=6,pval=0.20))
getsm(gum02, normality.JarqueB=0.025)
getsm(gum02, ar.LjungB=NULL, arch.LjungB=NULL)
getsm(gum02, ar.LjungB=NULL, arch.LjungB=NULL, normality.JarqueB=0.025)
getsm(gum02, info.method="hq")
getsm(gum02, keep=1:2)
getsm(gum02, include.gum=TRUE)
getsm(gum02, include.empty=TRUE)
getsm(gum02, max.paths=1)
getsm(gum02, max.paths=2)
getsm(gum02, print.searchinfo=FALSE)
getsm(gum02, print.searchinfo=FALSE, estimate.specific=FALSE,
  ar.LjungB=NULL, arch.LjungB=NULL, do.pet=TRUE) #should produce warning: "No estimated model...etc."

##extraction functions (both mean and variance equations):
getsm02 <- getsm(gum02)
print(getsm02)
summary(getsm02)
coef(getsm02)
coef(getsm02, spec="m")
coef(getsm02, spec="v")
coef(getsm02, spec="b")
plot(ES(getsm02, level=c(0.99,0.95,0.9)), #expected shortfall
  plot.type="single", col=c("blue","red","green4"))
plot(cbind(fitted(getsm02),
  fitted(getsm02, spec="m"),
  fitted(getsm02, spec="v")))
logLik(getsm02)
plot(getsm02)
plot(cbind(abs(residuals(getsm02)),
  sqrt(fitted(getsm02, spec="v"))), plot.type="single",
  col=c("blue","red"))
plot(cbind(residuals(getsm02),
  residuals(getsm02, std=FALSE),
  residuals(getsm02, std=TRUE)))
predict(getsm02) #should return the error: 'newvxreg' is NULL
predict(getsm02, newvxreg=matrix(0,12,4)) #should return the error: 'newmxreg' is NULL
predict(getsm02, newmxreg=matrix(0,12,2), newvxreg=matrix(0,12,4))
predict(getsm02, n.ahead=1, newmxreg=matrix(0,1,2),
  newvxreg=matrix(0,1,4))
predict(getsm02, newmxreg=matrix(0,12,2), newvxreg=matrix(0,12,4),
  plot.options=list(newmactual=rep(0,12)))
predict(getsm02, newmxreg=matrix(0,12,2), newvxreg=matrix(0,12,4),
  plot.options=list(keep=4))
predict(getsm02, spec="variance") #should return the error: 'newvxreg' is NULL
predict(getsm02, spec="variance", newvxreg=matrix(1,12,4))
predict(getsm02, spec="variance", newvxreg=matrix(1,12,4),
  n.sim=100)
predict(getsm02, spec="variance", n.ahead=1, newvxreg=matrix(1,1,4),
  n.sim=100)
paths(getsm02)
paths(mod01) #should return the error-message: object 'mod01' not found
terminals(getsm02)
terminals(mod01) #should return the error-message: object 'mod01' not found
plot(VaR(getsm02, level=c(0.99,0.95,0.9)), #value-at-risk
  plot.type="single", col=c("blue","red","green4"))
vcov(getsm02)
vcov(getsm02, spec="m")
vcov(getsm02, spec="v")

gum03 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX, arch=1:2, asym=1,
  log.ewma=list(length=3), vxreg=log(mX^2), vcov.type="w")
getsm(gum03)
gum04 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX, arch=1:2, asym=1,
  log.ewma=list(length=3), vxreg=log(mX^2), vcov.type="n")
getsm(gum04)


##################################################
## 3 TEST USER DEFINED ESTIMATION
##################################################

##user-defined estimator (too minimal for gets-modelling):
Gfun <- function(y, x, method=3){
  tmp <- ols(y, x, method=method)
  tmp <- list(coefficients=tmp$coefficients, df=tmp$df, vcov=tmp$vcov)
  return(tmp)
}
gum01 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX,
  user.estimator=list(name="Gfun"), plot=FALSE)
summary(gum01)
print(gum01)
getsm(gum01) #should not work; lacks logl, etc.

##user-defined estimator (usual, gets-modelling should work):
Gfun <- function(y, x, ...){
  tmp <- ols(y, x, method=3)
  tmp$mean.fit <- tmp$fit
  return(tmp)
}
gum01 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX,
  user.estimator=list(name="Gfun"), plot=FALSE)
summary(gum01)
print(gum01)
myspecific <- getsm(gum01)
myspecific <- suppressMessages(getsm(gum01))
summary(myspecific)


##################################################
## 4 SIMULATIONS (FOR THE FUTURE)
##################################################

