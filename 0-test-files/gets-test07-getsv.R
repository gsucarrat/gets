##################################################
## Test file for gets package. First created
## 23 September 2014, Oslo.
##
## 1 INITIATE
## 2 TEST getsv() ARGUMENTS
## 3 SIMULATIONS (FOR THE FUTURE)
##
##################################################

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
source("gets-base-source.R")
source("gets-isat-source.R")


##################################################
## 2 TEST getsv() ARGUMENTS
##################################################

set.seed(123)
n <- 500
k <- 5
mX <- matrix(rnorm(k*n), n, k)
mX <- ts(mX, frequency=12, end=c(2015,12))
vX <- log(mX^2)
vX <- ts(vX, frequency=12, end=c(2015,12))
library(lgarch)
y <- lgarchSim(n, arch=c(0.2), garch=0, verbose=FALSE)
y <- ts(y, frequency=12, end=c(2015,12))
getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##only variance spec:
vgum01 <- arx(y, arch=1:3, asym=1:2, vxreg=vX)
vgum01
getsv(vgum01)
getsv(vgum01, t.pval=0.10)
getsv(vgum01, wald.pval=0.15)
getsv(vgum01, do.pet=FALSE)
getsv(vgum01, ar.LjungB=list(lag=5, pval=0.05))
getsv(vgum01, arch.LjungB=list(lag=4, pval=0.1))
getsv(vgum01, ar.LjungB=NULL, arch.LjungB=NULL)
getsv(vgum01, normality.JarqueB=0.05)
getsv(vgum01, info.method="sc")
getsv(vgum01, info.method="hq")
getsv(vgum01, keep=1:3)
getsv(vgum01, keep=c(1,3))
##issue raised by Jonas Kurle/Moritz Schwarz in email
##24 October 2019 sent to F-bear. If only a single non-keep
##regressor, then no search is undertaken. Solved in getsFun
##by G in 0.24:
getsv(vgum01, keep=1:10)
getsv(vgum01, keep=NULL) #should return 'warning': "Regressor 1..."
getsv(vgum01, include.gum=TRUE)
getsv(vgum01, include.1cut=FALSE)
getsv(vgum01, include.empty=TRUE)
getsv(vgum01, include.gum=TRUE, include.1cut=TRUE, include.empty=TRUE)
getsv(vgum01, max.paths=1)
getsv(vgum01, max.paths=3)
getsv(vgum01, tol=1) #should return the error: "Error in qr.solve...
getsv(vgum01, turbo=TRUE)
getsv(vgum01, print.searchinfo=FALSE)
getsv(vgum01, plot=FALSE)

##extraction functions:
vgets01 <- getsv(vgum01)
print(vgets01)
summary(vgets01)
coef(vgets01)
coef(vgets01, spec="m") #should be NULL
coef(vgets01, spec="v")
coef(vgets01, spec="b")
plot(ES(vgets01, level=c(0.99,0.95,0.9))) #expected shortfall
fitted(vgets01)
fitted(vgets01, spec="m") #should be 0 all over
plot(fitted(vgets01))
plot(fitted(vgets01, spec="v"))
logLik(vgets01)
paths(vgets01)
paths(vgum01) #should return the error-message: The object does not belong to the 'gets' nor 'isat' class
plot(vgets01)
plot(residuals(vgets01))
plot(residuals(vgets01, std=TRUE))
plot(residuals(vgets01, std=FALSE))
predict(vgets01)
sigma(vgets01)
rsquared(vgets01) #should return NA
terminals(vgets01)
terminals(vgum01) #should return the error-message: The object does not belong to the 'gets' nor 'isat' class
plot(VaR(vgets01, level=c(0.99,0.95,0.9))) #value-at-risk
vcov(vgets01)
vcov(vgets01, spec="m")
vcov(vgets01, spec="v")

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$std.residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
vgum01 <- arx(y, arch=1:3, asym=1:2, vxreg=vX,
  user.diagnostics=list(name="SWtest", pval=0.025))
vgum01
getsv01 <- getsv(vgum01,
  user.diagnostics=list(name="SWtest", p.value=0.025))
getsv01

##both mean and variance specs:
vgum02 <- arx(y, mc=TRUE, ar=1:2, mxreg=mX,
  arch=1:3, asym=1:2, vxreg=vX)
vgum02
vgets02 <- getsv(vgum02)
print(vgets02)
summary(vgets02)
coef(vgets02)
coef(vgets02, spec="m") #should be NULL
coef(vgets02, spec="v")
coef(vgets02, spec="b")
plot(ES(vgets02, level=c(0.99,0.95,0.9))) #expected shortfall
fitted(vgets02, spec="m") #used to return zeros, should be non-zero
plot(fitted(vgets02))
plot(fitted(vgets02, spec="v"))
#plot(fitted(vgets02, spec="b")) #should result in error
logLik(vgets02)
paths(vgets02)
paths(vgum02) #should return error: "The object does not belong to the 'gets' or 'isat' class
plot(vgets02)
plot(residuals(vgets02))
plot(residuals(vgets02, std=TRUE))
plot(residuals(vgets02, std=FALSE))
predict(vgets02)
terminals(vgets02)
terminals(vgum02) #should return the error-message: The object does not belong to the 'gets' nor 'isat' class
plot(VaR(vgets02, level=c(0.99,0.95,0.9))) #value-at-risk
vcov(vgets02)
vcov(vgets02, spec="m") #should be NULL
vcov(vgets02, spec="v")


##################################################
## 3 SIMULATIONS (FOR THE FUTURE)
##################################################
