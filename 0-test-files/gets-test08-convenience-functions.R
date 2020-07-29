##################################################
## Test file for gets package. First created
## 14 November 2018, Oslo.
##
## 1 INITIATE
## 2 TEST eviews() AND stata()
## 3 TEST printtex()
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
## 2 TEST eviews() AND stata()
##################################################

set.seed(123)
y.n <- 40
inflation <- arima.sim(list(ar=0.3),y.n)
inflation <- ts(inflation, frequency=4, end=c(2015,4))
mX <- matrix(rnorm(4*y.n), y.n, 4)
mX <- ts(mX, frequency=4, end=c(2015,4))

mod01 <- arx(inflation, mc=TRUE, ar=1:4, mxreg=mX,
  vcov.type="w")
mod01
eviews(mod01)
suppressMessages(eviews(mod01)) #should not print anything
eviews(mod01, file='getsdata.csv')
eviews(mod01, file='getsdata.csv', print=FALSE)
eviews(mod01, print=FALSE, return=TRUE)
stata(mod01)
suppressMessages(stata(mod01))
stata(mod01, file='getsdata.csv')
stata(mod01, file='getsdata.csv', print=FALSE)
stata(mod01, print=FALSE, return=TRUE)

gets01 <- getsm(mod01)
gets01
eviews(gets01)
suppressMessages(eviews(gets01))
eviews(gets01, file='getsdata.csv')
eviews(gets01, print=FALSE)
eviews(gets01, print=FALSE, return=TRUE)
stata(gets01)
suppressMessages(stata(gets01))
stata(gets01, file='getsdata.csv')
stata(gets01, print=FALSE)
stata(gets01, print=FALSE, return=TRUE)

isat01 <- isat(inflation, iis=TRUE, sis=TRUE)
isat01
eviews(isat01)
suppressMessages(eviews(isat01))
eviews(isat01, print=FALSE, return=TRUE)
stata(isat01)
suppressMessages(stata(isat01))
stata(isat01, print=FALSE, return=TRUE)

mod02 <- arx(inflation, mc=TRUE, ar=1:4, mxreg=mX, vcov.type="n")
mod02
eviews(mod02)
suppressMessages(eviews(mod02))
eviews(mod02, file='getsdata.csv', print=FALSE)
eviews(mod02, print=FALSE, return=TRUE)
stata(mod02)
suppressMessages(stata(mod02))
stata(mod02, file='getsdata.csv', print=FALSE)
stata(mod02, print=FALSE, return=TRUE)


##################################################
## 3 TEST printtex()
##################################################

##arx-object:
printtex(mod01)
printtex(mod01, fitted.name="gvar")
printtex(mod01, xreg.names=paste0("gman", 1:9))
printtex(mod01, xreg.names=paste0("gman", 1:8))
printtex(mod01, digits=6)
printtex(mod01, intercept=FALSE)
printtex(mod01, intercept=2)
printtex(mod01, gof=FALSE)
printtex(mod01, diagnostics=FALSE)

##gets-object:
printtex(gets01)
printtex(gets01, fitted.name="gvar")
printtex(gets01, xreg.names=paste0("gman", 1:9))
printtex(gets01, digits=6)
printtex(gets01, intercept=FALSE)
printtex(gets01, gof=FALSE)
printtex(gets01, diagnostics=FALSE)

##isat-object:
printtex(isat01)
printtex(isat01, fitted.name="gvar")
printtex(isat01, xreg.names="gman")
printtex(isat01, xreg.names=paste0("gman", 1:9))
printtex(isat01, digits=6)
printtex(isat01, intercept=FALSE)
printtex(isat01, gof=FALSE)
printtex(isat01, diagnostics=FALSE)

###=========================
###test that as.lm() works:
#
set.seed(123)
y <- rnorm(30)
arxmod <- arx(y, mc=TRUE, ar=1:3)
as.lm(arxmod)

getsmod <- getsm(arxmod, keep=1)
as.lm(getsmod)

isatmod <- isat(y)
as.lm(isatmod)
