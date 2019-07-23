##################################################
## Test file for gets package. First created
## 23 September 2014, Oslo.
##
## 1 INITIATE
## 2 TEST ols()
## 3 TEST diagnostics()
## 4 TEST regressorsMean()
##
##################################################

##################################################
##1 INITIATE
##################################################

##remove everything in workspace:
rm(list=ls())

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd(choose.dir())

##load source code:
source("gets-base-source.R")
source("gets-isat-source.R")


##################################################
##2 test ols()
##################################################

set.seed(123)
vY <- rnorm(20)
mX <- matrix(rnorm(3*20), 20, 3)

##test methods 1-5 (coefficients), compare with lm:
mCoefs <- matrix(NA,6,NCOL(mX))
ols(vY, mX, method=1)
mCoefs[1,] <- ols(vY, mX, method=1)$coefficients
ols(vY, mX, method=2)
mCoefs[2,] <- ols(vY, mX, method=2)$coefficients
ols(vY, mX, method=3)
mCoefs[3,] <- ols(vY, mX, method=3)$coefficients
ols(vY, mX, method=4)
mCoefs[4,] <- ols(vY, mX, method=4)$coefficients
ols(vY, mX, method=5)
mCoefs[5,] <- ols(vY, mX, method=5)$coefficients
mCoefs[6,] <- lm(vY ~ mX-1)$coefficients #compare with lm
rownames(mCoefs) <- paste("method ",1:NROW(mCoefs), ":", sep="")
mCoefs

##compare log-likelihoods:
tmp <- ols(vY, mX, method=3)
tmp$logl == sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE))

## until version 0.9, this example failed for
## methods 2-5 (example by F-bear, see his
## email 14/11-2016):
set.seed(123)
y <- rnorm(100, 0, 1)
x <- rnorm(100, 0, 1)*10^8
ols(y,cbind(1,x), method=2, LAPACK=FALSE)$coefficients
ols(y,cbind(1,x), method=2, LAPACK=TRUE)$coefficients

##check that method=0 yields error:
ols(vY, mX, method=0)

##test tol argument (only used if LAPACK=FALSE):
ols(vY, mX, tol=1, LAPACK=FALSE) #should return 'Error...'
ols(vY, mX, tol=1, LAPACK=TRUE) #should work
ols(vY, mX, tol=0.96, LAPACK=FALSE) #should return 'Error...'
ols(vY, mX, tol=0.96, LAPACK=TRUE) #should work
ols(vY, mX, tol=0.95, LAPACK=FALSE) #should work
ols(vY, mX, tol=0.95, LAPACK=TRUE) #should work

##view the return values of methods 1 to 6:
names(ols(vY, mX, method=1))
names(ols(vY, mX, method=2))
names(ols(vY, mX, method=3))
names(ols(vY, mX, method=4))
names(ols(vY, mX, method=5))
names(ols(vY, mX, method=6))


##################################################
##3 TEST diagnostics()
##################################################

set.seed(123)
vY <- rnorm(20)
mX <- matrix(rnorm(3*20), 20, 3)

##check x for autocorrelation and ARCH, and return a
##data-frame with the results:
x <- ols(vY, mX, method=3)
diagnostics(x)

##add the Jarque-Bera normality test to the diagnostics:
diagnostics(x, normality.JarqueB=TRUE)
diagnostics(x, normality.JarqueB=0.9)
diagnostics(x, normality.JarqueB=0.9, verbose=FALSE)

##check x for autocorrelation and ARCH, and indicate
##whether it passes the check:
diagnostics(x, verbose=FALSE)
diagnostics(x, normality.JarqueB=0.8, verbose=FALSE)

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
SWtest(x)
diagnostics(x, user.fun=list(name="SWtest", pval=0.025))
diagnostics(x, user.fun=list(name="SWtest", pval=0.025),
  verbose=FALSE)

##user-defined Shapiro-Wilks test for normality in the residuals,
##but with user-specified row-name:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
  rownames(result) <- "G-test"
  return(result)
}
diagnostics(x, user.fun=list(name="SWtest", pval=0.025))
diagnostics(x, user.fun=list(name="SWtest", pval=0.025), verbose=FALSE)

##test the envir entry:
rm("SWtest") #make sure SWtest is not defined in the global environment
myenv <- new.env()
assign("SWtest",
  function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
  rownames(result) <- "G-test"
  return(result)
  },
  envir=myenv)
diagnostics(x, user.fun=list(name="SWtest")) #should not work
diagnostics(x, user.fun=list(name="SWtest", envir=myenv)) #should work


##################################################
##4 TEST regressorsMean()
##################################################

##test 1:
##=======

set.seed(123)
iT <- 10
y <- rnorm(iT)
y <- arima.sim(list(ar=0.3), iT)
y[1] <- NA; y[iT] <- NA
mxreg <- matrix(rnorm(5*iT), iT, 5)
#mxreg <- cbind(rep(1, iT)); colnames(mX) <- "mconst"
#mxreg[1:5,2] <- NA

regressorsMean(y)
regressorsMean(log(y^2))
regressorsMean(y, mc=TRUE)
regressorsMean(y, ar=c(1,3))
regressorsMean(y, ewma=list(length=c(2,4)))
regressorsMean(y, mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.as.zoo=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, na.trim=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE, return.as.zoo=FALSE,
  na.trim=FALSE)

##test 2:
##=======

set.seed(123)
iT <- 10 #60 or 100. If iT=60, then usually no specific
y <- arima.sim(list(ar=0.3),iT)
y <- ts(y, frequency=4, end=c(2015,4))
mxreg <- matrix(rnorm(4*iT), iT, 4)
mxreg <- ts(mxreg, frequency=4, end=c(2015,4))
y[1] <- NA; y[iT] <- NA

regressorsMean(y)
regressorsMean(log(y^2))
regressorsMean(y, mc=TRUE)
regressorsMean(y, ar=c(1,3))
regressorsMean(y, ewma=list(length=c(2,4)))
regressorsMean(y, mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.as.zoo=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, na.trim=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE, return.as.zoo=FALSE,
  na.trim=FALSE)
