##################################################
## Test file for the 'isat' source of the 'gets'
## package. First created 23 September 2014, Oslo.
##
## Current version: 0.24
##
## 1 INITIATE
## 2 TEST iim(), sim() AND tim()
## 3 TEST MAIN isat() ARGUMENTS
## 4 TEST USER-DEFINED DIAGNOSTICS
## 5 TEST USER-DEFINED ESTIMATION
## 6 TEST USER-DEFINED GOF FUNCTION
## 7 TEST PARALLEL COMPUTING
## 8 TEST ROBUST COEFFICIENT COVARIANCES
## 9 SIMULATIONS
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
## 2 TEST iim(), sim() AND tim()
##################################################

x <- 6
iim(x); sim(x); tim(x)
which.ones <- c(2,5)
iim(x, which.ones=which.ones)
sim(x, which.ones=which.ones)
tim(x, which.ones=which.ones)

set.seed(123)
x <- ts(rnorm(5), start=2010)
iim(x); sim(x); tim(x)

set.seed(123)
x <- ts(rnorm(6), frequency=4, end=c(2015,4))
iim(x); sim(x); tim(x)
which.ones <- c(2014.25,2014.75,2015,2015.50)
iim(x, which.ones=which.ones)
sim(x, which.ones=which.ones)
tim(x, which.ones=which.ones)

##the following commands should return the
##error-message: 'which.ones' not in index
which.ones <- c(2001.25,2002.50,2003.75,2004)
iim(x, which.ones=which.ones)
sim(x, which.ones=which.ones)
tim(x, which.ones=which.ones)

##used to yield error:
set.seed(123); y <- rnorm(30); x <- as.vector(sim(y, which.ones=15))
isat(y, mxreg=x, LAPACK=FALSE)

##some ideas for the future?:
#y <- ymon <- yqtr <- rnorm(24)
###zooreg:
#y <- zooreg(y, frequency=4, start=c(2000,1))
#index(y)
#
###quarterly:
#yqtr <- yearqtr(yqtr, start=c(2000,1))
#
#index(y)
#tmp <- as.yearqtr(index(y))
#format(tmp, format="%YQ%q")
#iim(y)


##################################################
## 3 TEST MAIN isat() ARGUMENTS
##################################################

##control the plotting:
##=====================

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##generate some data:
##===================

set.seed(123); dgp.n <- 50
y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
y[1:10] <- y[1:10] + 4 #step-shift
plot(as.zoo(y), col=4) #plot

##some basic tests:
##=================

isat(y, sis=TRUE)
isat(y, sis=FALSE) #should return the error: "Error in isat..."
isat(y, iis=TRUE, sis=TRUE)
isat(y, iis=TRUE, sis=FALSE)
isat(y, sis=FALSE, tis=TRUE)
isat(y, sis=TRUE, tis=TRUE)
isat(y, iis=TRUE, sis=FALSE, tis=TRUE)
isat(y, iis=TRUE, sis=TRUE, tis=TRUE)

isat(y, ar=1:2, sis=TRUE)
isat(y, ar=1:2, iis=TRUE, sis=TRUE)
isat(y, ar=1:2, iis=TRUE, sis=FALSE)
isat(y, ar=1:2, iis=FALSE, sis=FALSE, tis=TRUE)
isat(y, ar=1:2, iis=TRUE, sis=FALSE, tis=TRUE)
isat(y, ar=1:2, iis=FALSE, sis=TRUE, tis=TRUE)
isat(y, ar=1:2, iis=TRUE, sis=TRUE, tis=TRUE)

isat(y, ar=1:2, mxreg=mX)
isat(y, ar=1:2, mxreg=mX, iis=TRUE)
isat(y, ar=1:2, mxreg=mX, iis=FALSE, sis=TRUE)
isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE)
isat(y, ar=1:2, mxreg=mX, iis=FALSE, tis=TRUE)
isat(y, ar=1:2, mxreg=mX, iis=TRUE, tis=TRUE)
isat(y, ar=1:2, mxreg=mX, iis=FALSE, sis=TRUE, tis=TRUE)
isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE)

##yielded error in version 0.9 to 0.23:
isat(y, ar=1:2, mxreg=as.data.frame(mX))

##issue reported by F-bear regarding version 0.14 (email 30/3-2018).
##in version 0.14 the impulse dummies were not detected:
set.seed(123); y <- rnorm(100, 0, 1)
y[30] <- y[30]+10; y[40] <- y[40]+10; y[60] <- y[60]-10
plot(as.zoo(y))
isat(y, iis=TRUE, sis=FALSE, t.pval=0.05, plot=TRUE)
z <- rnorm(100, 0, 1)
z[1:30] <- z[1:30] + 10
isat(z, iis=TRUE,  sis=FALSE, t.pval=0.05, plot=TRUE)

##issue reported by F-bear regarding version 0.20 (email 26/9-2019).
## "...there seems to be a serious bug in the latest version of the package.
## ISnames seems to be null, even if there are impulses retained. This
## seems to break a lot of other functions building on the ISnames element.
## For example, this does not work and returns null: [fixed by G in 0.21]"
set.seed(123)
y <- rnorm(100, 0, 1)
my <- isat(y, iis=TRUE, sis=FALSE, t.pval=0.05)
my
my$ISnames


##test the extraction functions:
##==============================

##same data as earlier:
set.seed(123); dgp.n <- 50
y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
y[1:10] <- y[1:10] + 4 #step-shift
plot(as.zoo(y), col=4) #plot

isatmod <- isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE)
print(isatmod)
summary(isatmod)
coef(isatmod)
plot(cbind(fitted(isatmod),
  fitted(isatmod, spec="m"),
  fitted(isatmod, spec="v")))
logLik(isatmod)
plot(cbind(residuals(isatmod),
  residuals(isatmod, std=FALSE),
  residuals(isatmod, std=TRUE)))
paths(isatmod)
paths(mod01) #should return the error-message: object 'mod01' not found
plot(isatmod)
predict(isatmod) #should return the error-message: 'newmxreg' is NULL
predict(isatmod, newmxreg=matrix(0,12,5))
predict(isatmod, n.ahead=1, newmxreg=matrix(0,1,5)) #used to yield error
predict(isatmod, newmxreg=matrix(0,12,5), newindex=13:24)
predict(isatmod, newmxreg=matrix(0,12,5), return=FALSE)
predict(isatmod, newmxreg=matrix(0,12,5), plot=FALSE)
predict(isatmod, newmxreg=matrix(0,12,5), return=FALSE, plot=FALSE)
terminals(isatmod)
terminals(mod01) #should return the error-message: object 'mod01' not found
vcov(isatmod)

##additional test of predict.isat():
##==================================

predict(isatmod, newmxreg=matrix(0,12,5), ci.levels=seq(0.20,0.95,by=0.05),
  n.sim=20000)
predict(isatmod, newmxreg=matrix(0,12,5), ci.levels=seq(0.20,0.95,by=0.05),
  n.sim=20000, plot.options=list(shades=seq(20,95,by=5)))
predict(isatmod, newmxreg=matrix(0,12,5), ci.levels=seq(0.20,0.95,by=0.05),
  n.sim=20000, plot.options=list(shades=seq(95,20,by=-5)))
predict(isatmod, newmxreg=matrix(0,12,5), ci.levels=seq(0.20,0.95,by=0.05),
  n.sim=20000, plot.options=list(shades=seq(100,25,by=-5)))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(keep=1))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(line.at.origin=FALSE))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(start.at.origin=TRUE))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(start.at.origin=TRUE, fitted=FALSE))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(dot.at.origin=FALSE))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(hlines=c(-2,0,2,4,6,8)))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(col=c("darkred","green")))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(lty=c(3,2)))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(lwd=3))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(ylim=c(-8,16)))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(ylab="G-values"))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(main="Plot slightly lower when 'main' is specified"))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(legend.text=c("Prognose","Faktisk")))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(fitted=FALSE))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(newmactual=rep(0,6)))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(shades=c(95,50)))
predict(isatmod, newmxreg=matrix(0,12,5),
  plot.options=list(shades=c(50,95))) #invert shades

##In the following model (isatmod), the constant was not correctly
##named 'mconst' at one point. Instead, it was named 'mxreg',
##which created problems for predict.isat: predict(isatmod). Issue
##solved 31/7/2019.
set.seed(123)
y <- rnorm(30)
isatmod <- isat(y)
isatmod
predict(isatmod, plot=TRUE)

##same y, but slightly different model:
isatmod <- isat(y, ar=1)
isatmod
predict(isatmod, plot=TRUE)

##test further arguments:
##=======================

isat(y, t.pval=0.1) #default: 0.001
isat(y, do.pet=TRUE) #default: FALSE
isat(y, wald.pval=0.1) #default: 0.001
isat(y, ar.LjungB=list(lag=NULL, pval=0.01)) #default: NULL
#no search, because the gum does not pass arch-diagnostics:
isat(y, arch.LjungB=list(lag=NULL, pval=0.01)) #default: NULL
isat(y, normality.JarqueB=0.025) #default: NULL
isat(y, info.method="aic") #default: sc
#should return warning:
isat(y, include.gum=TRUE) #default: NULL
isat(y, include.1cut=TRUE) #default: FALSE
isat(y, include.empty=TRUE) #default: FALSE
isat(y, max.paths=3) #default: NULL (i.e. "multi-path")
isat(y, max.paths=5) #default: NULL (i.e. "multi-path")
isat(y, parallel.options=NULL) #default: NULL
isat(y, parallel.options=2) #default: NULL
isat(y, parallel.options=5) #default: NULL
isat(y, parallel.options=4) #default: NULL
isat(y, parallel.options=2, max.paths=2)
isat(y, turbo=TRUE)
isat(y, parallel.options=2, max.paths=2, turbo=TRUE) #default: NULL

##test uis argument:
##==================

##uis as matrix:
uis <- iim(dgp.n)
isat(y, sis=FALSE, uis=uis)
isat(y, sis=FALSE, uis=uis, max.paths=1)
uis <- sim(dgp.n)
##as of August 2020, these do not work (why?):
isat(y, sis=FALSE, uis=uis)
isat(y, sis=FALSE, uis=uis[,seq(1,dgp.n,2)])
isat(y, sis=FALSE, uis=uis, max.paths=1)
uis <- tim(dgp.n)
isat(y, sis=FALSE, uis=uis)
isat(y, sis=FALSE, uis=uis[,seq(1,dgp.n,2)])
isat(y, sis=FALSE, uis=uis, max.paths=1)

##used to crash (uis is a matrix):
set.seed(123); y <- rnorm(30); z <- rnorm(30)
mX <- matrix(rnorm(1*30), 30, 1)
isat(y, mxreg=z, iis=FALSE, sis=TRUE, uis=mX)

##used to crash (uis is a matrix):
##as of August 2020, these do not work (why?):
dgpN <- 50
set.seed(123)
y <- rnorm(dgpN); x <- rnorm(dgpN)
x_mis <- sim(dgpN)*x
colnames(x_mis) <- paste("mis", seq(2:(NCOL(x_mis)+1)), sep="")
isat(y, ar=1, mxreg=x, sis=FALSE, uis = x_mis, t.pval=0.05, plot=TRUE)
isat(y, ar=1, mxreg=x, sis=FALSE, uis = x_mis, t.pval=0.05,
  max.paths=1, plot=TRUE)

##used to yield error because NCOL(mX) > length(y):
set.seed(123); y <- rnorm(20); mX <- matrix(rnorm(20*40), 20, 40);
isat(y, sis=FALSE, uis=mX)

##uis as list:
uis <- list(sis=sim(y) , tis=tim(y))
##as of August 2020, these do not work (why?):
isat(y, sis=FALSE, uis=uis)
isat(y, sis=FALSE, uis=uis, max.paths=1)
isat(y, sis=FALSE, uis=uis, max.paths=2)

##test blocks argument:
##=====================

##same data as earlier:
set.seed(123); dgp.n <- 50; y <- rnorm(dgp.n)

myblocks <- list()
myblocks[[1]] <- list(10:20, 30:40)
isat(y, iis=TRUE, sis=FALSE, tis=FALSE, uis=FALSE,
  blocks=myblocks)
isat(y, iis=TRUE, sis=FALSE, tis=FALSE, uis=FALSE,
  blocks=myblocks, max.paths=1)
isat(y, iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE,
  blocks=myblocks)
isat(y, iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE,
  blocks=myblocks, max.paths=1)
uis <- list(sim(dgp.n), tim(dgp.n))
myblocks[[2]] <- list(7:19, 27:34, 40:45)
##as of August 2020, these do not work (why?):
isat(y, iis=FALSE, sis=FALSE, tis=FALSE, uis=uis,
  blocks=myblocks)
isat(y, iis=FALSE, sis=FALSE, tis=FALSE, uis=uis,
  blocks=myblocks, max.paths=1)

##further tests of predict.isat:
##==============================

##issue reported by Steven Sabol (email 27/01-2017). The
##code should work as of 0.11:
set.seed(123)
mxreg <- zooreg(matrix(rnorm(700),ncol=7),start = 2002 ,frequency = 12)
colnames(mxreg) <- c("c1","c2","c3","c4","c5","c6","c7")
y = zooreg(rnorm(88),start = 2002 ,frequency = 12)
isat_mod <- isat(y, mxreg = mxreg, mc =TRUE,ar = 4,
  sis = TRUE,t.pval = 0.01,vcov.type = "white")
newmxreg  <- tail(na.trim(mxreg),12)
new_index <- index(tail(na.trim(mxreg),12))
##as of 17 July 2019, does not work:
prediction_isat <- predict.isat(isat_mod,newmxreg=newmxreg,
  n.ahead=12,newindex = new_index, return = TRUE)
prediction_isat


##tests of biascorr, isattest, isatvar, ...etc.:
##==============================================

##issue reported by Gareth Thomas (EViews, email 22/5-2017):
##"The following code produces a blank graphics output.  It seems as
##though isattest calls the R graphics even if plot=FALSE is set".
##Should be solved as of 0.13.
set.seed(123)
d <- matrix(0,100,1)
d[35:55] <- 1
e <- rnorm(100, 0, 1)
y <- d*2 +e
##Static Test against hnull=0 using bias-correction:
ys <- isat(y, sis=TRUE, iis=FALSE, tis=FALSE, t.pval=0.01, plot=FALSE)
isattest(ys, hnull=0, lr=FALSE, ci.pval = 0.99, plot = FALSE,
  biascorr=TRUE)


##################################################
## 4 TEST USER-DEFINED DIAGNOSTICS
##################################################

##generate some data:
set.seed(123); dgp.n <- 50; y <- rnorm(dgp.n)
y[1:10] <- y[1:10] + 4
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
isat(y, user.diagnostics=list(name="SWtest", pval=1e-10))

##test the envir entry:
rm("SWtest") #make sure SWtest is not defined in the global environment
myenv <- new.env()
assign("SWtest",
  function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
  rownames(result) <- "SWtest"
  return(result)
  },
  envir=myenv) #close assign
isat(y, #should not work
  user.diagnostics=list(name="SWtest", pval=0.025)) 
isat(y, #should work
  user.diagnostics=list(name="SWtest", pval=1e-05, envir=myenv))
isat(y, #should work
  user.diagnostics=list(name="SWtest", pval=0.025, envir=myenv))


##################################################
## 5 TEST USER-DEFINED ESTIMATION
##################################################

##define estimator:
myEstimator <- function(y, x){ ols(y,x) }

##isat w/user-defined estimator:
isat(y, user.estimator=list(name="myEstimator"))

##isat w/user-defined estimator and test for normality:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
isat(y, user.estimator=list(name="myEstimator"),
  user.diagnostics=list(name="SWtest", pval=1e-10))

##faster ols?:
##There are packages and routines that make OLS faster in
##certain situations, e.g. the Matrix package. The code below
##creates a new function, ols2, which is essentially a copy
##of ols(y, x, method=3), but based on routines from the Matrix
##package.
library(Matrix)
ols2 <- function(y, x){
  out <- list()
  out$n <- length(y)
  if (is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
  out$df <- out$n - out$k
  if (out$k > 0) {
    x <- as(x, "dgeMatrix")
    out$xpy <- crossprod(x, y)
    out$xtx <- crossprod(x)
    out$coefficients <- as.numeric(solve(out$xtx,out$xpy))
    out$xtxinv <- solve(out$xtx)
    out$fit <- out$fit <- as.vector(x %*% out$coefficients)
  }else{
    out$fit <- rep(0, out$n)
  }
  out$residuals <- y - out$fit
  out$residuals2 <- out$residuals^2
  out$rss <- sum(out$residuals2)
  out$sigma2 <- out$rss/out$df
  if(out$k > 0){ out$vcov <- as.matrix(out$sigma2 * out$xtxinv) }
  out$logl <-
    -out$n * log(2 * out$sigma2 * pi)/2 - out$rss/(2 * out$sigma2)
  return(out)
}

##isat w/ols2:
isat(y, user.estimator=list(name="ols2"))

##compare speed 1:
system.time(isat(y))
system.time(isat(y, user.estimator=list(name="ols2")))
##Conclusion: here, ols is faster than ols2

##comparisons 2: w/microbenchmark, see
##https://nelsonareal.net/blog/2017/06/speeding_up_ols.html
library(microbenchmark)
microbenchmark( ols(y,mX), ols2(y,mX), times=10)
microbenchmark( isat(y),
  isat(y, user.estimator=list(name="ols2")),
  times=10)

##compare speed 2:
set.seed(123); dgp.n <- 1000; y <- rnorm(dgp.n)
system.time(isat(y))
system.time(isat(y, user.estimator=list(name="ols2")))
##Conclusion: sample size matters, additional experiments
##suggests the speed increase is increasing in sample size


##################################################
## 6 TEST USER-DEFINED GOF FUNCTION
##################################################

##generate some data:
set.seed(123); dgp.n <- 50; y <- rnorm(dgp.n)
y[1:10] <- y[1:10] + 4
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)

##user-defined gof-function:
myGof <- function(object){ infocrit(object) }

##do isat:
isat(y, gof.function=list(name="myGof"))
isat(y, user.estimator=list(name="myEstimator"),
  gof.function=list(name="myGof"))
isat(y, user.diagnostics=list(name="SWtest", pval=1e-10),
  user.estimator=list(name="myEstimator"),
  gof.function=list(name="myGof"))

##adjusted R-squared:
myGof <- function(object){
  yvar <- object$fit + object$residuals
  TSS <- sum( (yvar - mean(yvar))^2 )
  RSS <- sum(object$residuals^2)
  Rsquared <- 1 - RSS/TSS
  result <- 1 - (1-Rsquared)*(object$n-1)/(object$n-object$k)
  return(result)
}

##do isat while maximising R-squared:
isat(y, gof.function=list(name="myGof"), gof.method="max")

##minimise the number of parameters/regressors:
myGof <- function(x, ...){ return( x$k ) }
isat(y, gof.function=list(name="myGof"), gof.method="min")


##################################################
## 7 TEST PARALLEL COMPUTING
##################################################

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
isat(y, user.diagnostics=list(name="SWtest", pval=1e-10))
isat(y, user.diagnostics=list(name="SWtest", pval=1e-10),
  parallel.options=2)

##user-defined estimator:
myEstimator <- function(y, x){ ols(y,x) }
isat(y, user.estimator=list(name="myEstimator"))
isat(y, user.estimator=list(name="myEstimator"),
  parallel.options=2)

##user-defined gof:
myGof <- function(x, ...){ return( x$k ) }
isat(y, gof.function=list(name="myGof"), gof.method="min")
isat(y, gof.function=list(name="myGof"), gof.method="min",
  parallel.options=2)

##all three user-defined:
isat(y, user.diagnostics=list(name="SWtest", pval=1e-10),
  user.estimator=list(name="myEstimator"),
  gof.function=list(name="myGof"), gof.method="min",
  parallel.options=2)


##################################################
## 8 TEST ROBUST COEFFICIENT COVARIANCES
##################################################

##"white" and "newey-west" coefficient covariances generally lead to
##either outright errors, or strange results. Currently, therefore,
##until more numerically stable versions are derived, users are
##discouraged to use these robust coefficient covariances. The tests
##here therefore only tests whether the arguments work or not. The
##last set of tests suggests the source of the problem is near-zero
##residuals.

##test "white" vcov (all yield errors?):
isat(y, iis=TRUE, sis=FALSE, vcov.type="w")
isat(y, iis=FALSE, sis=TRUE, vcov.type="w")
isat(y, iis=TRUE, sis=TRUE, vcov.type="w")
isat(y, iis=FALSE, tis=TRUE, vcov.type="w")
isat(y, iis=TRUE, tis=TRUE, vcov.type="w")
isat(y, iis=FALSE, sis=TRUE, vcov.type="w", tis=TRUE, plot=plotarg)
isat(y, iis=TRUE, sis=TRUE, vcov.type="w", tis=TRUE, plot=plotarg)

##test "newey-west" vcov:
set.seed(123); dgp.n <- 100
y <- rnorm(dgp.n)
#y <- rt(dgp.n, df=4.1)
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
plotarg <- FALSE

##test "newey-west" vcov:
isat(y, iis=TRUE, sis=FALSE, vcov.type="n", plot=plotarg)
#these ones yields errors!:
isat(y, iis=FALSE, sis=TRUE, vcov.type="n")
isat(y, iis=TRUE, sis=TRUE, vcov.type="n")
isat(y, iis=FALSE, tis=TRUE, vcov.type="n")
isat(y, iis=TRUE, tis=TRUE, vcov.type="n")
isat(y, iis=FALSE, sis=TRUE, vcov.type="n", tis=TRUE,
  plot=plotarg)
isat(y, iis=TRUE, sis=TRUE, vcov.type="w", tis=TRUE,
  plot=plotarg)

##code that sheds light on what the possible source of
##the problem is (near-zero residuals)
set.seed(123)
y <- rnorm(100, 0, 1)
isat(y, vcov.type=c("newey-west"))
isat(y, sis=FALSE, uis=sim(y, which.ones=1:10))
isat(y, sis=FALSE, uis=sim(y, which.ones=1:10),
  vcov.type="newey-west")
isat(y, sis=FALSE, uis=sim(y, which.ones=seq(1,20,2)),
  vcov.type="newey-west")


##################################################
## 9 SIMULATIONS
##################################################

##For the future...