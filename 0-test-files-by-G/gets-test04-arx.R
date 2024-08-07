##################################################
## Test file for gets package. First created
## 23 September 2014, Oslo.
##
## 1 INITIATE
## 2 TEST MAIN arx() ARGUMENTS
## 3 MORE TESTS OF predict.arx()
## 4 TEST USER-DEFINED DIAGNOSTICS
## 5 TEST USER-DEFINED ESTIMATION
## 6 SIMULATIONS (FOR THE FUTURE)
##
##################################################

##################################################
##1 INITIATE
##################################################

##set working directory:
#
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd("C:/Users/sucarrat/Documents/R/gs/gets/20210902build/")
#setwd(choose.dir())

##load required packages:
require(parallel)
require(zoo)

##remove everything in workspace (.GlobalEnv):
rm(list=ls())

##load source:
source("./gets/R/gets-base-source.R")
source("./gets/R/gets-isat-source.R") #needed for sim()


##################################################
## 2 TEST MAIN arx() ARGUMENTS
##################################################

##bug detected in 0.35: name of regressand not
##correctly stored by arx(), solved in 0.36
set.seed(123)
x <- rnorm(50)
arx(x)

##set plot option:
##================

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##generate some data:
##===================

set.seed(123)
y.n <- 60
y <- arima.sim(list(ar=0.7),y.n)
y <- ts(y, frequency=4, end=c(2015,4))
mX <- matrix(rnorm(10*y.n), y.n, 10)
mX <- ts(mX, frequency=4, end=c(2015,4))
mX[1:5,2] <- NA
stepdum1 <- sim(y)[,floor(y.n/2)]
stepdum2 <- sim(y)[,floor(y.n/1.5)]
##DO NOT WORK!!:
#stepdum1 <- sim(y, which.ones=index(y)[c(floor(y.n/2))])
#stepdum2 <- sim(y, which.ones=index(y)[floor(y.n/1.5)])
mX <- cbind(mX, stepdum1, stepdum2)
y[1] <- NA; y[y.n] <- NA
colnames(mX) <- c(paste0("xvar",1:10), "stepdum1", "stepdum2")

##test each argument separately and together:
arx(y, mc=FALSE) #should return "Warning message: In plot.arx(out) : No estimated...etc."
arx(y, mc=FALSE, normality.JarqueB=TRUE)
arx(y, mc=TRUE)
arx(y, mc=FALSE, ar=c(1,3))
arx(y, mc=FALSE, ewma=list(length=c(2,4)))
arx(y, mc=FALSE, mxreg=mX)
arx(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)), mxreg=mX)
arx(y, mc=FALSE, vc=TRUE)
arx(y, mc=FALSE, arch=c(2,4))
arx(y, mc=FALSE, asym=c(2,4))
arx(y, mc=FALSE, log.ewma=list(length=c(3,5)))
arx(y, mc=FALSE, vxreg=cbind(log(mX[,1:2]^2), mX[,3:4]))
arx(y, mc=FALSE, vc=TRUE, arch=c(2,4), asym=c(1,3),
  log.ewma=list(length=c(3,5)),
  vxreg=cbind(log(mX[,1:2]^2), mX[,3:4]))
arx(y, mc=TRUE, ar=c(1,3), vcov.type="o")
arx(y, mc=TRUE, ar=c(1,3), vcov.type="w")
arx(y, mc=TRUE, ar=c(1,3), vcov.type="n")
arx(y, mc=TRUE, ar=c(1,3), qstat.options=c(5,5))
arx(y, mc=TRUE, ar=c(1,3), tol=1e-15)
arx(y, mc=TRUE, ar=c(1,3), tol=1, LAPACK=TRUE) #Returns error, but should it?: "Error in ols...: singular regressor-matrix"
arx(y, mc=TRUE, ar=c(1,3), tol=1, LAPACK=FALSE) #should give warning about removed regressors

##only mean specification:
mod01 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX)
print(mod01)
print(mod01, signif.stars=FALSE)
coef(mod01)
coef(mod01, spec="m")
coef(mod01, spec="v") #should be NULL
coef(mod01, spec="b") #should be the same as "m"
plot(ES(mod01, level=c(0.99,0.95,0.9))) #expected shortfall
plot(cbind(fitted(mod01),
  fitted(mod01, spec="m"),
  fitted(mod01, spec="v")))
fitted(mod01, spec="b") #should be NULL
logLik(mod01)
plot(mod01)
predict(mod01) #should return the error-message: 'newmxreg' is NULL
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)))
predict(mod01, n.ahead=30, newmxreg=matrix(0,30,NCOL(mX)))
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), spec="mean")
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),
  plot=FALSE, return=FALSE) #no plot, return nothing
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), plot=FALSE)
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),
  return=FALSE, plot=TRUE) #plot only
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),
  plot.options=list(fitted=TRUE))
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),
  plot.options=list(newmactual=rep(0,5)))
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),
  spec="variance") #should return "Set 'vc=TRUE'..etc."
predict(mod01, n.ahead=5, spec="variance") #should return "Set 'vc=TRUE'..etc."
recursive(mod01)
recursive(mod01, return=FALSE) #only plot
recursive(mod01, plot=FALSE) #only return (values)
recursive(mod01, plot=FALSE, return=FALSE) #return nothing
recursive(mod01, std.errors=FALSE)
recursive(arx(y, mc=FALSE)) #should return the error-message: No mean-equation
recursive(arx(y, mc=FALSE), spec="variance") #should return the error-message No variance-equation
recursive(arx(y, mc=TRUE, plot=FALSE))
recursive(arx(y, mc=TRUE, plot=FALSE), return=FALSE)
recursive(arx(y, mc=FALSE, ar=1, plot=FALSE))
recursive(arx(y, mc=FALSE, mxreg=mX, plot=FALSE))
plot(cbind(residuals(mod01),
  residuals(mod01, std=FALSE),
  residuals(mod01, std=TRUE)))
sigma(mod01)
rsquared(mod01)
summary(mod01)
plot(VaR(mod01, level=c(0.99,0.95,0.9))) #value-at-risk
vcov(mod01)
vcov(mod01, spec="m")
vcov(mod01, spec="v") #should return NULL

##both mean and variance specifications:
mX <- mX[,1:2]
mod02 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX,
  arch=1:4,asym=1:2, vxreg=log(mX^2))
print(mod02)
print(mod02, signif.stars=FALSE)
coef(mod02)
coef(mod02, spec="m")
coef(mod02, spec="v")
coef(mod02, spec="b")
plot(ES(mod02, level=c(0.99,0.95,0.9))) #expected shortfall
plot(cbind(fitted(mod02),
  fitted(mod02, spec="m"),
  fitted(mod02, spec="v"),
  fitted(mod02, spec="b")))
logLik(mod02)
plot(mod02)
predict(mod02) #should return the error-message 'newvxreg' is NULL
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX))) #should return the error-message 'newmxreg' is NULL
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),
  newmxreg=matrix(0,5,NCOL(mX)))
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),
  newmxreg=matrix(0,5,NCOL(mX)), plot=FALSE)
predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)))
predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),
  plot.options=list(fitted=TRUE))
predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),
  plot.options=list(newvactual=rep(1,5)))
predict(mod02, n.ahead=30, newvxreg=matrix(0,30,NCOL(mX)),
  newmxreg=matrix(0,30,NCOL(mX)))
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),
  newmxreg=matrix(0,5,NCOL(mX)), spec="mean")
predict(mod02, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),
  spec="variance") #should return the error-message 'newvxreg' is NULL
recursive(mod02)
recursive(mod02, spec="variance")
recursive(mod02, spec="variance", return=FALSE)
recursive(mod02, spec="variance", plot=FALSE, return=FALSE) #return nothing
recursive(mod02, spec="variance", std.errors=FALSE)
recursive(arx(y, vc=TRUE), spec="variance")
recursive(arx(y, arch=1), spec="variance")
recursive(arx(y, vxreg=mX), spec="variance")
recursive(arx(y, vc=TRUE), spec="variance", return=FALSE)
plot(cbind(residuals(mod02),
  residuals(mod02, std=FALSE),
  residuals(mod02, std=TRUE)))
summary(mod02)
plot(VaR(mod02, level=c(0.99,0.95,0.9))) #value-at-risk
vcov(mod02)
vcov(mod02, spec="m")
vcov(mod02, spec="v")
##test the generation of LaTeX code:
printtex(mod02)
printtex(mod02, fitted.name="g")
printtex(mod02, fitted.name="g", print.info=FALSE)
toLatex(mod02)
toLatex(mod02, fitted.name="g")
toLatex(mod02, fitted.name="g", print.info=FALSE)

##only variance specification:
mod03 <- arx(y, mc=FALSE, arch=1:4, asym=1:2, log.ewma=3, vxreg=log(mX^2))
print(mod03)
print(mod03, signif.stars=FALSE)
coef(mod03)
coef(mod03, spec="m") #should be NULL
coef(mod03, spec="v")
coef(mod03, spec="b")
plot(ES(mod01, level=c(0.99,0.95,0.9)), #expected shortfall
  plot.type="single", col=c("blue","red","green4"))
plot(cbind(fitted(mod03),
  fitted(mod03, spec="m"),
  fitted(mod03, spec="v")))
fitted(mod03, spec="b") #should be NULL
logLik(mod03)
plot(mod03)
plot(cbind(residuals(mod03),
  residuals(mod03, std=FALSE),
  residuals(mod03, std=TRUE)))
predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2))
predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2),
  plot.options=list(fitted=TRUE, newvactual=rep(1,24)))
recursive(mod03) #should return the error-message No mean-equation
recursive(mod03, spec="variance")
recursive(mod03, spec="variance", return=FALSE)
recursive(mod03, spec="variance", plot=FALSE, return=FALSE) #return nothing
recursive(mod03, spec="variance", std.errors=FALSE)
summary(mod03)
plot(VaR(mod03, level=c(0.99,0.95,0.9)), #value-at-risk
  plot.type="single", col=c("blue","red","green4"))
vcov(mod03)
##this used to erroneously be non-NULL,
##solved in version 0.28:
vcov(mod03, spec="m") #should return NULL
vcov(mod03, spec="v")


##################################################
## 3 TEST USER-DEFINED DIAGNOSTICS
##################################################

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
mod06 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX)
mod06 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.diagnostics=list(name="SWtest"))
print(mod06)
##the pval argument is ignored (as it should):
mod06 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.diagnostics=list(name="SWtest", pval=0.025))
mod06 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.diagnostics=list(name="SWtest", pval=0.99))
print(mod06)

##test the envir entry:
rm("SWtest") #make sure SWtest is not defined in the global environment
myenv <- new.env()
assign("SWtest",
  function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
  rownames(result) <- "myenv-test"
  return(result)
  },
  envir=myenv) #end function
mod06 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX, #should not work
  user.diagnostics=list(name="SWtest"))
mod06 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX, #should work
  user.diagnostics=list(name="SWtest", envir=myenv))
print(mod06)


##################################################
## 4 TEST USER-DEFINED ESTIMATION
##################################################

## Rules arx: The returned result should be a list with at least three items
## named "coefficients", "df" and "vcov". The item named "df" is used to
## compute the p-values associated with the t-statistics: coef/std.err. 

##user-defined estimator (minimal):
Gfun <- function(y, x, method=3){
  tmp <- ols(y, x, method=method)
  tmp <- list(coefficients=tmp$coefficients, df=tmp$df, vcov=tmp$vcov)
  return(tmp)
}

mod07 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.estimator=list(name="Gfun"), plot=FALSE)
summary(mod07)
print(mod07)
mod07 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX, #should work but produce warning
  user.estimator=list(name="Gfun"), plot=TRUE) 
summary(mod07)
print(mod07)
coef(mod07)
coef(mod07, spec="m")
coef(mod07, spec="v") #should be null
coef(mod07, spec="b") #mean coefs only
residuals(mod07) #should be null 
residuals(mod07, std=FALSE) #should be null
residuals(mod07, std=TRUE) #should be null
fitted(mod07) #should be NULL
fitted(mod07, spec="m") #should be NULL
fitted(mod07, spec="v") #should be NULL
fitted(mod07, spec="b") #should be NULL
logLik(mod07) #should produce warning: 'object$logl' is NULL
plot(mod07) #should return "...no plot produced"
recursive(mod07) #should return the error-message "...Not available..."
vcov(mod07)
vcov(mod07, spec="m") #should be NULL
vcov(mod07, spec="v") #should be NULL

##user-defined estimator (usual):
Gfun <- function(y, x, ...){
  tmp <- ols(y, x, method=3)
  tmp$mean.fit <- tmp$fit
  return(tmp)
}
mod08 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.estimator=list(name="Gfun"), plot=FALSE)
summary(mod08)
print(mod08)
mod08 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.estimator=list(name="Gfun"), plot=TRUE) #should return: "...no plot produced"
summary(mod08)
print(mod08)
coef(mod08)
coef(mod08, spec="m")
coef(mod08, spec="v") #should be null
coef(mod08, spec="b")
residuals(mod08) 
residuals(mod08, std=FALSE)
residuals(mod08, std=TRUE) ##should be NULL
fitted(mod08)
fitted(mod08, spec="m")
fitted(mod08, spec="v") #should be NULL
fitted(mod08, spec="b") #should be NULL
logLik(mod08)
plot(mod08) #should produce warning
##this does not work, fixing it (i.e. changing Gfun) requires some work!:
predict(mod08, n.ahead=24, newmxreg=matrix(0,24,2))
recursive(mod08) #should return the error-message "...Not available..."
vcov(mod08)
vcov(mod08, spec="m") #should return NULL
vcov(mod08, spec="v") #should return NULL

##user-defined estimator (fast ols in 'big' problems):
library(Matrix)
ols2 <- function(y, x){
  out <- list()
  out$n <- length(y)
  if (is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
  out$df <- out$n - out$k
  if (out$k > 0) {
    x <- as(as(as(x,"Matrix"), "generalMatrix"), "unpackedMatrix")
#OLD:  
#    x <- as(x, "dgeMatrix")
    out$xpy <- crossprod(x, y)
    out$xtx <- crossprod(x)
    out$coefficients <- as.numeric(solve(out$xtx,out$xpy))
    out$xtxinv <- solve(out$xtx)
    out$fit <- out$fit <- as.vector(x %*% out$coefficients)  
  }else{ out$fit <- rep(0, out$n)	}
	  out$residuals <- y - out$fit
	  out$residuals2 <- out$residuals^2
	  out$rss <- sum(out$residuals2)
	  out$sigma2 <- out$rss/out$df
	  if (out$k > 0) { out$vcov <- as.matrix(out$sigma2 * out$xtxinv) }
	  out$logl <-
	    -out$n * log(2 * out$sigma2 * pi)/2 - out$rss/(2 * out$sigma2)
	  return(out)            
}

mod09 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.estimator=list(name="ols2"), plot=FALSE)
summary(mod09)
print(mod09)
mod09 <- arx(y, mc=FALSE, ar=1:4, mxreg=mX,
  user.estimator=list(name="ols2"), plot=TRUE) #should produce warning
summary(mod09)
print(mod09)
coef(mod09)
coef(mod09, spec="m")
coef(mod09, spec="v") #should be null
coef(mod09, spec="b")
residuals(mod09)
residuals(mod09, std=FALSE)
residuals(mod09, std=TRUE) #should be null
fitted(mod09) #should be null
fitted(mod09, spec="m") #should be null
fitted(mod09, spec="v") #should be null
fitted(mod09, spec="b") #should be NULL
logLik(mod09)
plot(mod09) #should produce warning
recursive(mod09) #should return the error-message "Not available..."
vcov(mod09)
vcov(mod09, spec="m")
vcov(mod09, spec="v")

##user-defined estimator (gmm):
z <- matrix(rnorm(10*54), 54,10)
colnames(z) <- paste0("z", 1:NCOL(z))
mod10 <- arx(y, mc=TRUE, ar=1:4, mxreg=mX,
  user.estimator=list(name="gmm", z=z), plot=FALSE)
mod10


##################################################
## 5 SIMULATIONS
##################################################

##sigma() vs. log-arch estimates of intercept:
##--------------------------------------------

set.seed(123)
y <- rnorm(5000)

##for large T, they are approximately equal:
mymodel <- arx(y, vc=TRUE, plot=FALSE)
value1 <- round(sigma(mymodel)^2, digits=3)
value2 <- round(exp(coef(mymodel, spec="variance")[1]), digits=3)
value1 == value2

##for small T, they differ:
mymodel <- arx(y[1:50], vc=TRUE, plot=FALSE)
value1 <- round(sigma(mymodel)^2, digits=3)
value2 <- round(exp(coef(mymodel, spec="variance")[1]), digits=3)
value1 == value2
