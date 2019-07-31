##################################################
## Test file for gets package. First created
## 26 June 2017, Zaragoza.
##
## 1 INITIATE 
## 2 TEST MAIN getsFun ARGUMENTS
## 3 TEST getsFun() BOOKKEEPING
## 4 TEST USER-DEFINED DIAGNOSTICS
## 5 TEST USER-DEFINED ESTIMATION
## 6 TEST USER-DEFINED GOF FUNCTION
## 7 SIMULATIONS (FOR THE FUTURE)
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
## 2 TEST MAIN getsFun ARGUMENTS
##################################################

##dgp1:
set.seed(123)
dgpn <- 40 #default: 20. Others: 40, 100
dgpk <- 20 #default: 5. Others: 20, 60
vY <- rnorm(dgpn)
mX <- matrix(rnorm(dgpn*dgpk), dgpn, dgpk)
colnames(mX) <- paste("x",1:NCOL(mX),sep="")

##check if arguments work:
getsFun(vY, mX)
getsFun(vY, mX, turbo=TRUE)
getsFun(vY, mX, user.estimator=list(name="ols", method=4))
getsFun(vY, mX, user.estimator=list(name="ols", method=5))
gumResult <- ols(vY, mX, method=3)
getsFun(vY, mX, gum.result=gumResult)
getsFun(vY, mX, t.pval=0.1)
getsFun(vY, mX, wald.pval=0.01)
getsFun(vY, mX, do.pet=FALSE)
getsFun(vY, mX, ar.LjungB=c(1,0.025))
getsFun(vY, mX, arch.LjungB=c(1,0.025))
getsFun(vY, mX, normality.JarqueB=0.05)
getsFun(vY, mX, keep=c(1,6,10))
getsFun(vY, mX, include.gum=TRUE)
getsFun(vY, mX, include.1cut=TRUE)
getsFun(vY, mX, include.empty=TRUE)
getsFun(vY, mX, include.gum=TRUE, include.1cut=TRUE, include.empty=TRUE)
getsFun(vY, mX, max.paths=1)
getsFun(vY, mX, tol=1) #should give error
getsFun(vY, mX, LAPACK=TRUE)
getsFun(vY, mX, max.regs=5) #should give error
mod01 <- getsFun(vY, mX, print.searchinfo=FALSE)
getsFun(vY, mX, alarm=TRUE)


##################################################
## 3 TEST getsFun() BOOKKEEPING
##################################################

##plain:
mod01 <- getsFun(vY, mX)
mod01$no.of.estimations
mod01$terminals.results
mod02 <- getsFun(vY, mX, turbo=TRUE)
mod02$no.of.estimations
mod02$terminals.results
tmp <- NULL
for(m in 1:length(mod02$paths)){
  tmp[m] <- all( mod01$paths[[m]]==mod02$paths[[m]] )
}
tmp
mod01$specific.spec==mod02$specific.spec
mod01$terminals.results==mod02$terminals.results

##keep=1:
mod01 <- getsFun(vY, mX, keep=1)
mod01$no.of.estimations
mod01$terminals.results
mod02 <- getsFun(vY, mX, keep=1, turbo=TRUE)
mod02$no.of.estimations
mod02$terminals.results
tmp <- NULL
for(m in 1:length(mod02$paths)){
  tmp[m] <- all( mod01$paths[[m]]==mod02$paths[[m]] )
}
tmp
mod01$specific.spec==mod02$specific.spec
mod01$terminals.results==mod02$terminals.results

##keep=c(6,10):
mod01 <- getsFun(vY, mX, keep=c(6,10))
mod01$no.of.estimations
mod01$terminals.results
mod02 <- getsFun(vY, mX, keep=c(6,10), turbo=TRUE)
mod02$no.of.estimations
mod02$terminals.results
tmp <- NULL
for(m in 1:length(mod02$paths)){
  tmp[m] <- all( mod01$paths[[m]]==mod02$paths[[m]] )
}
tmp
mod01$specific.spec==mod02$specific.spec
mod01$terminals.results==mod02$terminals.results

##ar diagnostics:
mod01 <- getsFun(vY, mX, ar.LjungB=c(1, 0.025))
mod01$no.of.estimations
mod01$terminals.results
mod02 <- getsFun(vY, mX, ar.LjungB=c(1,0.025), turbo=TRUE)
mod02$no.of.estimations
mod02$terminals.results
tmp <- NULL
for(m in 1:length(mod02$paths)){
  tmp[m] <- all( mod01$paths[[m]]==mod02$paths[[m]] )
}
tmp
mod01$specific.spec==mod02$specific.spec
mod01$terminals.results==mod02$terminals.results

##arch diagnostics:
mod01 <- getsFun(vY, mX, arch.LjungB=c(1,0.025))
mod01$no.of.estimations
mod01$terminals.results
mod02 <- getsFun(vY, mX, arch.LjungB=c(1,0.025), turbo=TRUE)
mod02$no.of.estimations
mod02$terminals.results
tmp <- NULL
for(m in 1:length(mod02$paths)){
  tmp[m] <- all( mod01$paths[[m]]==mod02$paths[[m]] )
}
tmp
mod01$specific.spec==mod02$specific.spec
mod01$terminals.results==mod02$terminals.results

##normality diagnostics:
mod01 <- getsFun(vY, mX, normality.JarqueB=0.025)
mod01$no.of.estimations
mod01$terminals.results
mod02 <- getsFun(vY, mX, normality.JarqueB=0.025, turbo=TRUE)
mod02$no.of.estimations
mod02$terminals.results
tmp <- NULL
for(m in 1:length(mod02$paths)){
  tmp[m] <- all( mod01$paths[[m]]==mod02$paths[[m]] )
}
tmp
mod01$specific.spec==mod02$specific.spec
mod01$terminals.results==mod02$terminals.results


##################################################
## 4 TEST USER-DEFINED DIAGNOSTICS
##################################################

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
getsFun(vY, mX, user.diagnostics=list(name="SWtest", pval=0.025))

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
getsFun(vY, mX, #should not work
  user.diagnostics=list(name="SWtest", pval=0.025)) 
getsFun(vY, mX, #should work
  user.diagnostics=list(name="SWtest", pval=0.025, envir=myenv))
                            

##################################################
## 5 TEST USER-DEFINED ESTIMATION
##################################################

## Rules: The returned result, a list, should contain at least six items
## named "coefficients", "df", "vcov", "logl", "n" and "k". The first
## three items are used to compute the p-values associated with the
## t-statistics coef/std.err, and the PETs. "logl", "n" and "k" are
## needed in order to compute the information criterion. Finally, the
## estimator MUST be able to handle situations where the regressor-matrix
## x is NULL or NCOL(x)=0
##
## Rules getsm: In addition to the items for arx, "logl", "n" and "k" are
## needed in order to compute the information criterion. Also, the
## estimator MUST be able to handle situations where the regressor-matrix
## x is NULL or NCOL(x)=0


##define estimator:
myEstimator <- function(y, x){ ols(y,x) }
getsFun(vY,mX, user.estimator=list(name="myEstimator"))

##user-defined estimator and diagnostics:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
getsFun(vY, mX, user.estimator=list(name="myEstimator"),
  user.diagnostics=list(name="SWtest", pval=1e-10))
getsFun(vY, mX, user.estimator=list(name="myEstimator"),
  user.diagnostics=list(name="SWtest", pval=1e-10, envir=myenv))

##faster ols:
##There are packages and routines that can be used to make
##OLS faster, e.g. the Matrix package. The code below creates
##a new function, ols2, which is essentially a copy of
##ols(y, x, method=3), but based on routines from the Matrix
##package. The R package 'microbenchmark' suggests a speed
##improvement of 10%
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

##gets w/ols2:
getsFun(vY, mX, user.estimator=list(name="ols2"))

##compare speed 1 (small T):
system.time(getsFun(vY,mX))
system.time(getsFun(vY, mX, user.estimator=list(name="ols2")))
##Conclusion: ols is faster than ols2 when T is small

##comparisons 1: w/microbenchmark, see
##https://nelsonareal.net/blog/2017/06/speeding_up_ols.html
library(microbenchmark)
microbenchmark( ols(vY,mX), ols2(vY,mX), times=10)
microbenchmark( getsFun(vY, mX),
  getsFun(vY, mX, user.estimator=list(name="ols2")),
  times=10)

##compare speed 2 (large T):
set.seed(123); vY <- rnorm(10000)
mX <- matrix(rnorm(10000*40), 10000, 40)
colnames(mX) <- paste("x",1:NCOL(mX),sep="")

system.time(getsFun(vY,mX))
system.time(getsFun(vY,mX, user.estimator=list(name="ols2")))
##Conclusion: ols2 can be substantially faster than ols when
##T is large

##S3 method for lm:
gets.lm <- function(object, ...){

  ##regressand:
  y <- as.vector(object$model[,1])
  yName <- names(object$model)[1]

  ##regressors:
  x <- as.matrix(object$model[,-1])
  xNames <- colnames(x)
  if(NCOL(x)==0){
    x <- NULL; xNames <- NULL
  }else{
    if(is.null(xNames)){
      xNames <- paste0("X", 1:NCOL(x))
      colnames(x) <- xNames
    }
  }
  
  ##is there an intercept?
  if(length(coef(object))>0){
    cTRUE <- names(coef(object))[1] == "(Intercept)"
    if(cTRUE){
      x <- cbind(rep(1,NROW(y)),x)
      xNames <- c("(Intercept)", xNames)
      colnames(x) <- xNames
    }
  }

  ##do gets:
  myspecific <- getsFun(y, x, ...)

  ##which are the retained regressors?:
  retainedXs <- xNames[myspecific$specific.spec]
  cat("Retained regressors:\n ", retainedXs, "\n")

  ##return result
  return(myspecific)

} #close gets.lm function

##simulate, again, from dgp1:
set.seed(123)
dgpn <- 40 #default: 20. Others: 40, 100
dgpk <- 20 #default: 5. Others: 20, 60
vY <- rnorm(dgpn)
mX <- matrix(rnorm(dgpn*dgpk), dgpn, dgpk)
colnames(mX) <- paste("x",1:NCOL(mX),sep="")

##estimate lm model, do gets:
mymodel <- lm(vY ~ mX)
mymodel
myspecific <- gets(mymodel)
myspecific


##################################################
## 6 TEST USER-DEFINED GOF FUNCTION
##################################################

getsFun(vY, mX, gof.function=list(name="infocrit", method="aic"))
##user-defined gof-function:
myGofFun <- function(x, ...){ return( x$logl ) }
getsFun(vY, mX, gof.function=list(name="myGofFun"), gof.method="max")
getsFun(vY, mX, gof.method="min")
getsFun(vY, mX, gof.method="max")
##user-defined gof-function:
myGofFun <- function(x, ...){ return( x$k ) }
getsFun(vY, mX, gof.function=list(name="myGofFun"), gof.method="min")

#user-defined gof and diagnostics:
getsFun(vY, mX, gof.function=list(name="myGofFun"),
  user.diagnostics=list(name="SWtest", pval=0.05))

#user-defined gof, diagnostics and estimator:
getsFun(vY, mX, gof.function=list(name="myGofFun"),
  user.diagnostics=list(name="SWtest", pval=0.05),
  user.estimator=list(name="ols2"))

##adjusted R-squared:
myGofFun <- function(object){
  yvar <- object$fit + object$residuals
  TSS <- sum( (yvar - mean(yvar))^2 )
  RSS <- sum(object$residuals^2)
  Rsquared <- 1 - RSS/TSS
  result <- 1 - (1-Rsquared)*(object$n-1)/(object$n-object$k)
  return(result)
}
##do gets while maximising R-squared:
getsFun(vY, mX, gof.function=list(name="myGofFun"), gof.method="max")


##################################################
## 7 SIMULATIONS (FOR THE FUTURE)
##################################################
