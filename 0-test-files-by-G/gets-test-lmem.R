##################################################
## Test file for gets-test-lmem-source.R. First
## created 20 July 2024, Oslo.
##
## 1 INITIATE
## 2 TEST lmem() ARGUMENTS
## 3 CHECK PREDICTIONS
## 4 EXTRA TESTS OF gets.lmem() FUNCTION
## 5 SIMULATIONS (FOR THE FUTURE)
##
##################################################


##################################################
##1 INITIATE
##################################################

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd("C:/Users/sucarrat/Documents/R/gs/gets/20210902build/")
#setwd(choose.dir())

##load required packages:
require(parallel)
require(zoo)
library(gets)

##remove everything in workspace (.GlobalEnv):
rm(list=ls())

##load source:
#source("./gets/R/gets-base-source.R")
#source("./gets/R/gets-larch-source.R")
#source("./gets/R/gets-lmem-source.R")


##################################################
## 2 TEST lmem() ARGUMENTS
##################################################

##set plot option:
##================

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##generate some data:
##===================

##generate some data:
set.seed(123)
n <- 40
y <- rchisq(n, df=1)
y[1] <- NA; y[length(y)] <- NA
y[c(3,11)] <- 0
vX <- matrix(rnorm(length(y)*4), length(y), 4)

##test each argument separately and together:
lmem(y, intercept=FALSE) #should return "Error in lmem(y, intercept = FALSE) : argument 'intercept' must be TRUE"
lmem(y, intercept=TRUE)
lmem(y, arch=1)
lmem(y, arch=c(1,3))
lmem(y, log.ewma=c(2,4))
lmem(y, xreg=vX)
lmem(y, zero.adj=1)
lmem(y, vcov.type ="hac") 
lmem(y, qstat.options=3)
lmem(y, arch=1, tol=1e-15) 
lmem(y, arch=1, tol=1) 
lmem(y, arch=1, tol=1, singular.ok=FALSE) #should return "Error in ols(loge2, x, tol = tol, method = 2) : singular regressor-matrix" 
lmem(y, arch=1, log.ewma=c(3,5), xreg=vX, vcov.type="hac")

##test methods and extraction functions:
mod01 <- lmem(y, arch=1, log.ewma=c(3,5), xreg=vX, vcov.type="hac")
print(mod01)
print(mod01, signif.stars=FALSE)
coef(mod01)
fitted(mod01)
gets(mod01)
logLik(mod01)
model.matrix(mod01)
model.matrix(mod01, response=TRUE)
model.matrix(mod01, as.zoo=FALSE)
model.matrix(mod01, response=TRUE, as.zoo=FALSE)
nobs(mod01)
predict(mod01, newxreg=matrix(0,12,NCOL(vX)))
residuals(mod01)
mean(residuals(mod01))
summary(mod01)
toLatex(mod01)
vcov(mod01)

##test predict.lmem() on each argument:
predict(lmem(y))
predict(lmem(y, arch=1))
predict(lmem(y, log.ewma=5))
predict(lmem(y, xreg=vX), newxreg=matrix(0,12,NCOL(vX)))

##test predict.lmem() on each argument with n.ahead=1:
predict(lmem(y), n.ahead=1)
predict(lmem(y, arch=1), n.ahead=1)
predict(lmem(y, log.ewma=5), n.ahead=1)
predict(lmem(y, xreg=vX), newxreg=matrix(0,1,NCOL(vX)), n.ahead=1)


##################################################
## 3 CHECK VARIANCE PREDICTIONS
##################################################

##generate some data:
##===================

##small dgp:
set.seed(123)
vY <- rchisq(20, df=1)

##arch(0) model w/constant:
##=========================

mymodel <- lmem(vY)

##predictions of the variance:
functionVals <- predict(mymodel, n.ahead=3, n.sim=2)

##correct predictions:
yhat1 <- yhat2 <- yhat3 <- exp(coef(mymodel))
correctVals <- c(yhat1,yhat2,yhat3)

##do they correspond?:
all( functionVals == correctVals )


##################################################
## 4 EXTRA TESTS OF gets.lmem() FUNCTION
##################################################

##test extraction functions:
mod01 <- lmem(y, arch=1, log.ewma=c(3,5), xreg=vX, vcov.type="hac")
spec01 <- gets(mod01)
print(spec01)
print(spec01, signif.stars=FALSE)
coef(spec01)
fitted(spec01)
gets(spec01)
logLik(spec01)
model.matrix(spec01)
nobs(spec01)
predict(spec01) #should return error: "Error in predict.lmem(spec01) : 'newxreg' is NULL"
predict(spec01, newxreg=matrix(0,12,1))
residuals(spec01)
mean(residuals(spec01))
summary(spec01)
toLatex(spec01)
vcov(spec01)

##test predict.lmem():
predict(spec01, newxreg=matrix(0,1,1), n.ahead=1)
predict(spec01, newxreg=matrix(0,12,1), newindex=11:22)
predict(spec01, newxreg=matrix(0,12,1), n.sim=1)
predict(spec01, newxreg=matrix(0,12,1), n.sim=1, innov=rnorm(12)) #should return error
predict(spec01, newxreg=matrix(0,12,1), n.sim=1, innov=rchisq(12,df=1))
predict(spec01, newxreg=matrix(0,12,1), probs=c(0.25,0.75))
predict(spec01, newxreg=matrix(0,12,1), quantile.type=1, probs=c(0.25,0.75))
predict(spec01, newxreg=matrix(0,12,1), n.sim=10, verbose=TRUE)
predict(spec01, newxreg=matrix(0,12,1), n.sim=10, probs=c(0.25,0.75),
  verbose=TRUE)

##test 1-cut modelling:
specific <- gets(mod01, max.paths=0, do.pet=FALSE, include.1cut=TRUE)
print(specific)
specific <- gets(mod01, max.paths=0, do.pet=TRUE, include.1cut=TRUE)
print(specific)


##################################################
## 5 SIMULATIONS (FOR THE FUTURE)
##################################################

TBA