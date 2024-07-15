##################################################
## Test file for gets-test-larch-source.R. First
## created 27 May 2024, Oslo.
##
## 1 INITIATE
## 2 TEST larch() ARGUMENTS
## 3 CHECK VARIANCE PREDICTIONS
## 4 EXTRA TESTS OF gets.larch() FUNCTION
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

##remove everything in workspace (.GlobalEnv):
rm(list=ls())

##load source:
source("./gets/R/gets-base-source.R")
source("./gets/R/gets-larch-source.R")


##################################################
## 2 TEST larch() ARGUMENTS
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
eps <- rnorm(n)
eps[1] <- NA; eps[length(eps)] <- NA
eps[c(3,11)] <- 0
vX <- matrix(rnorm(length(eps)*4), length(eps), 4)

##test each argument separately and together:
larch(eps, vc=FALSE) #should return "Error in larch(eps, vc = FALSE) : argument 'vc' must be equal to TRUE"
larch(eps, vc=TRUE)
larch(eps, arch=1)
larch(eps, arch=c(1,3))
larch(eps, harch=1)
larch(eps, harch=c(2,4))
larch(eps, asym=1)
larch(eps, asymind=1)
larch(eps, log.ewma=c(2,4))
larch(eps, vxreg=vX)
larch(eps, zero.adj=1)
larch(eps, vcov.type ="hac") 
larch(eps, qstat.options=c(3,4))
larch(eps, normality.JarqueB=TRUE)
larch(eps, arch=1, tol=1e-15) 
larch(eps, arch=1, tol=1) 
larch(eps, arch=1, tol=1, singular.ok=FALSE) #should return "Error in ols(loge2, x, tol = tol, method = 2) : singular regressor-matrix" 
larch(eps, arch=1, harch=1, singular.ok=FALSE) #should return "Error in ols(loge2, x, tol = tol, method = 2) : singular regressor-matrix" 
larch(eps, arch=1, harch=c(2,4), asym=1, asymind=1, log.ewma=c(3,5),
  vxreg=vX, vcov.type="hac")

##test methods and extraction functions:
mod01 <- larch(eps, arch=1, harch=c(2,4), asym=1, asymind=1, log.ewma=c(3,5),
  vxreg=vX, vcov.type="hac")
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
predict(mod01, newvxreg=matrix(0,12,NCOL(vX)))
residuals(mod01)
summary(mod01)
toLatex(mod01)
vcov(mod01)

##test predict.larch() on each argument:
predict(larch(eps))
predict(larch(eps, arch=1))
predict(larch(eps, harch=1))
predict(larch(eps, asym=1))
predict(larch(eps, asymind=1))
predict(larch(eps, log.ewma=5))
predict(larch(eps, vxreg=vX), newvxreg=matrix(0,12,NCOL(vX)))

##test predict.larch() on each argument with n.ahead=1:
predict(larch(eps), n.ahead=1)
predict(larch(eps, arch=1), n.ahead=1)
predict(larch(eps, harch=1), n.ahead=1)
predict(larch(eps, asym=1), n.ahead=1)
predict(larch(eps, asymind=1), n.ahead=1)
predict(larch(eps, log.ewma=5), n.ahead=1)
predict(larch(eps, vxreg=vX), newvxreg=matrix(0,1,NCOL(vX)), n.ahead=1)


##################################################
## 3 CHECK VARIANCE PREDICTIONS
##################################################

##generate some data:
##===================

##small dgp:
set.seed(123)
vY <- rnorm(20)

##arch(0) model w/constant:
##=========================

mymodel <- larch(vY)

##predictions of the variance:
functionVals <- predict(mymodel, n.ahead=3, n.sim=2)

##correct predictions:
sd2hat1 <- sd2hat2 <- sd2hat3 <- exp(coef(mymodel))
correctVals <- c(sd2hat1,sd2hat2,sd2hat3)

##do they correspond?:
all( functionVals == correctVals )


##################################################
## 4 EXTRA TESTS OF gets.larch() FUNCTION
##################################################

##test extraction functions:
mod01 <- larch(eps, arch=1, harch=c(2,4), asym=1, asymind=1, log.ewma=c(3,5),
  vxreg=vX, vcov.type="hac")
spec01 <- gets(mod01)
print(spec01)
print(spec01, signif.stars=FALSE)
coef(spec01)
fitted(spec01)
gets(spec01)
logLik(spec01)
model.matrix(spec01)
nobs(spec01)
predict(spec01) #should return error: "Error in predict.larch(spec01) : 'newvxreg' is NULL"
predict(spec01, newvxreg=matrix(0,12,2))
residuals(spec01)
summary(spec01)
toLatex(spec01)
vcov(spec01)

##test predict.larch():
predict(spec01, newvxreg=matrix(0,1,2), n.ahead=1)
predict(spec01, newvxreg=matrix(0,12,2), newindex=11:22)
predict(spec01, newvxreg=matrix(0,12,2), n.sim=1)
predict(spec01, newvxreg=matrix(0,12,2), n.sim=1, innov=rnorm(12))
predict(spec01, newvxreg=matrix(0,12,2), probs=c(0.25,0.75))
predict(spec01, newvxreg=matrix(0,12,2), quantile.type=1, probs=c(0.25,0.75))
predict(spec01, newvxreg=matrix(0,12,2), n.sim=10, verbose=TRUE)
predict(spec01, newvxreg=matrix(0,12,2), n.sim=10, probs=c(0.25,0.75),
  verbose=TRUE)

##test 1-cut modelling (i.e. max.paths=0):
specific <- gets(mod01, max.paths=0, do.pet=FALSE)
print(specific)
specific <- gets(mod01, max.paths=0, do.pet=TRUE)
print(specific)


##################################################
## 5 SIMULATIONS (FOR THE FUTURE)
##################################################

TBA