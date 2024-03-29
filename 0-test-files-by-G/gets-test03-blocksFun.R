##################################################
## Test file for gets package. First created
## 22 July 2020, Zaragoza.
##
## 1 INITIATE
## 2 TEST MAIN getsFun ARGUMENTS
## 3 TEST SPECIAL CASES
## 7 SIMULATIONS (FOR THE FUTURE)
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
source("./gets/R/gets-base-source.R")
source("./gets/R/gets-isat-source.R")


##################################################
## 2 TEST MAIN blocksFun ARGUMENTS
##################################################

##dgp1:
set.seed(123)
vY <- rnorm(20)
mX <- matrix(rnorm(length(vY)*10), length(vY), 10)
vY <- mX[,10] + 0.1*rnorm(NROW(mX))

##note: in the tests that follows, it seems the name of the
##regressor matrix, mX, is used unnecessarily in naming objects
##of the value returned by blocksFun(), check?

##check if each argument works:
blocksFun(vY, mX)$specific.spec
blocksFun(vY, mX,
  untransformed.residuals=NA)$specific.spec #should have no effect
mXX <- list(); mXX[[1]] <- mX[,1:5]; mXX[[2]] <- mX[,6:10]
blocksFun(vY, mXX)$specific.spec
myblocks <- list(); myblocks[[1]] <- list(1:3, 4:6, 7:10)
blocksFun(vY, mX, blocks=myblocks)$blocks
blocksFun(vY, mX, no.of.blocks=4)$blocks
blocksFun(vY, mX, max.block.size=2)$blocks
blocksFun(vY, mX, ratio.threshold=0.4)$blocks
blocksFun(vY, mX, gets.of.union=FALSE)$specific.spec
blocksFun(vY, mX, force.invertibility=TRUE)$specific.spec
blocksFun(vY, mX,
  user.estimator=list(name="ols", method=4))$specific.spec
blocksFun(vY, mX,
  user.estimator=list(name="ols", method=5))$specific.spec
blocksFun(vY, mX, t.pval=0.5)$specific.spec
blocksFun(vY, mX, wald.pval=0.5)$specific.spec
blocksFun(vY, mX, do.pet=FALSE)$specific.spec
blocksFun(vY, mX, ar.LjungB=c(1,0.99))$specific.spec
blocksFun(vY, mX, arch.LjungB=c(1,0.99))$specific.spec
blocksFun(vY, mX, normality.JarqueB=0.99)$specific.spec
blocksFun(vY, mX, keep=1)$specific.spec
blocksFun(vY, mX, include.gum=TRUE)$specific.spec
blocksFun(vY, mX, include.1cut=TRUE)$specific.spec
blocksFun(vY, mX, include.empty=TRUE)$specific.spec
blocksFun(vY, mX,
  include.gum=TRUE, include.1cut=TRUE, include.empty=TRUE)$specific.spec
blocksFun(vY, mX, max.paths=1)$specific.spec
blocksFun(vY, mX, parallel.options=2)$specific.spec
blocksFun(vY, mX, parallel.options=9)$specific.spec #should return error if 9 > detectCores()
blocksFun(vY, mX, turbo=TRUE)$specific.spec
blocksFun(vY, mX, force.invertibility=FALSE)$specific.spec
blocksFun(vY, mX, tol=1)$specific.spec
blocksFun(vY, mX, LAPACK=TRUE)$specific.spec
blocksFun(vY, mX, max.regs=1)$specific.spec #should give error
blocksFun(vY, mX, print.searchinfo=FALSE)$specific.spec
blocksFun(vY, mX, alarm=TRUE)$specific.spec


##################################################
## 3 TEST SPECIAL CASES
##################################################

##a dgp:
set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(length(y)*10), length(y), 10)
mXX <- list()
mXX[[1]] <- x[,1:5]
mXX[[2]] <- x[,6:10]

##list without names:
blocksFun(y, mXX)

##a single block:
myblocks <- list(); myblocks[[1]] <- list(1:10)
blocksFun(y, x, blocks=myblocks)

##keep:
blocksFun(y, mXX, keep=1)
mykeep <- list(); mykeep[[1]] <- 2; mykeep[[2]] <- 3
blocksFun(y, mXX, keep=mykeep)

##compare isat() and blocksFun():
y <- rnorm(50)
mconst <- rep(1,length(y))
x <- list(IIS=cbind(mconst,coredata(iim(y))),
  SIS=cbind(mconst,coredata(sim(y))))
blocksFun(y, x, keep=list(1,1))
isat(y, iis=TRUE, sis=TRUE)
##compare speed (not much difference really...):
system.time(blocksFun(y,x))
system.time(isat(y, iis=TRUE, sis=TRUE))


##################################################
## 7 SIMULATIONS (FOR THE FUTURE)
##################################################
