##################################################
## Test file for logit function package. First
## created 8 December 2020, Oslo.
##
## 1 INITIATE
## 2 TEST logitxSim()
## 3 TEST logit()
## 4 TEST logitx()
## 5 TEST gets.logitx()
##
##################################################


##################################################
##1 INITIATE
##################################################

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd(choose.dir()) #interactively

##load required packages:
require(parallel)
require(zoo)

##remove everything in workspace (.GlobalEnv):
rm(list=ls())

##load source:
source("./gets/R/gets-base-source.R")
source("./gets/R/gets-logitx-source.R")


##################################################
## 2 TEST logitxSim()
##################################################

##generate some data:
set.seed(123)
n <- 20
z <- rnorm(n)

##for visual inspection:
logitxSim(n)
logitxSim(n, intercept=1)
logitxSim(n, ar=c(0.2,0.1))
#note: ewma not an argument in logitxSim() yet!
logitxSim(n, xreg=0.5*z)
logitxSim(n, verbose=TRUE)
logitxSim(n, as.zoo=FALSE)
logitxSim(n, intercept=1, ar=c(0.2,0.1), xreg=0.5*z, verbose=TRUE,
  as.zoo=FALSE)

##test the alias dlogitxSim():
dlogitxSim(n)
dlogitxSim(n, intercept=1)
dlogitxSim(n, ar=c(0.2,0.1))
dlogitxSim(n, xreg=0.5*z)
dlogitxSim(n, verbose=TRUE)
dlogitxSim(n, as.zoo=FALSE)
dlogitxSim(n, intercept=1, ar=c(0.2,0.1), xreg=0.5*z, verbose=TRUE,
  as.zoo=FALSE)

    
##################################################
## 3 TEST logit()
##################################################

##check arguments:
##================

set.seed(123)
n <- 20
y <- logitxSim(n, ar=1, as.zoo=FALSE)
z <- rnorm(n)
x <- cbind(1,z)

##for visual inspection:
logit(y, NULL)
logit(y,x)
logit(y,x, initial.values=c(0,0))
logit(y,x, lower=0)$coefficients
logit(y,x, upper=0.5)$coefficients
logit(y,x, method=1)$vcov #should be NULL
logit(y,x, method=2)$vcov
logit(y,x, method=3)$vcov
logit(y,x, method=3)$lag.length
logit(y,x, method=3, lag.length=5)$lag.length
logit(y,x, control=list(trace=1))$iterations
logit(y,x, eps.tol=1) 
logit(y,x, solve.tol=0.99) #should return error

    
##basic checks of consistency:
##============================

set.seed(123)
z <- rnorm(10000)
x <- cbind(1,z)

##experiment 1:
y <- logitxSim(10000, intercept=0.5, xreg=z)
result <- logit(y, x)
result$coefficients #should be approx 0.5 and 1

##experiment 2:
y <- logitxSim(10000, intercept=-1, xreg=2*z)
result <- logit(y, x)
result$coefficients #should be approx -1 and 2

##experiment 3:
y <- logitxSim(10000, intercept=0, xreg=1*z)
result <- logit(y, x)
result$coefficients #should be approx 0 and 1

##experiment 4:
y <- logitxSim(10000, intercept=1, xreg=0*z)
result <- logit(y, x)
result$coefficients #should be approx 1 and 0

##simulations:
##============

## dpg: ar(1)
##-----------

##set simulation parameters:
set.seed(123)
theta1 <- 0
theta2 <- 0.9 #0.5
#rho <- 0.8 #values:0, 0.9, -0.9
iReps <- 1000 #no. of replications
n <- 10000

##matrix to store results in:
resultnames <- c("mconst", "ar1")
resultCoefs <- matrix(NA, iReps, length(resultnames))
colnames(resultCoefs) <- resultnames
resultPvals <- resultCoefs

##loop (takes 2-3 minutes):
for(i in 1:iReps){

  ##simulate from dgp:
  y <- logitxSim(n, intercept=theta1, ar=theta2)
  
  ##estimate and record results:
  mydata <- regressorsMean(y, mc=TRUE, ar=1)
  y <- coredata(mydata[,1])
  x <- coredata(mydata[,-1])
  mymod <- logit(y, x, method=3) #method=3: Newey West (1987)
  resultCoefs[i,] <- mymod$coefficients
  tstat <-
    (mymod$coefficients-c(theta1,theta2))/sqrt(diag(mymod$vcov))
  resultPvals[i,] <- pt(tstat, mymod$df, lower.tail=FALSE) 
  
}

##check consistency:
##------------------

colMeans(resultCoefs)
round(colMeans(resultCoefs), digits=3) == c(theta1,theta2)

##check empirical size:
##---------------------

##10%:
colMeans(resultPvals < 0.1)

##5%:
colMeans(resultPvals < 0.05)

##1%
colMeans(resultPvals < 0.01)


##################################################
## 4 TEST logitx()
##################################################

set.seed(123)
n <- 20
y <- logitxSim(n, ar=1)
z <- rnorm(n)

##test arguments (visual inspection):
logitx(y)
logitx(y, intercept=FALSE)
logitx(y, ar=1)
logitx(y, ewma=list(length=c(2,4)))
logitx(y, ewma=c(2,4))
##should drop EqWMA(2):
logitx(y, ar=c(1,2), ewma=c(2,4)) 
##ISSUE: the z-variable is labelled as "mxreg", but should probably be
##labelled "xreg" or "z":
logitx(y, xreg=z)
logitx(y, vcov.type="ordinary")
logitx(y, vcov.type="robust")
logitx(y, vcov.type="robust", lag.length=5)
logitx(y, ar=1, initial.values=c(0,0))
logitx(y, ar=1, lower=0)
logitx(y, ar=1, upper=0.5)
logitx(y, control=list(trace=1))$coefficients
logitx(y, upper=0.5, control=list(trace=1))$coefficients
logitx(y, eps.tol=1) 
##should return error:
logitx(y, solve.tol=2)

##check methods:
x <- logitx(y, ar=1, xreg=z)
print(x)
coef(x)
fitted(x)
fitted(x, zero.prob=TRUE)
logLik(x)
summary(x)
toLatex(x)
vcov(x)

##test the alias dlogitx:
dlogitx(y)
dlogitx(y, intercept=FALSE)
dlogitx(y, ar=1)
dlogitx(y, ewma=list(length=c(2,4)))
dlogitx(y, xreg=z)
dlogitx(y, vcov.type="ordinary")
dlogitx(y, vcov.type="robust")
dlogitx(y, vcov.type="robust", lag.length=5)


##################################################
## 5 TEST gets.logitx()
##################################################

set.seed(123)
n <- 30
z <- rnorm(n)
y <- logitxSim(n, xreg=1*z)
xreg <- cbind(z, matrix(rnorm(n*9), n, 9))

mymod <- logitx(y, xreg=xreg)
gets(mymod)
gets(mymod, t.pval= 0.2)
gets(mymod, wald.pval = 0.0001)
gets(mymod, do.pet = FALSE)
gets(mymod, keep = 1)
gets(mymod, include.gum = TRUE)
gets(mymod, include.1cut = FALSE)
gets(mymod, include.empty = TRUE)
gets(mymod, max.paths = 3)
gets(mymod, turbo=FALSE)
gets(mymod, print.searchinfo = FALSE)
gets(mymod, plot = TRUE)
gets(mymod, alarm = TRUE)

##a tougher test (many variables):
set.seed(123)
n <- 200
z <- rnorm(n)
y <- logitxSim(n, ar=c(0.3,0.2,0.1), xreg=1*z)
xreg <- cbind(z, matrix(rnorm(n*29), n, 29))

mymod <- logitx(y, ar=1:5, xreg=xreg)
mymod
gets(mymod, max.paths=5)
