##################################################
## This file tests some of the base functions of
## the gets package.
##
## First created 23 September 2014, Oslo.
##
## 1 INITIATE
## 2 TEST ols()
## 3 TEST gmm()
## 4 TEST diagnostics()
## 5 TEST eqwma() AND leqwma()
## 6 TEST regressorsMean()
## 7 TEST regressorsVariance()
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

##remove everything in workspace (.GlobaleEnv):
rm(list=ls())

##load source:
source("./gets/R/gets-base-source.R")

##load library used for some of the tests:
library(testthat)
library(microbenchmark)


##################################################
##2 TEST ols()
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

##print mCoefs:
rownames(mCoefs) <- paste("method ",1:NROW(mCoefs), ":", sep="")
mCoefs

##check log-likelihoods:
test_that("Log-Likelihoods are the same",{
  tmp <- ols(vY, mX, method=3)
  expect_equal(tmp$logl,sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)))
  tmp <- ols(vY, mX, method=4)
  expect_equal(tmp$logl,sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)))
  tmp <- ols(vY, mX, method=5)
  expect_equal(tmp$logl,sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)))
})

##check variance.spec argument:
vX <- log(mX^2)
tmp <- ols(vY, mX, method=3, variance.spec=list(vc=TRUE, arch=1, asym=1,
  log.ewma=2, vxreg=vX))
test_that("Log-Likelihoods now differ with variance.spec argument",{
  expect_false(tmp$logl == sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)))
})

##check that length(y)!=NROW(vxreg) fails:
test_that("check that length(y)!=NROW(vxreg) fails",{
  expect_error(ols(vY, mX, method=3, variance.spec=list(vxreg=vX[-1,])))
})
  
## until version 0.9, this example failed for
## methods 2-5 (example by F-bear, see his
## email 14/11-2016):
set.seed(123)
y <- rnorm(100, 0, 1)
x <- rnorm(100, 0, 1)*10^8
ols(y,cbind(1,x), method=2, LAPACK=FALSE)$coefficients
ols(y,cbind(1,x), method=2, LAPACK=TRUE)$coefficients

##verify that method=0 returns error:
test_that("verify that ols with method=0 yields error",{
  expect_error(ols(vY, mX, method=0))
})

##verify that method=0 returns the info "method = 0 has been deprecated":
ols(vY, mX, method=0)

##verify that singularity fails:
ols(vY, mX, method=1, tol=1) #Should return "Error in ols(vY, mX, method = 1, tol = 1) : singular regressor-matrix" 
ols(vY, mX, method=2, tol=1) #Should return "Error in..."
ols(vY, mX, method=3, tol=1) #Should return "Error in..."
ols(vY, mX, method=4, tol=1) #Should return "Error in..."
ols(vY, mX, method=5, tol=1) #Should return "Error in..."
ols(vY, mX, method=6, tol=1) #Should return "Error in..."

##view the return values of methods 1 to 6:
names(ols(vY, mX, method=1))
names(ols(vY, mX, method=2))
names(ols(vY, mX, method=3))
names(ols(vY, mX, method=4))
names(ols(vY, mX, method=5))
names(ols(vY, mX, method=6))


##speed comparisons:
##==================

##20 observations:
##----------------

set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(20*5),20,5)

microbenchmark(
  ols(y,x,method=1),
  ols(y,x,method=2),
  ols(y,x,method=3),
  ols(y,x,method=4),
  ols(y,x,method=5),
  .lm.fit(x,y),
  lm.fit(x,y),
  lsfit(x,y),
  lm(y~x-1)
)

#Unit: microseconds
#                  expr   min     lq    mean median     uq   max neval
# ols(y, x, method = 1)   5.4   7.20   9.976   9.70  11.60  21.3   100
# ols(y, x, method = 2)   9.9  13.30  17.981  16.40  20.30  63.4   100
# ols(y, x, method = 3)  18.5  26.70  30.431  29.10  34.15  48.3   100
# ols(y, x, method = 4)  23.8  30.00  35.135  34.50  39.35  49.8   100
# ols(y, x, method = 5)  42.0  54.00  59.820  59.40  65.50  83.3   100
#         .lm.fit(x, y)   3.1   4.55   6.410   5.45   7.45  16.9   100
#          lm.fit(x, y)  27.1  37.70  41.910  41.40  45.70  62.3   100
#           lsfit(x, y)  53.9  68.65  77.486  73.65  79.20 368.3   100
#         lm(y ~ x - 1) 508.3 525.00 566.123 537.15 610.70 818.4   100


##1000 observations, many variables:
##----------------------------------

set.seed(123)
y <- rnorm(1000)
x <- matrix(rnorm(1000*30),1000,30)

microbenchmark(
  ols(y,x,method=1),
  ols(y,x,method=2),
  ols(y,x,method=3),
  ols(y,x,method=4),
  ols(y,x,method=5),
  .lm.fit(x,y),
  lm.fit(x,y),
  lsfit(x,y),
  lm(y~x-1)
)

#Unit: microseconds
#                  expr     min       lq      mean   median       uq     max neval
# ols(y, x, method = 1)   766.1   834.75   855.626   848.15   866.80  1011.3   100
# ols(y, x, method = 2)   828.0   874.55   899.553   891.05   916.45  1087.2   100
# ols(y, x, method = 3)   824.6   895.90   980.789   915.20   936.05  6795.2   100
# ols(y, x, method = 4)  1771.6  1902.05  2003.275  1945.90  1995.35  7780.7   100
# ols(y, x, method = 5) 12865.7 14107.20 14676.436 14242.95 14466.30 20852.7   100
#         .lm.fit(x, y)   770.1   826.60   903.578   839.30   863.65  6660.1   100
#          lm.fit(x, y)   812.4   884.75   915.807   910.10   941.90  1076.6   100
#           lsfit(x, y)  1021.6  1146.55  1303.762  1181.55  1211.10  7385.7   100
#         lm(y ~ x - 1)  1661.9  2030.95  2234.919  2134.30  2210.65  8359.4   100


##################################################
##3 TEST gmm()
##################################################

##in this experiment the aim is simply to verify
##that the arguments work
##----------------------------------------------

set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(20*2), 20, 2)
z <- matrix(rnorm(20*2), 20, 2)

gmm(y,x,z) #basic test
gmm(y,x,z, tol=1) ##should return error
gmm(y,x,z, tol=0.01) ##should return error
gmm(y,x,z, tol=0.001) ##should work
gmm(y,x,z, weighting.matrix="efficient", vcov.type="ordinary")
gmm(y,x,z, weighting.matrix="efficient", vcov.type="robust")
gmm(y,x,z, weighting.matrix="2sls", vcov.type="ordinary")
gmm(y,x,z, weighting.matrix="2sls", vcov.type="robust")
gmm(y,x,z, weighting.matrix="identity", vcov.type="ordinary")
gmm(y,x,z, weighting.matrix="identity", vcov.type="robust")


##in this experiment NCOL(x) < NCOL(z), so there
##are more instruments than regressors; the aim is
##to check whether 2sls and efficient gmm work.
##------------------------------------------------

set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(20*2), 20, 2)
z <- matrix(rnorm(20*3), 20, 3)

gmm(y,x,z, weighting.matrix="efficient", vcov.type="ordinary")
gmm(y,x,z, weighting.matrix="efficient", vcov.type="robust")
gmm(y,x,z, weighting.matrix="2sls", vcov.type="ordinary")
gmm(y,x,z, weighting.matrix="2sls", vcov.type="robust")

##in this experiment x=z, so the results from
##ols() and gmm() should be identical
##-------------------------------------------

set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(20*2), 20, 2)
z <- x

##vcov.type="ordinary":
olsEst <- ols(y,x, method=3)
ivEst <- gmm(y,x,x, weighting.matrix="identity", vcov.type="ordinary")
round(ivEst$coefficients, digits=10) == round(olsEst$coefficients, digits=10)
all( round(ivEst$fit, digits=10) == round(olsEst$fit, digits=10) )
all( round(ivEst$residuals, digits=10) == round(olsEst$residuals, digits=10))
round(ivEst$rss, digits=10) == round(olsEst$rss, digits=10)
round(ivEst$sigma2, digits=10) == round(olsEst$sigma2, digits=10)
round(ivEst$vcov, digits=10) == round(olsEst$vcov, digits=10)
round(ivEst$logl, digits=10) == round(olsEst$logl, digits=10)

##vcov.type="robust":
olsEst <- ols(y,x, method=4)
ivEst <- gmm(y,x,x, weighting.matrix="identity", vcov.type="robust")
round(ivEst$coefficients, digits=10) == round(olsEst$coefficients, digits=10)
all( round(ivEst$fit, digits=10) == round(olsEst$fit, digits=10) )
all( round(ivEst$residuals, digits=10) == round(olsEst$residuals, digits=10))
round(ivEst$rss, digits=10) == round(olsEst$rss, digits=10)
round(ivEst$sigma2, digits=10) == round(olsEst$sigma2, digits=10)
round(ivEst$vcov, digits=10) == round(olsEst$vcov, digits=10)
round(ivEst$logl, digits=10) == round(olsEst$logl, digits=10)


##in this experiment NCOL(x)==NCOL(z), but x!=z;
##also, the regressor x is correlated with the
##error; the aim is to obtain an idea of whether
##estimates are consistent
##------------------------------------------------

set.seed(123)
n <- 10000
z1 <- rnorm(n)
eps <- rnorm(n) #ensures cor(z,eps)=0
x1 <- 0.5*z1 + 0.5*eps #ensures cor(x,eps) is strong
y <- 0.4 + 0.8*x1 + eps #the dgp
cor(x1, eps) #check correlatedness
cor(z1, eps) #check uncorrelatedness

x <- cbind(1,x1) #regressor matrix
z <- cbind(1,z1) #instrument matrix

##ols as benchmark (should be inconsistent):
tmp <- ols(y,x)
tmp$coefficients #true values are 0.4. and 0.8

##check consistency of iv estimator:
tmp <- gmm(y,x,z, weighting.matrix="identity")
tmp$coefficients #should be approx 0.4 and 0.8
tmp$vcov
tmp <- gmm(y,x,z, weighting.matrix="identity", vcov.type="robust")
tmp$vcov

##check consistency of 2sls estimator:
tmp <- gmm(y,x,z, weighting.matrix="2sls")
tmp$coefficients #should be approx 0.4 and 0.8
tmp$vcov
tmp <- gmm(y,x,z, weighting.matrix="2sls", vcov.type="robust")
tmp$vcov

##check consistency of efficient gmm estimator:
tmp <- gmm(y,x,z, weighting.matrix="efficient")
tmp$coefficients #should be approx 0.4 and 0.8
tmp$vcov
tmp <- gmm(y,x,z, weighting.matrix="efficient", vcov.type="robust")
tmp$vcov

##compare log-likelihoods (should return TRUE):
tmp <- gmm(y,x,z, weighting.matrix="identity")
round(tmp$logl, digits=10) ==
  round(sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)),
  digits=10)

##compare log-likelihoods (should return TRUE):
tmp <- gmm(y,x,z, weighting.matrix="2sls")
round(tmp$logl, digits=10) ==
  round(sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)),
  digits=10)

##compare log-likelihoods (should return TRUE):
tmp <- gmm(y,x,z, weighting.matrix="efficient")
round(tmp$logl, digits=10) ==
  round(sum(dnorm(tmp$residuals, sd=sqrt(tmp$sigma2), log=TRUE)),
  digits=10)


##the aim of this experiment is to test whether
##gmm() works with getsFun() and blocksFun()
##------------------------------------------------

set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(20*10), 20, 10)
colnames(x) <- paste0("x",1:10)
z <- matrix(rnorm(20*10), 20, 10)

getsFun(y,x, user.estimator=list(name="gmm", z=z))
getsFun(y,x, user.estimator=list(name="gmm", z=z),
  keep=6)

blocksFun(y,x, user.estimator=list(name="gmm", z=z))
blocksFun(y,x, user.estimator=list(name="gmm", z=z),
  keep=6)


##################################################
##4 TEST diagnostics()
##################################################

set.seed(123)
vY <- rnorm(20)
mX <- matrix(rnorm(3*20), 20, 3)

##check x for autocorrelation and ARCH, and return a
##data-frame with the results:
x <- ols(vY, mX, method=3)
diagnostics(x)

##check x for autocorrelation and ARCH, and indicate
##whether it passes the check:
diagnostics(x, verbose=FALSE) ##should be TRUE
diagnostics(x, ar.LjungB=c(1,0.96), verbose=FALSE) ##should be FALSE
diagnostics(x, arch.LjungB=c(1,0.10), verbose=FALSE) ##should be FALSE

##add the Jarque-Bera normality test to the diagnostics:
diagnostics(x, normality.JarqueB=TRUE)
diagnostics(x, normality.JarqueB=0.8) #should have no effect
diagnostics(x, normality.JarqueB=0.8, verbose=FALSE) #should return TRUE
diagnostics(x, normality.JarqueB=0.9, verbose=FALSE) #should return FALSE

##use test_that to check:
test_that("Test the Diagnostics matrix",{
  expect_true(is.matrix(diagnostics(x)))
  
  ##check x for autocorrelation and ARCH, and indicate
  ##whether it passes the check:
  expect_true(diagnostics(x, verbose=FALSE))
  expect_false(diagnostics(x, ar.LjungB=c(1,0.96), verbose=FALSE))
  expect_false(diagnostics(x, arch.LjungB=c(1,0.10), verbose=FALSE))
  
  ##add the Jarque-Bera normality test to the diagnostics:
  expect_true(is.matrix(diagnostics(x, normality.JarqueB=TRUE)) & nrow(diagnostics(x, normality.JarqueB=TRUE))==3)
  expect_true(is.matrix(diagnostics(x, normality.JarqueB=0.9)) & nrow(diagnostics(x, normality.JarqueB=TRUE))==3) #should have no effect
  expect_false(diagnostics(x, normality.JarqueB=0.9, verbose=FALSE))
  
})

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}
SWtest(x)
diagnostics(x, user.fun=list(name="SWtest", pval=0.025))
diagnostics(x, user.fun=list(name="SWtest", pval=0.025),
  verbose=FALSE) #should return TRUE
diagnostics(x, user.fun=list(name="SWtest", pval=0.85),
  verbose=FALSE) #should return FALSE

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
diagnostics(x, user.fun=list(name="SWtest", pval=0.025), verbose=FALSE) #should return TRUE
diagnostics(x, user.fun=list(name="SWtest", pval=0.85), verbose=FALSE) #should return FALSE

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
diagnostics(x, user.fun=list(name="SWtest")) #should not work ("...could not find...")
diagnostics(x, user.fun=list(name="SWtest", envir=myenv)) #should work
rm("myenv") #clean up

##check whether variance.spec works (creates NAs in std.residuals):
x <- ols(vY, mX, method=3, variance.spec=list(vc=TRUE, arch=1))
diagnostics(x)
x <- ols(vY, mX, method=3, variance.spec=list(vc=TRUE, arch=1))
test_that("check whether variance.spec works (creates NAs in std.residuals)",{
  expect_silent(diagnostics(x))
})

##check no. 1 of "is.rejection.bad" entry:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
  rownames(result) <- "Shapiro-test"
  return(result)
}
diagnostics(x, user.fun=list(name="SWtest", pval=0.025))
diagnostics(x,
  user.fun=list(name="SWtest", pval=0.025),
  verbose=FALSE) #should return TRUE
diagnostics(x,
  user.fun=list(name="SWtest", pval=0.85),
  verbose=FALSE) #should return FALSE
diagnostics(x,
  user.fun=list(name="SWtest", pval=0.85, is.reject.bad=TRUE),
  verbose=FALSE) #should return FALSE
diagnostics(x,
  user.fun=list(name="SWtest", pval=0.85, is.reject.bad=FALSE),
  verbose=FALSE) #should return TRUE

##check no. 2 of "is.rejection.bad" entry:
SWtest <- function(x, ...){
  result <- matrix(NA, 2, 3)
  tmp <- shapiro.test(x$residuals)
  result[1,] <- c(tmp$statistic, NA, tmp$p.value)
  tmp <- shapiro.test(x$residuals^2)
  result[2,] <- c(tmp$statistic, NA, tmp$p.value)
  rownames(result) <- c("Test 1", "Test 2")
  return(result)
}
diagnostics(x, user.fun=list(name="SWtest", pval=0.025))
diagnostics(x,
  user.fun=list(name="SWtest", pval=0.025),
  verbose=FALSE) #should return FALSE
##SHOULD THIS ONE RETURN FALSE??:
diagnostics(x,
  user.fun=list(name="SWtest", pval=0.85, is.reject.bad=c(TRUE,FALSE)),
  verbose=FALSE) #should return TRUE


##################################################
## 5 TEST eqwma AND leqwma
##################################################

set.seed(123)
x <- rnorm(10)

##test eqwma:
##===========

eqwma(x)
eqwma(as.zoo(x))
eqwma(x, length=2)
eqwma(x, length=c(2,3))
eqwma(x, k=2 )
eqwma(x, p=2)
eqwma(x, abs=TRUE)
eqwma(x, log=TRUE)
eqwma(x, as.vector=TRUE)
eqwma(x, start=1) #should return error: 'start' has been deprecated
eqwma(x, lag=1) #should return error: 'lag' has been deprecated, use 'k' instead

test_that("test EQMA",{
  expect_silent(eqwma(x))
  expect_silent(eqwma(as.zoo(x)))
  expect_silent(eqwma(x, length=2))
  expect_silent(eqwma(x, length=c(2,3)))
  expect_silent(eqwma(x, k=2 ))
  expect_silent(eqwma(x, p=2))
  expect_silent(eqwma(x, abs=TRUE))
  expect_silent(eqwma(x, log=TRUE))
  expect_silent(eqwma(x, as.vector=TRUE))
  
  expect_error(eqwma(x, start=1)) #should return error
  expect_error(eqwma(x, lag=1)) #should return error
})

##test leqwma:
##============

leqwma(x)
leqwma(as.zoo(x))
leqwma(x, length=2)
leqwma(x, length=c(2,3))
leqwma(x, k=2 )
leqwma(x, p=1)
leqwma(x, as.vector=TRUE)
leqwma(x, start=1) #should return error: 'start' has been deprecated
leqwma(x, lag=1) #should return error: 'lag' has been deprecated, use 'k' instead

test_that("test leqwma",{
  expect_silent(leqwma(x))
  expect_silent(leqwma(as.zoo(x)))
  expect_silent(leqwma(x, length=2))
  expect_silent(leqwma(x, length=c(2,3)))
  expect_silent(leqwma(x, k=2 ))
  expect_silent(leqwma(x, p=1))
  expect_silent(leqwma(x, as.vector=TRUE))
  expect_error(leqwma(x, start=1)) #should return error
  expect_error(leqwma(x, lag=1)) #should return error
})


##################################################
## 6 TEST regressorsMean()
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

##for visual inspection:
regressorsMean(y)
regressorsMean(log(y^2))
regressorsMean(y, mc=TRUE)
regressorsMean(y, ar=c(1,3))
regressorsMean(y, ewma=list(length=c(2,4)))
regressorsMean(y, mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, prefix="G")
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.as.zoo=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, na.trim=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE, return.as.zoo=FALSE,
  na.trim=FALSE)
##erroneous until version 0.23:
colnames(mxreg) <- c("a", "", "c", "", "e")
regressorsMean(y, mxreg=mxreg)
regressorsMean(y, mxreg=mxreg, prefix="")
 
##unit test:
test_that("Testing the regressorsMean() function - Test 1",{
  expect_silent(regressorsMean(y))
  expect_silent(regressorsMean(log(y^2)))
  expect_silent(regressorsMean(y, mc=TRUE))
  expect_silent(regressorsMean(y, ar=c(1,3)))
  expect_silent(regressorsMean(y, ewma=list(length=c(2,4))))
  expect_silent(regressorsMean(y, mxreg=mxreg))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, return.regressand=FALSE))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, return.as.zoo=FALSE))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, na.trim=FALSE))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, return.regressand=FALSE, return.as.zoo=FALSE,na.trim=FALSE))
  
  
  ##erroneous until version 0.23:
  colnames(mxreg) <- c("a", "", "c", "", "e")
  expect_silent(regressorsMean(y, mxreg=mxreg))
})

##test 2:
##=======

set.seed(123)
iT <- 10 #60 or 100. If iT=60, then usually no specific
y <- arima.sim(list(ar=0.3),iT)
y <- ts(y, frequency=4, end=c(2015,4))
mxreg <- matrix(rnorm(4*iT), iT, 4)
mxreg <- ts(mxreg, frequency=4, end=c(2015,4))
y[1] <- NA; y[iT] <- NA

##for visual inspection:
regressorsMean(y)
regressorsMean(log(y^2))
regressorsMean(y, mc=TRUE)
regressorsMean(y, ar=c(1,3))
regressorsMean(y, ewma=list(length=c(2,4)))
regressorsMean(y, mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, prefix="G")
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.as.zoo=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, na.trim=FALSE)
regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),
  mxreg=mxreg, return.regressand=FALSE, return.as.zoo=FALSE,
  na.trim=FALSE)
##used to yield error:
set.seed(123)
vY <- rnorm(20)
regressorsMean(vY, mc=TRUE, ar=1, ewma=list(length=2))

##unit test:
test_that("Testing the regressorsMean() function - Test 2",{
  expect_silent(regressorsMean(y))
  expect_silent(regressorsMean(log(y^2)))
  expect_silent(regressorsMean(y, mc=TRUE))
  expect_silent(regressorsMean(y, ar=c(1,3)))
  expect_silent(regressorsMean(y, ewma=list(length=c(2,4))))
  expect_silent(regressorsMean(y, mxreg=mxreg))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, return.regressand=FALSE))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, return.as.zoo=FALSE))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, na.trim=FALSE))
  expect_silent(regressorsMean(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)),mxreg=mxreg, return.regressand=FALSE, return.as.zoo=FALSE,na.trim=FALSE))
  ##used to yield error:
  set.seed(123)
  vY <- rnorm(20)
  expect_silent(regressorsMean(vY, mc=TRUE, ar=1, ewma=list(length=2)))
  
})

##test 3:
##=======

##from 0.24: mxreg can be a data.frame
set.seed(123)
y <- rnorm(10)
mxreg <- matrix(rnorm(10*5), 10, 5)
mxreg <- as.data.frame(mxreg)

regressorsMean(y, mxreg=mxreg)

test_that("Testing the regressorsMean() function - Test 3 on Dataframe",{
  expect_silent(regressorsMean(y, mxreg=mxreg))
})


##################################################
## 7 TEST regressorsVariance()
##################################################

##test 1:
##=======

set.seed(123)
iT <- 10
eps <- rnorm(iT)
eps <- arima.sim(list(ar=0.3), iT)
eps[1] <- NA; eps[iT] <- NA
eps[3] <- 0
vxreg <- matrix(rnorm(5*iT), iT, 5)
#vxreg <- cbind(rep(1, iT)); colnames(mX) <- "mconst"
#vxreg[1:5,2] <- NA
    
##for visual inspection:
regressorsVariance(eps)
regressorsVariance(eps, vc=FALSE)
regressorsVariance(eps, arch=1)
regressorsVariance(eps, arch=c(2,3))
regressorsVariance(eps, harch=1)
regressorsVariance(eps, harch=c(2,3))
regressorsVariance(eps, asym=c(1))
regressorsVariance(eps, asymind=c(2,4))
regressorsVariance(eps, log.ewma=5)
regressorsVariance(eps, log.ewma=c(2,4))
regressorsVariance(eps, log.ewma=list(length=c(2,4)))
regressorsVariance(eps, vxreg=vxreg)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, prefix="G")
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, return.regressand=FALSE)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, return.as.zoo=FALSE)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, na.trim=FALSE)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, return.regressand=FALSE, return.as.zoo=FALSE,
  na.trim=FALSE)
##check naming when "" in colnames:
colnames(vxreg) <- c("a", "", "c", "", "d")
regressorsVariance(eps, vxreg=vxreg)
##check that unique naming works:
colnames(vxreg) <- paste0("arch", 1:5)
regressorsVariance(eps, arch=1, vxreg=vxreg)


##test 2:
##=======

set.seed(123)
iT <- 10
eps <- arima.sim(list(arch=0.3),iT)
eps <- ts(eps, frequency=4, end=c(2015,4))
vxreg <- matrix(rnorm(4*iT), iT, 4)
vxreg <- ts(vxreg, frequency=4, end=c(2015,4))
eps[1] <- NA; eps[iT] <- NA
eps[3] <- 0

##for visual inspection:
regressorsVariance(eps)
regressorsVariance(eps, vc=FALSE)
regressorsVariance(eps, arch=c(1,3))
regressorsVariance(eps, harch=c(2,4))
regressorsVariance(eps, asym=c(1,3))
regressorsVariance(eps, asymind=c(2,4))
regressorsVariance(eps, log.ewma=c(2,4))
regressorsVariance(eps, vxreg=vxreg)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, return.regressand=FALSE)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, return.as.zoo=FALSE)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, na.trim=FALSE)
regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),
  vxreg=vxreg, return.regressand=FALSE, return.as.zoo=FALSE,
  na.trim=FALSE)

##unit test:
test_that("regressorsVariance() - Test 2",{
  expect_silent(regressorsVariance(eps))
  expect_silent(regressorsVariance(eps, vc=FALSE))
  expect_silent(regressorsVariance(eps, arch=c(1,3)))
  expect_silent(regressorsVariance(eps, log.ewma=c(2,4)))
  expect_silent(regressorsVariance(eps, vxreg=vxreg))
  expect_silent(regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),vxreg=vxreg))
  expect_silent(regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),vxreg=vxreg, return.regressand=FALSE))
  expect_silent(regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),vxreg=vxreg, return.as.zoo=FALSE))
  expect_silent(regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),vxreg=vxreg, na.trim=FALSE))
  expect_silent(regressorsVariance(eps, vc=TRUE, arch=c(1,3), log.ewma=c(2,4),vxreg=vxreg, return.regressand=FALSE, return.as.zoo=FALSE,na.trim=FALSE))
})


##test 3:
##=======

##from 0.24: vxreg can be a data.frame
set.seed(123)
eps <- rnorm(10)
vxreg <- matrix(rnorm(10*5), 10, 5)
vxreg <- as.data.frame(vxreg)

regressorsVariance(eps, vxreg=vxreg)

test_that("regressorsVariance() - Test 3 vxreg as data.frame",{
  expect_silent(regressorsVariance(eps, vxreg=vxreg))
})
