##################################################
## Test file for the lm functions. First created
## 22 July 2021, Zaragoza.
##
## 1 INITIATE
## 2 TEST as.arx()
## 3 TEST gets.lm()
## 3 TEST isat.lm()
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

##remove everything in workspace (.GlobalEnv):
rm(list=ls())

##load source:
source("./gets/R/gets-base-source.R")
source("./gets/R/gets-lm-source.R")


##################################################
## 2 TEST as.arx.lm()
##################################################

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##generate some data:
set.seed(123) #for reproducibility
y <- rnorm(40) #generate Y
x <- matrix(rnorm(40*20), 40, 20) #create matrix of Xs

##typical situation:
mymodel <- lm(y ~ x)
as.arx(mymodel)

##only intercept)
mymodel <- lm(y ~ NULL)
as.arx(mymodel)

##empty model:
mymodel <- lm(y ~ -1)
as.arx(mymodel)

##use hetero-robust vcov:
mymodel <- lm(y ~ x)
as.arx(mymodel, vcov.type="white")

##add ar-dynamics:
mymodel <- lm(y ~ x)
as.arx(mymodel, ar=1:2)

##add log-variance specification:
mymodel <- lm(y ~ x)
as.arx(mymodel, arch=1:2)


##################################################
## 3 TEST gets.lm()
##################################################

##generate some data:
set.seed(123) #for reproducibility
v <- rnorm(30) #generate Y
w <- matrix(rnorm(30*10), 30, 10) #create matrix of Xs
colnames(w) <- paste0("var", 1:NCOL(w))

##basic tests of selected getsFun arguments:
mymodel <- lm(v ~ w)
mymodel
gets(mymodel)
gets(mymodel, t.pval=0.2)
gets(mymodel, t.pval=0.9999)
gets(mymodel, wald.pval=0.2)
gets(mymodel, keep=c(1,5))
gets(mymodel, keep=c(1,5), do.pet=FALSE)
gets(mymodel, print.searchinfo=FALSE)
gets(mymodel, include.1cut=TRUE) 
gets(mymodel, include.empty=TRUE) 
gets(mymodel, max.paths=1)

##test 'subset' argument:
mymodel <- lm(v ~ w, subset=1:20)
mymodel
gets(mymodel)

##test 'weigths' argument:
weights <- 0.1 + runif(length(v))
mymodel <- lm(v ~ w, weights=weights)
mymodel
gets(mymodel, keep=c(1,3,6))


##################################################
## 4 TEST isat.lm()
##################################################

##load 'isat' source:
source("./gets/R/gets-isat-source.R")

##use same data:
mymodel <- lm(v ~ w)
mymodel

##basic test:
isat(mymodel)
