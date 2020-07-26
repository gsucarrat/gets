  ####################################################
## This file contains a prototype of the gets.plm()
## function.
##
## Created 8 July 2020 by Genaro Sucarrat.
##
## CONTENTS:
##
## 1 INITIATE
## 2 CREATE ARTIFICIAL DATA
## 3 PROTOTYPE OF gets.plm()
## 4 CODE EXAMPLES
##                       
## Workflow: First estimate a model with plm(), then
## apply GETS to it. Example:
##
## mygum <- plm(y ~ x1 + x2 + x3)
## gets(mygum)
##
####################################################


####################################################
## 1 INITIATE
####################################################

##required packages:
library(gets)
library(plm)

##clean workspace:
rm(list = ls())


###########################################################
## 2 CREATE ARTIFICIAL DATA
###########################################################

##create some artificial data:
##============================

iN <- 20 #no. of firms
iT <- 4 #no. of time periods (e.g. year)
iNiT <- iN*iT
set.seed(123)
Z <- rnorm(iNiT)
x <- matrix(rnorm(iNiT*10), iNiT, 10)
colnames(x) <- letters[1:10]
firm <- as.vector( t( 1:iN*matrix(rep(1,iNiT), iN, iT) ) )
year <- rep(2001:2004, iN)
mydata <- data.frame(firm, year, Z, x)
head(mydata)

##delete unnecessary stuff from workspace:
rm(iN, iT, iNiT, Z, x, firm, year)

##check that plm works on data:
##=============================

object <- 
  plm(Z ~ a + b + c + d + e + f + g + h + i + j,
  data=mydata)
summary(object)


###########################################################
## 3 PROTOTYPE OF gets.plm()
###########################################################

gets.plm <- function(object, ...)
{

  callAsList <- as.list(object$call)
  listOfArgs <- callAsList[-1]
  plmdata <- as.data.frame( get(as.character(callAsList$data)) )
  idAndTime <- plmdata[,1:2]
  
  ##extract names:
  formulaNames <- all.vars(callAsList$formula)
  yName <- formulaNames[1]
  xNames <- formulaNames[-1]
  
  ##make copies of y and x:
  y <- plmdata[,yName]
  x <- plmdata[,xNames]
  colnames(x) <- paste0("x", 1:NCOL(x))
  
  ##create GUMdata myEstimator:
  GUMdata <- as.data.frame(cbind(idAndTime, y, x))
  
  ##do some clean-up:
  rm(idAndTime)
  
  ##re-define x:
  x <- colnames(x)
  x <- rbind(x,x)
  
  ##create user-specified estimator:
  ##================================
  
  myEstimator <- function(y, x, data=NULL, listOfArgs=NULL)
  {
    ##handle NULL-matrices (obligatory):
    if(is.null(x) || NCOL(x)==0){
  
      result <- list()
    	eps <- data[,"y"]
    	result$logl <- sum(dnorm(eps, sd=sd(eps), log=TRUE))
    	result$n <- NROW(eps)
    	result$k <- 0
    	result$df <- result$n - result$k
  
    }else{ ##if x is not NULL:
  
      ##make formula:
      myformula <- paste0("y ~ ", x[1,])
      if(NCOL(x)>1){
        for(i in 2:NCOL(x)){
          myformula <- paste0(myformula, " + ", x[1,i])
        }
      }
      myformula <- as.formula(myformula)

      ##estimate:
      listOfArgs$formula <- myformula
      listOfArgs$data <- GUMdata
      tmp <- do.call("plm", listOfArgs)
  
    	#rename and re-organise:
     	result <- list()
    	result$coefficients <- coef(tmp)
    	result$vcov <- vcov(tmp)
    	eps <- residuals(tmp)
    	result$logl <- sum(dnorm(eps, sd=sd(eps), log=TRUE))
    	result$n <- NROW(eps)
    	result$k <- NCOL(x)
    	result$df <- result$n - result$k
  
    }
  
    ##final output:
    return(result)

  } #close myEstimator()

  ##do gets:
  ##========

  ##do the gets:
  getsResult <- getsFun(y, x,
    user.estimator=list(name="myEstimator", data=GUMdata,
    listOfArgs=listOfArgs, envir=environment()), ...)
##for testing purposes:
##    listOfArgs=listOfArgs, envir=environment()), keep=5)
  bestTerminal <- getsResult$terminals[[ getsResult$best.terminal ]]
  getsResult$call <- NULL

  ##final model:
  ##============
  
  plmResult <- NULL
  plmClass <- NULL
  
  ##if empty:
  if( length(bestTerminal)==0 ){ cat("\n   Final model is empty \n\n") }
  
  ##if non-empty (estimate final model):
  if( length(bestTerminal)>0 ){

    ##create formula:
    bestXs <- xNames[bestTerminal]    
    myformula <- paste0(yName, " ~ ", bestXs)
    if( length(bestTerminal)>1 ){
      for(i in 2:length(bestTerminal)){
        myformula <- paste0(myformula, " + ", bestXs[i])
      }
    }
    myformula <- as.formula(myformula)
  
    ##estimate final model:
    listOfArgs$formula <- myformula
    plmResult <- do.call("plm", listOfArgs)

    ##extract class:
    plmClass <- class(plmResult)
  }
  
  ##return result:
  ##==============
  
  result <- c(getsResult, plmResult)
  if(!is.null(plmClass)){ class(result) <- plmClass }
  return(result)
  
} #close gets.plm


###########################################################
## 4 CODE EXAMPLES
###########################################################

##estimate gum, do gets:
##======================

mygum <- 
  plm(Z ~ a + b + c + d + e + f + g + h + i + j,
  data=mydata)
summary(mygum)

myspecific <- gets(mygum) #101 estimations
summary(myspecific)

myspecific <- gets(mygum, turbo=TRUE) #56 estimations
summary(myspecific)

myspecific <- gets(mygum, keep=2)
summary(myspecific)

myspecific <- gets(mygum, t.pval=0.4)
summary(myspecific)

##new gum, do gets:
##=================

mygum <- 
  plm(Z ~ a + b + c + d + e + f + g + h + i + j,
  data=mydata, effect="twoways")
summary(mygum)

myspecific <- gets(mygum) #101 estimations
summary(myspecific)

myspecific <- gets(mygum, turbo=TRUE) #56 estimations
summary(myspecific)

myspecific <- gets(mygum, keep=2)
summary(myspecific)

myspecific <- gets(mygum, t.pval=0.4)
summary(myspecific)

##new gum, do gets:
##=================

mygum <- 
  plm(Z ~ a + b + c + d + e + f + g + h + i + j,
  data=mydata, model="random")
summary(mygum)

myspecific <- gets(mygum) #101 estimations
summary(myspecific)

myspecific <- gets(mygum, turbo=TRUE) #56 estimations
summary(myspecific)

myspecific <- gets(mygum, keep=2)
summary(myspecific)

myspecific <- gets(mygum, t.pval=0.4)
summary(myspecific)
