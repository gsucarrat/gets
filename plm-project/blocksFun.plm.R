####################################################
## This file contains a draft of how to combine
## blocksFun() with plm()
##
## First created 25 July 2020.
##
## CONTENTS:
##
## 1 INITIATE
## 2 CREATE ARTIFICIAL DATA
## 3 DRAFT OF HOW TO COMBINE blocksFun() WITH plm()
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

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd(choose.dir())


####################################################
## 2 CREATE SOME ARTIFICIAL DATA
####################################################

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
coef(object)
vcov(object)


####################################################
## 3 DRAFT OF HOW TO COMBINE blocksFun() WITH plm()
####################################################

##note: the code/draft is almost identical to the
##gets.plm() function.

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

  ##create GUMdata for myEstimator:
  GUMdata <- as.data.frame(cbind(idAndTime, y, x))

  ##do some clean-up:
  rm(idAndTime)

#  ##G-man idea for how to create 'individual effects' indicators:
#  NROWplmdata <- NROW(plmdata)
#  IDs <- as.character(union(plmdata[,1], NULL))
#  IIS <- matrix(0, NROWplmdata, length(IDs))
#  colnames(IIS) <- paste0("iis", IDs)
#  for(i in 1:length(IDs)){
#    whichRows <- which( plmdata[,1]==IDs[i] )
#    IIS[ whichRows , i] <- 1
#  }
#
#  ##add indicators to GUMdata:
#  GUMdata <- as.data.frame(cbind(GUMdata, IIS))
#
#  ##this example illustrates, I think, that iis1 is automatically
#  ##removed (due to exact colinearity):
#  deleteme <-
#    plm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + iis1,
#    data=GUMdata)
#  summary(deleteme)

  ##re-define x:
  x <- colnames(x)
  x <- rbind(x,x) #x must be a matrix
  colnames(x) <- xNames

  ##create user-specified estimator:
  ##================================

  myEstimator <- function(y, x, data=NULL, listOfArgs=NULL)
  {
    ##handle empty matrices (obligatory):
    if(is.null(x) || NCOL(x)==0){

      result <- list()
    	eps <- data[,"y"]
    	result$logl <- sum(dnorm(eps, sd=sd(eps), log=TRUE))
    	result$n <- NROW(eps)
    	result$k <- 0
    	result$df <- result$n - result$k

    }else{ ##if x is not empty:

      ##make formula:
      
      # replace following lines:
      # myformula <- paste0("y ~ ", x[1,])
      # if(NCOL(x)>1){
      #   for(i in 2:NCOL(x)){
      #     myformula <- paste0(myformula, " + ", x[1,i])
      #   }
      # }
      # myformula <- as.formula(myformula)
      myformula <- paste(x[1,], collapse = " + ")
      myformula <- paste("y ~", myformula, sep = " ")
      myformula <- as.formula(myformula)

      ##estimate:
      listOfArgs$formula <- myformula
      listOfArgs$data <- GUMdata
      tmp <- do.call("plm", listOfArgs)

    	##rename and re-organise:
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

  ##do block-based gets:
  ##====================
  
  ##block-based gets:
  myblocks <- list()
  myblocks[[1]] <- list(1:5, 6:10)
  blocksResult <- blocksFun(y, x, blocks=myblocks,
    user.estimator=list(name="myEstimator", data=GUMdata,
    listOfArgs=listOfArgs, envir=environment()),
    gets.of.union=TRUE) #this needs to be added, when inside a function: , ...)
##note: to ensure variables are kept, add (say) 'keep=1'
##isat(), by contrast, does not work:
##  isat(y, uis=x, user.estimator=list(name="myEstimator", data=GUMdata,
##    listOfArgs=listOfArgs, envir=environment()))
  bestXs <- integer(0)
  for(i in 1:length(blocksResult$specific.spec)){
    bestXs <-
      union(bestXs, blocksResult$specific.spec[[i]])
  }
  blocksResult$call <- NULL

  ##final model:
  ##============

  plmResult <- NULL
  plmClass <- NULL

  ##if empty:
  if( length(bestXs)==0 ){ cat("\n   Final model is empty \n\n") }

  ##if non-empty (estimate final model):
  if( length(bestXs)>0 ){

    ##create formula:
    
    # replace following code:
    # myformula <- paste0(yName, " ~ ", bestXs)
    # if( length(bestXs)>1 ){
    #   for(i in 2:length(bestXs)){
    #     myformula <- paste0(myformula, " + ", bestXs[i])
    #   }
    # }
    # myformula <- as.formula(myformula)
    myformula <- paste(bestXs, collapse = " + ")
    myformula <- paste(yName, "~", myformula, sep = " ")
    myformula <- as.formula(myformula)

    ##estimate final model:
    listOfArgs$formula <- myformula
    plmResult <- do.call("plm", listOfArgs)

    ##extract class:
    plmClass <- class(plmResult)
  }

  ##return result:
  ##==============

  result <- c(blocksResult, plmResult)
  if(!is.null(plmClass)){ class(result) <- plmClass }
##ad hoc visualisation of results:
summary(result)
#  return(result)
#
#} #close
