####################################################
## This file contains the lm-source of the gets
## package.
##
## CONTENTS:
##
## as.arx
## gets.lm
##
####################################################


##==================================================
## convert 'lm' to 'arx':
as.arx <- function(object, ...)
{ 
  ##'lm' object?:
  if( class(object)!="lm" ){
    objectName <- deparse(substitute(object))
    stop(paste0("'", objectName, "' not of class 'lm'"))
  }
  
  ##make y:
  y <- as.vector(object$model[,1])
  yName <- names(object$model)[1]

  ##make x, intercept?:
  coefs <- coef(object)
  if( length(coefs)==0 ){
    x <- NULL
    mc <- FALSE
  }else{
    x <- model.matrix(object)
    attr(x, "assign") <- NULL #remove attribute
    mc <- ifelse(colnames(x)[1] == "(Intercept)", TRUE, FALSE)
    if(mc){ x <- cbind(x[,-1]) }
    if(NCOL(x)==0){ x <- NULL }
  }
      
  ##estimate arx-model:
  result <- arx(y, mc=mc, mxreg=x, ...)
  return(result)
  
} #close as.arx.lm()  


##==================================================
## gets modelling of 'lm' objects:
gets.lm <- function(x, keep=NULL, print.searchinfo=TRUE, ...)
{
  ## contents:
  ## 1 initiate
  ## 2 print start (gum) model
  ## 3 create interface function
  ## 4 gets modelling
  ## 5 estimate final model
  ## 6 result

  ##-----------------------------
  ## 1 initiate
  ##-----------------------------

  ##re-name x:
  object <- x
  
  ##re-define x:
  x <- model.matrix(object)
  attr(x, "assign") <- NULL #remove attribute
  xNames <- colnames(x)

  ##make copy of y:
  y <- as.vector(object$model[,1])
  yName <- names(object$model)[1]
  
  ##startmodel empty?:
  if( NCOL(x)==0 ){ stop("startmodel is empty") }

  ##call arguments
  callAsList <- as.list(object$call)
  listOfArgs <- callAsList[-1]
  listOfArgs$data <- NULL
  if( length(listOfArgs)==0 ){ listOfArgs <- list() }
    
  ##-----------------------------
  ## 2 print start (gum) model
  ##-----------------------------
  
  if( print.searchinfo ){
  
    ##create gum matrix:
    coefs <- coef(object)
    stderrs <- sqrt(diag(vcov(object)))
    tstats <- coefs/stderrs
    dfs <- length(y) - length(coefs)
    colnamesgum <-
      c("reg.no", "keep", "coef", "std.error", "t-stat", "p-value")
    gum <- matrix(NA, length(coefs), length(colnamesgum))
    colnames(gum) <- colnamesgum
    rownames(gum) <- xNames
    gum[,"reg.no"] <- 1:length(xNames)
    gum[,"keep"] <- 0
    if( !is.null(keep) ){ gum[keep,"keep"] <- 1}
    gum[,"coef"] <- coefs
    gum[,"std.error"] <- stderrs
    gum[,"t-stat"] <- tstats
    gum[,"p-value"] <- 2*pt(abs(tstats), df=dfs, lower.tail=FALSE)
    gum <- as.data.frame(gum)
    
    ##print gum:
    message("")
    message("Start model:")
    message("")
    printCoefmat(gum, cs.ind=c(3,4), tst.ind=c(5), signif.stars=TRUE)
    message("")

  } #end if( print.searchinfo )
     
  ##-----------------------------
  ## 3 create interface function
  ##-----------------------------

  ##modify formula:
  listOfArgs$formula <- as.formula( paste0("y ~ x - 1") )
  
  ##create interface estimator:
  lmFun <- function(y, x, listOfArgs=NULL)
  {
    ##empty x:
    if(is.null(x) || NCOL(x)==0){
  
      result <- list()
    	eps <- y
    	result$logl <- sum(dnorm(eps, sd=sd(eps), log=TRUE))
    	result$n <- NROW(eps)
    	result$k <- 0
    	result$df <- result$n - result$k
  
    }else{ ##non-empty x:
  
      tmp <- do.call("lm", listOfArgs)
  
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

  } #close lmFun()

  ##-----------------------------
  ## 4 gets modelling
  ##-----------------------------

  ##do gets modelling:
  getsResult <- getsFun(y, x,
    user.estimator=list(name="lmFun", listOfArgs=listOfArgs,
    envir=environment()), keep=keep,
    print.searchinfo=print.searchinfo, ...)
  getsResult$call <- NULL
  getsResult$start.model <- gum

  ##print the retained regressors:
  if( print.searchinfo ){
    message("")
    message("Retained regressors:")
    message("")
    if( length(getsResult$specific.spec)==0 ){
      message("  none")
    }else{
      message(paste0("  ", xNames[as.numeric(getsResult$specific.spec)]))
    }
  }
  
  ##-----------------------------
  ## 5 estimate final model
  ##-----------------------------
  
  ##create re-named y-variable:
  assign(yName, y)
  
  ##if empty final model:
  ##---------------------
  if( length(getsResult$specific.spec)==0 ){
    lmformula <- paste0(yName, " ~ NULL - 1")
  }
  
  ##if non-empty final model:
  ##-------------------------
  if( length(getsResult$specific.spec)>0 ){
    
    ##is the intercept retained?
    retainedXs <- getsResult$specific.spec
    whereIntercept <- which( xNames[retainedXs] == "(Intercept)" )
    if( length(whereIntercept)>0 ){
      retainedXs <- retainedXs[ -whereIntercept ]
    }

    ##create remaining X's and lm-formula:
    xNames <- make.names(xNames, unique=TRUE)
    lmformula <- paste0(yName, " ~ ")
    if( length(retainedXs)==0 ){
      lmformula <- paste0(lmformula, " NULL")
    }else{
      for(i in 1:length(retainedXs)){
        if( i>1 ){ lmformula <- paste0(lmformula, " + ") } 
        lmformula <- paste0(lmformula, xNames[retainedXs[i]])
        assign(xNames[retainedXs[i]], x[,retainedXs[i]])
      }
    }
              
    ##no intercept?:
    if( length(whereIntercept)==0 ){
      lmformula <- paste0(lmformula, "- 1")
    }    

  } #end if(non-empty)
      
  ##estimate final model:
  listOfArgs$formula <- as.formula(lmformula)
  result <- do.call("lm", listOfArgs)
  result <- c(getsResult, result)  
  class(result) <- "lm"
    
  ##-----------------------------
  ## 6 return result
  ##-----------------------------
  
  return(result)

} #close gets.lm()
