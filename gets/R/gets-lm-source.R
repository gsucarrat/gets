####################################################
## This file contains the lm-source of the gets
## package.
##
## FUNCTIONS:
##
## as.arx.lm()
## gets.lm()
## isat.lm()
##
####################################################

##==================================================
## S3 method: convert 'lm' object to 'arx' object:
as.arx.lm <- function(object, ...)
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
    attr(x, "assign") <- NULL #remove 'assign' attribute
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
gets.lm <- function(x, keep=NULL, include.1cut=TRUE,
  print.searchinfo=TRUE, ...)
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

  ##re-name x, create start model object:
  object <- x
  xGUM <- x
  
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
    
  ##------------------------------
  ## 2 print start model (the gum)
  ##------------------------------
  
  gum <- NULL #make sure the 'gum' exists
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
    cat("\n")
    cat("Start model (GUM):\n")
    cat("\n")
    printCoefmat(gum, cs.ind=c(3,4), tst.ind=c(5), signif.stars=TRUE,
      P.values=TRUE)
    cat("\n")

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
##for testing purposes:
#    print.searchinfo=print.searchinfo)
  getsResult$call <- NULL
  getsResult$start.model <- gum

  ##print paths, terminals and retained regressors:
  if( print.searchinfo && !is.null(getsResult$terminals.results) ){

    ##paths:
    if( length(getsResult$paths)>0 ){
      cat("\n")
      for(i in 1:length(getsResult$paths)){
        txt <- paste0(getsResult$paths[[i]], collapse=" ")
        txt <- paste0("  Path ", i, ": ", txt)    
        cat(txt, "\n")
      }
    }

    ##terminals:
    cat("\n")
    cat("Terminal models:\n")
    cat("\n")
    print(getsResult$terminals.results)    
      
    ##retained regressors:
    cat("\n")
    cat("Retained regressors (final model):\n")
    cat("\n")
    if( length(getsResult$specific.spec)==0 ){
      cat("  none\n")
    }else{
      cat(paste0("  ", xNames[as.numeric(getsResult$specific.spec)]), "\n")
    }

  } #end if( print.searchinfo )

  ##messages:
  if( print.searchinfo && !is.null(getsResult$messages) ){
    cat("\n")
    cat("Messages:\n")
    cat("\n")
    cat(getsResult$messages)
    cat("\n")
  }
  
  ##-----------------------------
  ## 5 estimate final model
  ##-----------------------------

  ##search not undertaken:
  ##----------------------

  if( is.null(getsResult$terminals.result) ){
    result <- xGUM
    result <- c(getsResult, xGUM)  
  }

  ##search undertaken:
  ##------------------
  
  if( !is.null(getsResult$terminals.result) ){
      
    ##create re-named y-variable:
    assign(yName, y)
    
    ##if empty final model:
    if( length(getsResult$specific.spec)==0 ){
      lmformula <- paste0(yName, " ~ NULL - 1")
    }
    
    ##if non-empty final model:
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

  } #end if( !is.null(getsResult$terminals.result) )    
  
  ##-----------------------------
  ## 6 return result
  ##-----------------------------
  
  class(result) <- "lm"
  return(result)

} #close gets.lm()

##==================================================
## isat modelling of 'lm' objects:
isat.lm <- function(y, ar=NULL, ewma=NULL, 
                    iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=NULL,
                    ratio.threshold=0.8, max.block.size=30, t.pval=0.001,
                    wald.pval=t.pval, vcov.type=c("ordinary", "white", "newey-west"),
                    do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL,
                    normality.JarqueB=NULL, info.method=c("sc", "aic", "hq"), 
                    user.diagnostics=NULL, user.estimator=NULL, gof.function=NULL, 
                    gof.method=c("min","max"), include.gum=NULL,
                    include.1cut=FALSE, include.empty=FALSE, max.paths=NULL,
                    parallel.options=NULL, turbo=FALSE, tol=1e-07, LAPACK=FALSE,
                    max.regs=NULL, print.searchinfo=TRUE, plot=NULL, alarm=FALSE, ...){

  # Checks
  if(!is.null(y$weights)){stop("Usage of weights is not yet implemented in isat. Please estimate the lm object without weights.")}

  data <- stats::model.frame(y)
  dep_var <- data.frame(y = data[,1])
  names(dep_var) <- names(data)[1]

  mxreg <- stats::model.matrix(y$terms, data)

  # Deal with intercept
  mc <- ifelse(attr(y$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept

  out <- isat(y = dep_var, mxreg = mxreg, mc = mc,
              ar, ewma, iis, sis, tis, uis, blocks,
              ratio.threshold, max.block.size, t.pval,
              wald.pval, vcov.type,
              do.pet, ar.LjungB, arch.LjungB,
              normality.JarqueB, info.method,
              user.diagnostics, user.estimator, gof.function,
              gof.method, include.gum,
              include.1cut, include.empty, max.paths,
              parallel.options, turbo, tol, LAPACK,
              max.regs, print.searchinfo, plot, alarm)
  return(out)
  
} #close isat.lm() function
