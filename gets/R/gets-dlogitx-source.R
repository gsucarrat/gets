###########################################################
## This file contains the source of the dlogitx()
## functions.
##
## Created 8 December 2020, Oslo.
##
## Contents:
##
## dlogitxSim()
## logit()
## dlogitx()
## coef.dlogitx         #extraction functions
## fitted.dlogitx       #(all are S3 methods)
## gets.dlogitx
## logLik.dlogitx
## plot.dlogitx
## WIP: predict.dlogitx
## print.dlogitx
## summary.dlogitx
## toLatex.dlogitx
## vcov.dlogitx
##
##
###########################################################

###########################################################
## The function dlogitxSim() simulates from an
## autoregressive logit model with covariates.
##
## Function arguments:
##
## n:           integer, the number of observations to
##              simulate from
## intercept:   the value of the intercept in the logit
##              specification
## ar:          numeric vector w/value(s) of the AR-parameters
## xreg:        numeric vector, e.g. theta1*x1+theta2*x2
## verbose:     logical, whether to return vI only or all
##              the info from the simulations
## as.zoo:      logical, whether the returned object should
##              be of class 'zoo' or not
##
## The conditional probability of yt = 1:
##
##    pi1t = 1/( 1 + exp(-ht) )
##
## The logit specification:
##
##    ht = intercept + ar-terms + covariates
##       = intercept + theta1*y_1 + ... + thetap*y_p + covariates
##
###########################################################

dlogitxSim <- function(n, intercept=0, ar=NULL, xreg=NULL,
  verbose=FALSE, as.zoo=TRUE)
{

  ##orders:
  p.order <- ifelse(is.null(ar), 0, length(ar))
  #q.order <- length(ma) (for the future)
  maxpq <- max(1, p.order)
  nplusMaxpq <- n + maxpq

  ##xreg:
  if( is.null(xreg) ){
    xreg <- rep(0, n)
    xregMean <- 0
  }else{
    xreg <- as.vector(xreg)
    xregMean <- mean(xreg)
  }

  ##p.order = 0 (recursion not necessary):
  if( p.order == 0 ){
    ht <- intercept + xreg
    pi1 <- 1/(1 + exp(-ht))
    vItmp <- runif(n)-1
    vI <- as.numeric( vItmp + pi1 > 0)
  }
  
  ##if p.order > 0 (recursion necessary):
  if( p.order > 0 ){

    ##ht:
    htmean <- (intercept + xregMean)/(1 - sum(ar))
    ht <- rep(htmean, n)

    #pi1, vI:
    vI <- rep(0, n)
    pi1 <- 1/(1 + exp(-ht))
    vItmp <- runif(n)-1

    ##recursion:
    xregsum <- intercept + xreg
    arsum <- 0
    for(i in c(maxpq+1):n){
      arsum <- sum( ar*vI[ c(i-1):c(i-p.order) ] )
      ht[i] <- xregsum[i] + arsum
      pi1[i] <- 1/(1+exp(-ht[i]))
      vI[i] <- as.numeric( vItmp[i] + pi1[i] > 0)
    }

  } #end if p.order > 0
  
  ## result:
  if( verbose ){
    result <- cbind(vI, pi1, ht)
    colnames(result) <- c("I", "pi1", "h")
  }else{
    result <- vI
  }
  if( as.zoo ){
    result <- as.zoo(result)
  }
  return(result)

} #end dlogitxSim

###########################################################
## The function logit() estimates a logit model
##
## Function arguments:
##
## y:               vector with the binary data
## x:               matrix of covariates
## initial.values:  initial values of the parameters
## lower:           lower bounds of the permissible parameter values
## upper:           upper bounds of the permissible parameter values
## method:          estimation method, 1=only estimates,
##                  2=ordinary vcov, 3=hac vcov
## lag.length:      integer controlling the bandwidth when method=3
## control:         options passed on to nlminb()
## eps.tol:         pi0 is computed as pmax(pi0, eps.tol) to
##                  ensure log(pi0) is not too small
## solve.tol:       tolerance for singularity when inverting
##                  a matrix with solve()
##
## The conditional probability of yt = 1:
##
##    pi1t = 1/( 1 + exp(-ht) )
##
## The logit specification:
##
##    ht = intercept + ar-terms + covariates
##       = intercept + theta1*y_1 + ... + thetap*y_p + covariates
##
############################################################

logit <- function(y, x, initial.values=NULL, lower=-Inf, upper=Inf,
  method=2, lag.length=NULL, control=list(), eps.tol=.Machine$double.eps,
  solve.tol=.Machine$double.eps)
{
  ##initiate:
  ##=========
  
  out <- list()
  out$n <- length(y)
  if( is.null(x) ){
    out$k <- 0
  }else {
    out$k <- NCOL(x)
  }
  out$df <- out$n - out$k

  ##initial values: 
  if( is.null(initial.values) && out$k>0 ){
    out$initial.values <- rep(0.1, out$k)
  }else{
    out$initial.values <- initial.values
  }

  ##bounds, tolerance:
  out$lower <- lower
  out$upper <- upper
  out$eps.tol <- eps.tol
  
  ##create log-likelihood function:
  ##===============================
    
  logitObjective <- function(pars, aux, minimise=TRUE){

    ##check parameters:
    parsOK <- TRUE
    if( any(is.na(pars)) || any(pars<aux$lower) || any(pars>aux$upper) ){
      parsOK <- FALSE
    }

    ##compute log-likelihood:
    if( parsOK ){
      ht <- as.vector( aux$x %*% pars )
      pi1that <- 1/(1+exp(-ht))
      pi0that <- pmax(1-pi1that, aux$eps.tol)
      logl <- sum( aux$y*log(pi1that) + aux$oneminy*log(pi0that) )
      if(is.nan(logl) || is.na(logl) || abs(logl) == Inf){
        logl <- aux$logl.penalty
      }
    }else{
      logl <- aux$logl.penalty
    }
    if(minimise){ logl <- -logl }

    ##result:
    return(logl)

  }

  ##create aux:
  ##===========

  aux <- out
  aux$y <- y
  aux$x <- x
  aux$oneminy <- 1-y
  aux$penalty.value <- out$n*log(0.5)
  if( out$k>0 ){
    aux$penalty.value <- logitObjective(out$initial.values, aux)
  }
  
  ##estimate:
  ##=========

  ##w/regressors:
  if( out$k>0 ){
    result <- nlminb(out$initial.values, logitObjective, aux=aux, lower=lower,
      upper=upper, control=control)
    names(result)[1] <- "coefficients"
    result$objective <- -result$objective #due to minimisation
    names(result)[2] <- "logl"
    ht <- as.vector(x %*% result$coefficients)
    result$fit <- 1/(1+exp(-ht)) #pi1hat
  }
  
  ##if there are no regressors:
  if( out$k==0 ){
    result <- list()
    result$coefficients <- NULL
    result$logl <- out$n*log(0.5)
    result$convergence <- 0
    result$iterations <- 0
    result$evaluations <- c(1,0)
    names(result$evaluations) <- c("function", "gradient")
    #result$message <- "convergence (no iterations)"
    result$fit <- rep(0.5, out$n)  
  }

  ##merge out with result:
  result <- c(out, result)

  ##vcov:
  ##=====

  ##compute the hessian:
  if( method>1 && result$k>0 ){
    xadj <- (result$fit^2 - result$fit)*x
    hessian <- crossprod(xadj, x)
    hessianinv <- solve(hessian, tol=solve.tol)
  }  

  ##ordinary vcov:
  if( method==2 && result$k>0 ){
    result$vcov <- -hessianinv
  }  

  ##hac coefficient covariance:
  if( method==3 && result$k>0 ){
  
    ##lag length (bandwidth):
    if( is.null(lag.length) ){
      ##EViews, see Wooldridge (2009, p. 430):
      iL <- as.integer(4*(result$n/100)^(2/9)) 
    }else{
      iL <- as.integer(lag.length)
    }
    result$lag.length <- iL
    vW <- 1 - 1:iL/(iL+1) #weights
  
    ##compute the "meat":
    mS <- cbind((y - result$fit) * x)
    mSigmahat <- crossprod(mS)
    for(l in 1:iL){
      mS0 <- cbind(mS[-c(1:l),])
      mSj <- cbind(mS[-c(c(result$n-l+1):result$n),])
      mSigmahat <-
        mSigmahat + vW[l] * (crossprod(mS0,mSj) + crossprod(mSj,mS0))
    }

    ##compute vcov:
    result$vcov <- hessianinv %*% mSigmahat %*% hessianinv
  }  
  
  ##return result:
  ##==============

  return(result)

} #close logit function


############################################################
## estimate model of class "dlogitx":

dlogitx <- function(y, intercept=TRUE, ar=NULL, ewma=NULL,
  xreg=NULL, vcov.type=c("ordinary", "hac"), lag.length=NULL,
  initial.values=NULL, lower=-Inf, upper=Inf, control=list(),
  eps.tol=.Machine$double.eps, solve.tol=.Machine$double.eps,
  plot=NULL)
{
  ##auxiliary list:
  aux <- list()
  aux$call <- sys.call()  
  aux$date <- date()
  aux$control <- control
  
  ##regressand, regressors:
  tmp <- regressorsMean(y, mc=intercept, ar=ar, ewma=ewma, mxreg=xreg,
    return.regressand=TRUE, return.as.zoo=TRUE,
    na.trim=TRUE,
    na.omit=FALSE)
  whereMconst <- which( colnames(tmp)=="mconst" )
  if( length(whereMconst)>0 ){
    colnames(tmp)[ whereMconst ] <- "intercept"
  }
  aux$y <- coredata(tmp[,1])
  aux$y.name <- colnames(tmp)[1]
  aux$y.index <- index(tmp)
  if( NCOL(tmp)>1 ){
    aux$mX <- cbind(coredata(tmp[,-1]))
    aux$mXnames <- colnames(tmp)[-1]
    colnames(aux$mX) <- NULL
    aux$mXncol <- NCOL(aux$mX)
  }

  ##determine estimation method/vcov type:
  types <- c("ordinary", "hac")
  whichType <- charmatch(vcov.type[1], types)
  aux$logit.method <- whichType + 1 

  ##estimate:
  result <- logit(aux$y, aux$mX, initial.values=initial.values,
    lower=lower, upper=upper, method=aux$logit.method,
    lag.length=lag.length, control=control, eps.tol= eps.tol,
    solve.tol=solve.tol)
  names(result$coefficients) <- aux$mXnames
  colnames(result$vcov) <- aux$mXnames
  rownames(result$vcov) <- aux$mXnames
  result$fit <- zoo(result$fit, order.by=aux$y.index)
  result <- c(aux, result)
  
  ##plot:
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.dlogitx(result) }

  ##return result:
  
  class(result) <- "dlogitx"
  return(result)
  
} #close dlogitx()
    
############################################################
## extract coefficients
coef.dlogitx <- function(object, ...)
{
  object$coefficients
} #close coef.dlogitx

############################################################
## extract fitted probabilities
fitted.dlogitx <- function(object, zero.prob=FALSE, ...)
{
  if( zero.prob ){
    result <- 1-object$fit
  }else{
    result <- object$fit
  }
  return(result)
} #close fitted.dlogitx

############################################################
## do gets on a logitx object
gets.dlogitx <- function(x, t.pval=0.05, wald.pval=t.pval,
  do.pet=TRUE, keep=NULL, include.gum=FALSE, include.1cut=TRUE,
  include.empty=FALSE, max.paths=NULL, turbo=TRUE,
  print.searchinfo=TRUE, plot=NULL, alarm=FALSE, ...)
{
  ##logit() arguments:
  ##------------------
  
  callArgs <- names(x$call)
  userEstArgs <- list()
  userEstArgs$name <- "logit"
  if( "lower" %in% callArgs ){ userEstArgs$lower <- x$lower }
  if( "upper" %in% callArgs ){ userEstArgs$lower <- x$upper }
  userEstArgs$method <- x$logit.method
  if( "lag.length" %in% callArgs ){ userEstArgs$lag.length <- x$lag.length }
  userEstArgs$control <- x$control
  if( "eps.tol" %in% callArgs ){ userEstArgs$eps.tol <- eval(x$call$eps.tol) }
  if( "solve.tol" %in% callArgs ){
    userEstArgs$solve.tol <- eval(x$call$solve.tol)
  }

  ##do gets:
  ##--------

  vY <- x$y
  mX <- x$mX
  result1 <- getsFun(vY, mX, user.estimator=userEstArgs, 
    gum.result=x, t.pval=t.pval, wald.pval=wald.pval, do.pet=do.pet, 
    keep=keep, include.gum=include.gum, include.1cut=include.1cut,
    include.empty=include.empty, max.paths=max.paths, turbo=turbo,
    print.searchinfo=print.searchinfo, alarm=alarm)
  result1$call <- NULL

  ##create gum result:
  ##------------------

  if( x$k>0 ){
    gum.result <- matrix(0, x$k, 6)
    colnames(gum.result) <- c("reg.no.", "keep", "coef", "std.error",
      "t-stat", "p-value")
    rownames(gum.result) <- x$mXnames
    gum.result[,"reg.no."] <- 1:x$k
    if( !is.null(keep) ){ gum.result[keep,"keep"] <- 1 }
    gum.result[,"coef"] <- x$coefficients
    gum.result[,"std.error"] <- sqrt(diag(x$vcov))
    gum.result[,"t-stat"] <- gum.result[,"coef"]/gum.result[,"std.error"]
    gum.result[,"p-value"] <-
      pt(abs(gum.result[,"t-stat"]), df=x$df, lower.tail=FALSE)    
    result1$gum.result <- gum.result
  }
    
  ##estimate specific:
  ##------------------

  vY <- cbind(zoo(vY, order.by=x$y.index))
  colnames(vY) <- x$y.name
  xfinal <- result1$terminals[[ result1$best.terminal ]]
  vcov.type <- c("none", "ordinary", "hac")[ x$logit.method ]
  solve.tol <- ifelse( is.null(userEstArgs$solve.tol),
    .Machine$double.eps, userEstArgs$solve.tol)

  ##empty final model:
  if( length(xfinal)==0 ){
    result2 <- dlogitx(vY, intercept=FALSE, vcov.type=vcov.type,
      lag.length=x$lag.length, lower=x$lower, upper=x$upper,
      control=x$control, eps.tol=x$eps.tol, solve.tol=solve.tol, 
      plot=plot)
  }
  
  ##non-empty final model:
  if( length(xfinal)>0 ){
    mX <- cbind(mX[,xfinal])
    colnames(mX) <- x$mXnames[ xfinal ]
    mX <- zoo(mX, order.by=x$y.index)
    result2 <- dlogitx(vY, intercept=FALSE, xreg=mX, vcov.type=vcov.type,
      lag.length=x$lag.length, lower=x$lower, upper=x$upper,
      control=x$control, eps.tol=x$eps.tol, solve.tol=solve.tol, 
      plot=plot)
  }
  
  ##return result:
  ##--------------
  
  result2$call <- NULL
  result <- c(result1, result2)
  class(result) <- "dlogitx"
  return(result)  
 
} #close gets.dlogitx() 

############################################################
## extract log-likelihood
logLik.dlogitx <- function(object, ...)
{
  result <- object$logl
  attr(result, "df") <- length(object$coefficients)
  attr(result, "nobs") <- object$n
  class(result) <- "logLik"
  return(result)
}

############################################################
## plot fitted probabilities
plot.dlogitx <- function(x, ...)
{
  plot(x$fit, ylab="probability", xlab="", col="blue")
} #close plot.dlogitx

############################################################
## print estimation result of model of class "dlogitx":
print.dlogitx <- function(x, signif.stars=TRUE, ...)
{
  ##header:
  ##-------
  
  ##first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Dependent var.:", x$y.name, "\n")
  cat("Method: Maximum Likelihood (logit) \n")
  cat("Variance-Covariance:",
    switch(x$logit.method, "1" = "Ordinary", "2" = "Ordinary",
    "3" = "HAC, Newey and West (1987)"), "\n")
  cat("No. of observations:", x$n, "\n")

  ##sample info:
  isRegular <- is.regular(x$fit, strict=TRUE)
  isCyclical <- frequency(x$fit) > 1
  indexTrimmed <- index(x$fit)
  if(isRegular && isCyclical){
    cycleTrimmed <- cycle(x$fit)
    startYear <- floor(as.numeric(indexTrimmed[1]))
    startAsChar <- paste(startYear,
      "(", cycleTrimmed[1], ")", sep="")
    endYear <- floor(as.numeric(indexTrimmed[length(indexTrimmed)]))
    endAsChar <- paste(endYear,
      "(", cycleTrimmed[length(indexTrimmed)], ")", sep="")
  }else{
    startAsChar <- as.character(indexTrimmed[1])
    endAsChar <- as.character(indexTrimmed[length(indexTrimmed)])
  }
  cat("Sample:", startAsChar, "to", endAsChar, "\n")

  ##if gets modelling:
  ##------------------

  if( !is.null(x$gum.result) ){

    cat("\n")
    cat("Start model (GUM):\n")
    cat("\n")
    printCoefmat(x$gum.result, tst.ind=c(1,2), signif.stars=signif.stars)

    ##paths:
    cat("\n")
    cat("Paths searched: \n")
    cat("\n")
    if( is.null(x$paths) || length(x$paths)==0 ){
      print(NULL)
    }else{
      for(i in 1:length(x$paths)){
        cat("path",i,":",x$paths[[i]],"\n")
      }
    } #end if(is.null(x$paths))
  
    ##terminal models and results:
    if( !is.null(x$terminals) && length(x$terminals)>0 ){
      cat("\n")
      cat("Terminal models: \n")
      if(!is.null(x$terminals)){
        cat("\n")
        for(i in 1:length(x$terminals)){
          cat("spec",i,":",x$terminals[[i]],"\n")
        }
      }
    }
    if( !is.null(x$terminals.results) ){
      cat("\n")
      printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
        signif.stars=FALSE)
    }

  } #close if gets modelling

  ##results matrix:
  ##---------------
  
  if( length(x$coefficients)>0 ){  

    resultsmat <- matrix(NA, length(x$coefficients), 4)
    colnames(resultsmat) <- c("coef", "std.error", "t-stat", "p-value")
    rownames(resultsmat) <- names(x$coefficients)
    resultsmat[,"coef"] <- x$coefficients
    resultsmat[,"std.error"] <- sqrt(diag(x$vcov))
    resultsmat[,"t-stat"] <- resultsmat[,1]/resultsmat[,2]
    resultsmat[,"p-value"] <-
      pt(abs(resultsmat[,3]), df=x$df, lower.tail=FALSE) 
    
    cat("\n")
    if( is.null(x$gum.result) ){
      cat("Estimation results:\n")
    }else{
      cat("Final model:\n")
    }
    cat("\n")
    printCoefmat(resultsmat, signif.stars=signif.stars)
  
  }else{
    cat("\n")
    cat("   The empty model\n")
  }

  ##goodness-of-fit:
  ##----------------
  
  gof <- matrix(NA, 1, 1)
  rownames(gof) <- paste0("Log-lik.(n=", x$n, ")")
  colnames(gof) <- ""
  gof[1,1] <- x$logl
  printCoefmat(gof, digits=6, signif.stars=signif.stars)

} #end print.dlogitx


############################################################
## summarise output
summary.dlogitx <- function(object, ...)
{
  summary.default(object)
} #end summary.arx

############################################################
### LaTeX code (equation form)
toLatex.dlogitx <- function(object, digits=4, gof=TRUE,
  nonumber=FALSE, nobs="T", ...)
{

  ##probability equation:
  ##---------------------
  
  txtAddNonumber <- ifelse(nonumber, " \\nonumber ", "")
  probtxt <- paste0("  Pr(", object$y.name,
    "_t = 1| ...) &=& \\frac{1}{1 + \\exp(-\\widehat{h}_t)}",
    txtAddNonumber, " \\\\[2mm]", " \n")

  
  ##logit equation:
  ##---------------

  ##y name:
  hName <- paste0("\\widehat{h}_t")

  ##coefs, coef names, std.errors:
  coefs <- coef.dlogitx(object)
  coefsNames <- names(coefs)
  whereIntercept <- which( coefsNames=="intercept" )
  if( whereIntercept > 0 ){ coefsNames[ whereIntercept ] <- "" }
  coefs <- as.numeric(coefs)
  stderrs <- as.numeric(sqrt(diag(vcov.dlogitx(object))))

  ##equation (main part):
  eqtxt <- NULL
  if(length(coefs) > 0){
    for(i in 1:length(coefs) ){
      ifpluss <- ifelse(i==1, "", " + ")
      eqtxt <- paste(eqtxt,
        ifelse(coefs[i] < 0, " - ", ifpluss), "\\underset{(",
        format(round(stderrs[i], digits=digits), nsmall=digits), ")}{",
        format(round(abs(coefs[i]), digits=digits), nsmall=digits), "}",
        coefsNames[i], sep="")
    }
  }

  ##equation (put parts together):
  txtAddNonumber <- ifelse(nonumber, " \\nonumber ", "")
  eqtxt <- paste0("  ", hName, " &=& ", eqtxt, "", txtAddNonumber,
    " \\\\[2mm]", " \n")

  ##goodness of fit:
  ##----------------

  goftxt <- paste("   && LogL=",
    format(round(as.numeric(logLik.dlogitx(object)), digits=digits), nsmall=digits),
    " \\qquad ", nobs, " = ", object$n, " \\nonumber \n", sep="")

  ##print code:
  ##-----------

  cat("\\begin{eqnarray}\n")
  cat(probtxt)
  cat(eqtxt)
  cat(goftxt)
  cat("\\end{eqnarray}\n")

} #end toLatex.dlogitx

############################################################
## variance-covariance extraction function
vcov.dlogitx <- function(object, ...)
{
  object$vcov
} #end vcov.dlogitx
  