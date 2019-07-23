####################################################
## This file contains the base-source of the gets
## package.
##
## Current version: 0.19
##
## CONTENTS:
##
## 1 INITIATE
## 2 BASE FUNCTIONS
## 3 ARX FUNCTIONS
## 4 GETS FUNCTIONS
## 5 ADDITIONAL CONVENIENCE FUNCTIONS
##
####################################################
##1 INITIATE
####################################################
##
## print startup message
## load required packages
## create new S3 generics/methods
##
####################################################
##2 BASE FUNCTIONS
####################################################
##
## diagnostics
## dropvar
## eqwma
## infocrit
## info.criterion
## leqwma
## ols
## getsFun
##
####################################################
## 3 ARX FUNCTIONS
####################################################
##
## arx
## coef.arx         #extraction functions
## ES
## fitted.arx
## logLik.arx
## plot.arx
## predict.arx
## print.arx
## recursive
## residuals.arx
## sigma.arx
## rsquared
## summary.arx
## VaR
## vcov.arx
##
####################################################
## 4 GETS FUNCTIONS
####################################################
##
## getsm
## getsv
## coef.gets        #extraction functions
## fitted.gets      #(some are S3 methods)
## logLik.gets
## paths
## plot.gets
## predict.gets
## print.gets
## residuals.gets
## summary.gets
## terminals
## vcov.gets
##
####################################################
## 5 ADDITIONAL CONVENIENCE FUNCTIONS
####################################################
##
## periodicdummies #use regular times series to make periodic dummies
## eviews     #export function
## stata      #export function
## printtex   #print latex-code (equation-form)
##
####################################################


####################################################
## 1 INITIATE
####################################################

##==================================================
## print startup message/code contained in the
## gets-internal.R file:

  ##set start-up message:
  txt <- c("\n",
    paste(sQuote("gets"), "version 0.19\n"),
    "\n",
    paste("General-to-Specific (GETS) and indicator saturation methods, see", sQuote("help(gets)"), "for details"),
    "\n",
    paste("CRAN website: https://CRAN.R-project.org/package=gets"),
    paste("An introduction (PDF): https://www.jstatsoft.org/article/view/v086i03"),
    "\n",
    paste("For automatic plotting, set plot=TRUE in options: options(plot = TRUE)"),
    "\n")
  
  ##print message:
  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }

rm(txt) #remove txt from global environment


##==================================================
##load data :
#hpdata <- read.csv("hpdata.csv")
#so2data <- read.csv("so2.csv")
#infldata <- read.csv("inflation-quarterly.csv")
#sp500data <- read.csv("sp500.csv")

##==================================================
##required packages:
require(parallel)
require(zoo)

##==================================================
##create S3 method 'gets':
gets <- function(x, ...){ UseMethod("gets") }
###test the method:
#set.seed(123); y <- arima.sim(list(ar=0.4, ma=0.1), 100)
#mod01 <- arima(y, order=c(1,0,1))
#gets.Arima <- function(x){ print("It works - cool!") }
#gets(mod01)

##for the future?:
##==================================================
##create S3 method 'isat':
#isat <- function(x, ...){ UseMethod("isat") }
###test the method:
#set.seed(123); y <- arima.sim(list(ar=0.4, ma=0.1), 100)
#mod01 <- arima(y, order=c(1,0,1))
#isat.Arima <- function(x){ print("Cool!") }
#isat(mod01)


####################################################
##2 BASE FUNCTIONS
####################################################

##==================================================
##diagnostics checking of residuals
diagnostics <- function(x, ar.LjungB=c(1,0.025), arch.LjungB=c(1,0.025),
  normality.JarqueB=NULL, verbose=TRUE, user.fun=NULL, ...)
{
  ##initiate:
  if( is.null(x$std.residuals) ){
    zhat <- x$residuals
  }else{
    zhat <- x$std.residuals
  }
  diagnosticsGood <- TRUE

  ##serial correlation:
  if( !is.null(ar.LjungB) ){
    ar.LjungBox <- Box.test(zhat, lag=ar.LjungB[1], type="L")
    if( ar.LjungBox$p.value <= ar.LjungB[2] ){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end serial correlation

  ##arch:
  if( diagnosticsGood && !is.null(arch.LjungB) ){
    zhat2 <- zhat^2
    arch.LjungBox <- Box.test(zhat2, lag=arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end arch

  ##normality:
  if( diagnosticsGood && !is.null(normality.JarqueB)){
    n <- length(zhat)
    avgzhat <- mean(zhat) #do I really need this?
    zhat.avgzhat <- zhat-avgzhat #do I really need this?
    zhat.avgzhat2 <- zhat.avgzhat^2
    K <- n*sum(zhat.avgzhat^4)/(sum(zhat.avgzhat2)^2)
    S <- (sum(zhat.avgzhat^3)/n)/(sum(zhat.avgzhat2)/n)^(3/2)
    JB <- (n/6)*(S^2 + 0.25*((K-3)^2))
    JBpval <- pchisq(JB, df = 2, lower.tail=FALSE)
    if(JBpval <= normality.JarqueB){
        diagnosticsGood <- FALSE
        diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  } #end normality

  ##user.fun:
  if( diagnosticsGood && !is.null(user.fun) ){
    userVals <- do.call(user.fun$name, x, envir=.GlobalEnv)
    userVals <- rbind(userVals)
    if( length(user.fun)>=2 ){
      userFunPval <- as.numeric(userVals[,3])
      if( any(userFunPval <= user.fun[[2]]) ){
        diagnosticsGood <- FALSE
      }
    }
  }

  ##if(verbose): return diagnostics table
  if(verbose){
    result <- NULL
    resultRowNames <- NULL
    if(exists("ar.LjungBox")){
      tmp <- as.numeric(ar.LjungBox[1:3])
      resultRowNames <- c(resultRowNames,
        paste("Ljung-Box AR(", tmp[2], ")", sep=""))
      result <- rbind(result, tmp)
    }
    if(exists("arch.LjungBox")){
      tmp <- as.numeric(arch.LjungBox[1:3])
      resultRowNames <- c(resultRowNames,
        paste("Ljung-Box ARCH(", tmp[2], ")", sep=""))
      result <- rbind(result, tmp)
    }
    if(exists("JBpval")){
      tmp <- c(JB, 2, JBpval)
      resultRowNames <- c(resultRowNames,
        paste("Jarque-Bera", sep=""))
      result <- rbind(result, tmp)
    }
    if(exists("userVals")){
      result <- rbind(result, userVals)
      userValsNames <- rownames(userVals)
      if( identical(userValsNames, "userVals") ){
        userValsNames <- user.fun$name
      }
      if( is.null(userValsNames) ){
        userValsNames <- rep(user.fun$name, NROW(userVals))
      }
      resultRowNames <- c(resultRowNames, userValsNames)
    }
    if(!is.null(result)){
      rownames(result) <- resultRowNames
      colnames(result) <- c("Chi-sq", "df", "p-value")
      if(!is.null(user.fun)){ colnames(result)[1] <- "Statistic" }
    }
  } #end verbose

  ##if(!verbose): return logical
  if(!verbose){ result <- diagnosticsGood }

  ##return result:
  return(result)

} #close diagnostics function

##==================================================
##drop variables that cause exact multicolinearity
dropvar <- function(x, tol=1e-7, LAPACK=FALSE, silent=FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0
{
  ## test and match arguments:
  stopifnot(is.matrix(x))
  silent <- as.logical(silent)[1]
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(x, tol=tol, LAPACK=LAPACK)
  if(qr.X$rank == NCOL(x))
    return(x) ## return x if x has full column rank
  if(!silent){ ## message the no. of dropped columns:
    message("regressor-matrix is column rank deficient, so dropping ",
      NCOL(x) - qr.X$rank, " regressors", appendLF=TRUE)
    message("\n", appendLF=FALSE)
  }
#OLD:
#    message(gettextf("regressor-matrix is column rank deficient, so dropping %d regressors",
#                     NCOL(x) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of x:
  newX <- x[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  ## did we succeed? stop-if-not:
  if(qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}

##==================================================
##generate eqwma regressors
eqwma <- function(x, length=5, lag=1, start=1, p=1,
  log=FALSE, abs=FALSE, as.vector=TRUE)
{
#zoo related:
if(is.zoo(x)){
  xindex <- index(x)
  x <- coredata(as.vector(x))
  zoochk <- TRUE
}else{
  x <- coredata(as.vector(x))
  zoochk <- FALSE
}

#compute eqwma:
xn <- length(x)
x <- na.trim(x, sides="left")
xnadj <- length(x)
if(xn > xnadj){nachk <- TRUE}else{nachk <- FALSE}
if(abs){xabs <- abs(x)}else{xabs <- x}
if(p!=1){xabsp <- xabs^p}else{xabsp <- xabs}

out <- NULL
for(j in 1:length(length)){
  movavg <- rollmean(xabsp, length[j], na.pad=TRUE, align="r")
  if(start <= I(length[j]-1)){
    for(i in start:I(length[j]-1)){
      movavg[i] <- sum(xabsp[1:i])/i
    }
  }
  movavg <- c(rep(NA,lag), as.vector(movavg[1:I(length(movavg)-lag)]))
  out <- cbind(out, movavg)
}

##alternative/future/more efficient approach:
#out <- NULL
#for(j in 1:length(length)){
#  movavg <- rollmean(xabsp, length[j], na.pad=TRUE, align="r")
#  if(length[j] > 1){
#    movavgmean <- mean(movavg, na.rm=TRUE)
#    xabsp.pre <- c( rep(movavgmean, I(length[j]-1)), xabsp[1:I(length[j]-1)] )
#    movavg.pre <- rollmean(xabsp.pre, length[j], na.pad=FALSE, align="r")
#    movavg[1:I(length[j]-1)] <- movavg.pre
#  }
#  movavg <- c(rep(NA,lag), as.vector(movavg[1:I(length(movavg)-lag)]))
#  out <- cbind(out, movavg)
#}

#apply log?:
if(log){
  out <- log(out)
  colnames(out) <- paste("logEqWMA(", length, ")", sep="")
}else{
  colnames(out) <- paste("EqWMA(", length, ")", sep="")
}
if(nachk) out <- rbind( matrix(NA, I(xn-xnadj), dim(out)[2]), out)
if(as.vector){
  if( dim(out)[2]==1 ){ out <- as.vector(out) }
}
if(zoochk){ out <- zoo(out, order.by=xindex) }

return(out)
} #end eqwma

##==================================================
##Compute information criterion:
infocrit <- function(x, method=c("sc", "aic", "aicc", "hq"))
{
  method <- match.arg(method)
  ##Schwarch criterion:
  if(method == "sc") infoVal <- -2*x$logl/x$n + x$k*log(x$n)/x$n
  #Akaike criterion:
  if(method == "aic") infoVal <- -2*x$logl/x$n + 2*x$k/x$n
  #AICC of Hurvich and Tsai (1989):
  if(method == "aicc") infoVal <- -2*x$logl/x$n + 2*x$k/(x$n-x$k-1)
  #Hannan-Quinn criterion:
  if(method == "hq") infoVal <- -2*x$logl/x$n + 2*x$k*log(log(x$n))/x$n
  #out:
  return(infoVal)
} #end infocrit

##==================================================
##Compute information criterion:
info.criterion <- function(logl, n = NULL, k = NULL,
  method = c("sc", "aic", "aicc", "hq"))
{
#check arguments:
if (!is.numeric(logl))
  stop(" Non-numeric argument to logl")
if(is.null(n))
  stop(" Value n is NULL")
if(is.null(k))
  stop(" Value k is NULL")

#initiate:
method <- match.arg(method)
#info.val <- NULL

#Schwarch criterion
if(method == "sc") info.val <- -2*logl/n + k*log(n)/n

#Akaike criterion
if(method == "aic") info.val <- -2*logl/n + 2*k/n

#AICC of Hurvich and Tsai (1989):
if(method == "aicc") info.val <- -2*logl/n + 2*k/(n-k-1)

#Hannan-Quinn criterion
if(method == "hq") info.val <- -2*logl/n + 2*k*log(log(n))/n

#out values:
out <- list()
out$method <- method
out$n <- n
out$k <- k
out$value <- info.val

return(out)

} #end info.criterion

##==================================================
##Estimate log(EqWMA) model
leqwma <- function(x, length=5, lag=1, start=1, p=2, as.vector=FALSE)
{
  eqwma(x, length=length, lag=lag, start=start, p=p,
    log=TRUE, abs=TRUE, as.vector=as.vector)
} #end leqwma

##==================================================
##OLS estimation using the QR decomposition
ols <- function(y, x, untransformed.residuals=NULL, tol=1e-07,
  LAPACK=FALSE, method=3, user.fun=NULL, user.options=NULL)
{

  ##outlook:
  ## - remove the method=0 option (and therefore also the
  ## user.fun argument)?
  ## - rename ols to estFun? Split estFun into two functions,
  ## estFun and vcovFun`?

  ##user-specified:
  if(method==0){
    message("NOTE: method=0 will be deprecated in future versions")
    user.options <- c(list(y = y, x = x), user.options)
    out <- do.call(user.fun, user.options, envir = .GlobalEnv)
  }

  ##fastest (usually only for estimates):
  if(method==1){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK)
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
  }

  ##second fastest (slightly more output):
  if(method==2){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$residuals <- y - out$fit
  }

  ##ordinary vcov:
  if(method==3){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k>0){
      out$vcov <- out$sigma2 * out$xtxinv
    }
    out$logl <- -out$n*log(2*out$sigma2*pi)/2 - out$rss/(2*out$sigma2)
  }

  ##white (1980) vcov:
  if(method==4){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k>0){
      out$omegahat <- crossprod(x, x*out$residuals2)
      out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
    }
    out$logl <- -out$n*log(2*out$sigma2*pi)/2 - out$rss/(2*out$sigma2)
  }

  ##newey-west(1987) vcov:
  if(method==5){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k>0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df

    if(out$k>0){
      iL <- round(out$n^(1/4), digits=0)
      vW <- 1 - 1:iL/(iL+1)
      vWsqrt <- sqrt(vW)
      mXadj <- out$residuals*x
      mS0 <- crossprod(mXadj)

      mSum <- 0
      for(l in 1:iL){
        mXadjw <- mXadj*vWsqrt[l]
        mXadjwNo1 <- mXadjw[-c(1:l),]
        mXadjwNo2 <- mXadjw[-c(c(out$n-l+1):out$n),]
        mSum <- mSum + crossprod(mXadjwNo1, mXadjwNo2) + crossprod(mXadjwNo2, mXadjwNo1)
      }

      out$omegahat <- mS0 + mSum
      out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
    } #end if(out$k>0)
    
    out$logl <- -out$n*log(2*out$sigma2*pi)/2 - out$rss/(2*out$sigma2)
  }

  ##log-variance w/ordinary vcov (note: y = log(e^2)):
  if(method==6){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit #residuals of AR-X representation
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k>0){
      out$vcov <- out$sigma2 * out$xtxinv
    }
    ##log-variance part:
    out$Elnz2 <- -log(mean(exp(out$residuals)))
    out$var.fit <- exp(out$fit - out$Elnz2)
    out$std.residuals <- untransformed.residuals/sqrt(out$var.fit)
    out$logl <- -out$n*log(2*pi)/2 - sum(log(out$var.fit))/2 - sum(untransformed.residuals^2/out$var.fit)/2
  }

  ##result:
  return(out)

} #close ols function

##==================================================
##do gets fast, but unfriendly
getsFun <- function(y, x, untransformed.residuals=NULL,
  user.estimator=list(name="ols", tol=1e-07, LAPACK=FALSE, method=3),
  gum.result=NULL, t.pval=0.05, wald.pval=t.pval, do.pet=TRUE,
  ar.LjungB=NULL, arch.LjungB=NULL, normality.JarqueB=NULL,
  user.diagnostics=NULL, gof.function=list(name="infocrit", method="sc"),
  gof.method=c("min","max"), keep=NULL, include.gum=FALSE,
  include.1cut=FALSE, include.empty=FALSE, max.paths=NULL, turbo=FALSE,
  tol=1e-07, LAPACK=FALSE, max.regs=NULL, print.searchinfo=TRUE,
  alarm=FALSE)
{
  ## TO DO:
  ## - rename turbo argument to estimation.turbo?
  ## - turbo: replace length(regsDeleteList) with regsDeleteList.n?
  ## - turbo: redefine regsFun function (careful!: setequal is delicate)
  ## - add checks on y (e.g. is.vector) and x (e.g. is.matrix)?
  ## - let out$specific.spec be equal to the GUM in the case where
  ##   all regressors are significant in the GUM?

  ### ARGUMENTS: ###########

#  info.method <- match.arg(info.method)
  gof.method <- match.arg(gof.method)

  ##y, x, make auxiliary list:
  if( is.null(x) || NCOL(x)==0 ){ stop("GUM regressor matrix is empty") }
  x <- cbind(x) #ensure x is a matrix
  aux <- list()
  aux$y.n <- NROW(y)
  aux$xNCOL <- NCOL(x)

  ##make user-estimator argument:
  if( is.null(user.estimator$envir) ){ user.estimator$envir <- .GlobalEnv }
  userEstArg <- user.estimator
  userEstArg$name <- NULL
  userEstArg$envir <- NULL
  if( length(userEstArg)==0 ){ userEstArg <- NULL }

  ##make gof.function argument:
  if( is.null(gof.function$envir) ){ gof.function$envir <- .GlobalEnv }
  if( gof.function$name=="infocrit" && is.null(gof.function$method) ){
    gof.function$method <- "sc"
  }
  gofFunArg <- gof.function
  gofFunArg$name <- NULL
  gofFunArg$envir <- NULL
  if( length(gofFunArg)==0 ){ gofFunArg <- NULL }

  ##max.paths argument:
  if( !is.null(max.paths) && max.paths < 1 ){
    stop("'max.paths' cannot be smaller than 1")
  }

  ##do diagnostics?:
  if( !is.null(ar.LjungB) || !is.null(arch.LjungB)
    || !is.null(normality.JarqueB) || !is.null(user.diagnostics) ){
      doDiagnostics <- TRUE
  }else{ doDiagnostics <- FALSE }

  ## max.regs:
  if(is.null(max.regs)){ max.regs <- 10*aux$y.n }

  ### INITIALISE ##########

  ##add to auxiliary list:
  aux$mR <- matrix(0, aux$xNCOL, aux$xNCOL)
  diag(aux$mR) <- 1 #restriction matrix for PETs

  ##make out list, add to out list:
  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  out$no.of.estimations <- 0
  out$messages <- NULL
  out$paths <- list() #the paths
  out$terminals <- list() #terminal specifications
  out$terminals.results <- NULL #matrix w/terminals results
  row.labels <- NULL #row labels for out$terminals.results matrix

  ##deletable, non-deletable regressors, re-organise:
  keep <- as.integer(keep) #do not change to as.numeric(NULL) nor as.vector(NULL), since this may affect setequal/!anyNA...etc. inside the turbo
  keep.n <- length(keep)
  gum <- 1:aux$xNCOL
  delete <- setdiff(gum, keep) #integer(0) if empty
  delete.n <- length(delete)

  ## GUM: ##################

  ##estimate GUM:
  if( is.null(gum.result) ){
    est <- do.call(user.estimator$name,
      c(list(y=y,x=x), userEstArg), envir=user.estimator$envir)
    out$no.of.estimations <- out$no.of.estimations + 1
  }else{ est <- gum.result }

  ##do diagnostics:
  if(doDiagnostics){
    gumDiagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
      arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
      verbose=FALSE, user.fun=user.diagnostics)
  }else{ gumDiagnosticsOK <- TRUE }

  ## if GUM passes diagnostic checks:
  if(gumDiagnosticsOK){

    ##record data for Wald-tests (pets) against gum:
    gum.regs <- gum
    gum.coefs <- est$coefficients
    gum.varcovmat <- est$vcov

    ##compute stderrs, t-stats, p-vals:
    stderrs <- sqrt(diag(est$vcov))
    gum.tstat <- est$coefficients/stderrs
    gum.pval <- pt(abs(gum.tstat), est$df, lower.tail=FALSE)*2

#    ##these two lines are repeated later under 1-cut and in
#    ##the multi-path search:
#    insig.regs <- setdiff( which(gum.pval > t.pval), keep)
#    n.paths <- length(insig.regs)

    ##include gum as terminal?:
    #NEW?: if(include.gum || n.paths==0){
    if(include.gum){

      out$terminals[[1]]  <- gum #add gum to list of terminal specs

      ##specification results
      gofValue <- do.call(gof.function$name,
        c(list(x=est),gofFunArg), envir=user.estimator$envir)
      out$terminals.results <- rbind(out$terminals.results,
        c(gofValue, est$logl, est$n, est$k))
      row.labels <- c(row.labels, "spec 1 (gum):")

    } #end if(include.gum)

  }else{
    out$messages <- paste(out$messages,
      "- MGUM does not pass one or more diagnostic checks", sep="")
  }


  ## 1-CUT MODEL #########

  if( gumDiagnosticsOK && delete.n>0 && include.1cut ){

    ##these two lines are repeated later in the multi-path search:
    ##move these two up under gum (see "NEW..")?
    insig.regs <- setdiff( which(gum.pval > t.pval), keep)
    n.paths <- length(insig.regs)

    ##all non-keep regressors significant:
    if(n.paths==0){
      out$messages <- paste(out$messages,
        "- 1-CUT not included (all non-keep regressors are significant)",
        sep="")
    }

    ##one or more non-keep regressor insignificant:
    if(n.paths > 0){

      ##estimate 1cut:
      mXadj <- cbind(x[,-insig.regs])
      est <- do.call(user.estimator$name,
        c(list(y=y,x=mXadj), userEstArg), envir=user.estimator$envir)
      out$no.of.estimations <- out$no.of.estimations + 1

      ##do diagnostics:
      if(doDiagnostics){
        diagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          verbose=FALSE, user.fun=user.diagnostics)
      }else{ diagnosticsOK <- TRUE }

      ## if 1cut passes diagnostic checks:
      if(diagnosticsOK){

        ## do pet (i.e. wald-test):
        if(do.pet){
          mR <- rbind(aux$mR[insig.regs,])
          mRestq <- mR %*% cbind(gum.coefs)
          wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=tol) %*% mRestq
          petOK <- as.logical(wald.pval < pchisq(wald.stat, n.paths, lower.tail = FALSE))
        }else{
          petOK <- TRUE
        } #end if(do.pet)else..

        ##add 1-cut to terminals?:
        if(petOK){

          #add 1cut to list of terminal specs:
          spec.1cut <- setdiff(gum, insig.regs)
          out$terminals[[ length(out$terminals)+1 ]] <- spec.1cut

          ##specification results
          gofValue <- do.call(gof.function$name,
            c(list(x=est),gofFunArg), envir=user.estimator$envir)
          out$terminals.results <- rbind(out$terminals.results,
            c(gofValue, est$logl, est$n, est$k))
          row.labels <- c(row.labels,
            paste("spec ", length(out$terminals), " (1-cut):", sep=""))

        } #end if(petOK)

      } ##end if(diagnosticsOK)

    } ###end if(n.paths > 0)

  } ####end if(1-cut model)


  ## EMPTY MODEL: ################

  if( gumDiagnosticsOK && delete.n>0 && include.empty ){

    ## DO NOT do pet in order to enable reality check!

    ##check if empty = 1-cut:
    if( include.1cut && exists("spec.1cut") ){
      emptyEqualTo1cut <- identical(keep, spec.1cut)
    }else{ emptyEqualTo1cut <- FALSE }

    ##empty equal to 1cut?:
    if( emptyEqualTo1cut ){

        out$messages <- paste(out$messages,
          "- The empty model is equal to the 1-cut model",
          sep="")

    }else{

      ## estimate model:
      mXadj <- cbind(x[,keep])
      #OLD:
      #if( keep.n==0 ){ mXadj <- NULL }else{ mXadj <- cbind(x[,keep]) }
      est <- do.call(user.estimator$name,
        c(list(y=y,x=mXadj), userEstArg), envir=user.estimator$envir)
      out$no.of.estimations <- out$no.of.estimations + 1

      ##do diagnostics:
      if(doDiagnostics){
        diagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          verbose=FALSE, user.fun=user.diagnostics)
      }else{ diagnosticsOK <- TRUE }

      ## if diagnostics are OK:
      if(diagnosticsOK){

        out$terminals[[ length(out$terminals)+1 ]] <- keep #note: integer(0) if keep=NULL
  #OLD:
  #      out$terminals[[ length(out$terminals)+1 ]] <- if( is.null(keep) ){ 0 }else{ keep }

        ##specification results
        gofValue <- do.call(gof.function$name,
          c(list(x=est),gofFunArg), envir=user.estimator$envir)
        out$terminals.results <- rbind(out$terminals.results,
          c(gofValue, est$logl, est$n, est$k))
        row.labels <- c(row.labels,
          paste("spec ", length(out$terminals), " (empty):", sep=""))

      }else{

          out$messages <- paste(out$messages,
            "- Empty model not included (it does not pass one or more diagnostics)",
            sep="")

      } #end if(empty passes diagnostics==TRUE){..}else{..}

    } ##end if( emptyEqualTo1cut )else(...)

  } ###end if(include empty model==TRUE)


## MULTI-PATH SEARCH: #################

insig.regs <- NULL
pathsTerminals <- list()
if( gumDiagnosticsOK && delete.n>1 ){

  ##number of paths:
  insig.regs <- setdiff( which(gum.pval > t.pval), keep)
  if( !is.null(max.paths) ){
    if(max.paths < length(insig.regs)){
      pvalRanksInv <- rank( 1-gum.pval[insig.regs] )
      insig.regs <- insig.regs[ pvalRanksInv <= max.paths ]
    }
  }
  n.paths <- length(insig.regs)

  ## if paths = 0:
  if(n.paths == 0){
    out$messages <- paste(out$messages,
      "- All regressors significant in GUM mean equation", sep="")
  }

  ## if paths > 0:
  if(n.paths > 0){

    if(print.searchinfo){
      message(n.paths, " path(s) to search")
      message("Searching: ", appendLF=FALSE)
    }

    ##initiate bookkeeping of paths:
    #add if(turbo){...}?
    regsDeleteList <- list()
    regsKeepList <- list()
    regsMat <- NULL

    ## paths:
    for(i in 1:n.paths){

      ## print searchinfo:
      if(print.searchinfo){
        newLine <- ifelse(i==n.paths, TRUE, FALSE)
        message(i, " ", appendLF=newLine)
      }

      ## prepare single-path search:
      path <- insig.regs[i]
      delete.adj <- setdiff(delete, insig.regs[i])
      keep.adj <- keep

      ## single-path search of path i:
      for(j in 1:max.regs){

        ##begin turbo:
        if(turbo && j>1){

          ##bookkeeping of paths:
          regsDeleteList.n <- length(regsDeleteList)
          if( regsDeleteList.n==0 || i==1 ){
#          if( length(regsDeleteList)==0 || i==1 ){

            counter <- regsDeleteList.n + 1
#            counter <- length(regsDeleteList)+1
            regsDeleteList[[ counter ]] <- delete.adj
            regsKeepList[[ counter ]] <- keep.adj
            regsMat <- rbind(regsMat, c(i,length(path)))

          }else{

            ##delete list:
            whichOnesInDelete <- which( sapply(regsDeleteList,
              setequal, delete.adj) )
            #OLD:
            #regsFun <- function(x){ setequal(x,delete.adj) }
            #whichOnesInDelete <- which( sapply(regsDeleteList, regsFun) )
            if( length(whichOnesInDelete)==0 ){
              counter <- regsDeleteList.n + 1
#              counter <- length(regsDeleteList)+1
              regsDeleteList[[ counter ]] <- delete.adj
              regsKeepList[[ counter ]] <- keep.adj
              regsMat <- rbind(regsMat, c(i,length(path)))
              regsDeleteAlreadyDone <- FALSE
            }else{
              regsDeleteAlreadyDone <- TRUE
            }

            ##keep list:
            if( regsDeleteAlreadyDone ){

              ##keep already done?
              whichOnesInKeep <- which( sapply(regsKeepList,
                setequal, keep.adj) )
              #OLD:
              #regsFun <- function(x){ setequal(x, keep.adj) }
              #whichOnesInKeep <- which( sapply(regsKeepList, regsFun) )
              whichOne <- intersect(whichOnesInDelete, whichOnesInKeep)
              #faster version of intersect:
              #y[match(as.vector(x), y, 0L)]
              if( length(whichOne) == 1 ){
                regsKeepAlreadyDone <- TRUE
              }else{
                counter <- regsDeleteList.n + 1
#                counter <- length(regsDeleteList)+1
                regsDeleteList[[ counter ]] <- delete.adj
                regsKeepList[[ counter ]] <- keep.adj
                regsMat <- rbind(regsMat, c(i,length(path)))
                regsKeepAlreadyDone <- FALSE
              }

              ##both delete and keep already done:
              if( regsKeepAlreadyDone ){
                spec.adj <- pathsTerminals[[ regsMat[whichOne,1] ]]
                pathtmp <- out$paths[[ regsMat[whichOne,1] ]]
                pathtmp <- pathtmp[ -c(1:regsMat[whichOne,2]) ]
                path <- c(path, pathtmp)
                break # stop single path search
              }

            } #end regsDeleteAlreadyDone

          } ##end bookkeeping of paths

        } ### end turbo

        ## estimate model:
        #regsAdj <- union(delete.adj, keep.adj)
        mXadj <- cbind(x[, union(delete.adj,keep.adj) ])
        est <- do.call(user.estimator$name,
          c(list(y=y,x=mXadj), userEstArg), envir=user.estimator$envir)
        out$no.of.estimations <- out$no.of.estimations + 1

        ##do diagnostics:
        if(doDiagnostics){
          diagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
            arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
            verbose=FALSE, user.fun=user.diagnostics)
        }else{ diagnosticsOK <- TRUE }

        ## move regressor to keep.adj?:
        if(!diagnosticsOK){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*c(-1))
          next #next j
        }

        ## if empty model passes diagnostic checks:
        if(diagnosticsOK){

          ## stop if no more deletable regressors:
          if( length(delete.adj)==0 ){
            spec.adj <- sort(keep.adj)
            break
          } #end if(length(..)==0)

          #for the future?:
          #if( is.null(est$vcov) ){
          #  est$vcov <- vcovFun(est, method="ordinary")
          #}
          #this will speed up estimation whenever diagnosticsOK
          #turns out to be FALSE. Also, it will provide the user
          #with more flexibility in choosing the covariance matrix

          ##compute stderrs, t-stats, p-vals:
          stderrs <- sqrt(diag(est$vcov))
          t.stat <- est$coefficients/stderrs
          p.val <- pt(abs(t.stat), est$df, lower.tail=FALSE)*2

          ## try deleting a regressor:
          if( any( p.val[1:c(length(delete.adj))] > t.pval) > 0 ){

            reg.no <- which.max( p.val[1:c(length(delete.adj))] )

            ## do pet test (i.e. wald-test):
            if(do.pet){
              deleted <- setdiff(delete, delete.adj[-reg.no])
              deleted <- sort(deleted) #sort() needed for correct restrictions
              n.deleted <- length(deleted)
              mR <- rbind(aux$mR[deleted,])
              mRestq <- mR %*% cbind(gum.coefs)
              wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=tol) %*% mRestq
              petOK <- as.logical(wald.pval < pchisq(wald.stat, n.deleted, lower.tail = FALSE))
            }else{
              petOK <- TRUE
            } #end if(do.pet)else..

            ## delete regressor if(petOK), else move to keep:
            if( petOK ){
              path <- union(path, delete.adj[reg.no])
              delete.adj <- delete.adj[-reg.no]
            }else{
              path <- union(path, delete.adj[reg.no]*c(-1))
              keep.adj <- union(delete.adj[reg.no], keep.adj)
              delete.adj <- delete.adj[-reg.no]
            } #end if( petOK )else{..}

          }else{
            spec.adj <- sort(union(delete.adj, keep.adj))
            break
          } #end if( any p-value > t.pval )else(..)

        } ##end if diagnostics are ok

      } ### end single-path search: for(j in..


      #add path to the paths list:
      counter <- length(out$paths)+1
      out$paths[[ counter ]] <- path
      pathsTerminals[[ counter ]] <- spec.adj

      ##check if spec.adj (terminal) is already in out$terminals:
      if( length(out$terminals)==0 ){
        chk.spec <- FALSE
      }else{
        for(l in 1:length(out$terminals)){
          chk.spec <- setequal(spec.adj, out$terminals[[l]])
          if(chk.spec==TRUE){break} #stop for(l in..)
        }
      } #end check

      ##if spec.adj not in out$terminals (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to out$terminals:
        out$terminals[[ length(out$terminals)+1 ]] <- spec.adj
        gofValue <- do.call(gof.function$name,
          c(list(x=est),gofFunArg), envir=user.estimator$envir)
        out$terminals.results <- rbind(out$terminals.results,
          c(gofValue, est$logl, est$n, est$k))
        row.labels <- c(row.labels, paste("spec ", length(out$terminals), ":", sep=""))

      } #end if(chk.spec==FALSE)

    } ##end multi-path search: for(i in 1:n.paths) loop

  } ###end if paths > 0

} #####end if( gumDiagnosticsOK && delete.n>1 )


  ## FIND THE BEST MODEL: ################

  if( !is.null(out$terminals.results) ){

    ##which is the best model(s):
    if( gof.method=="min" ){
      out$best.terminal <- which.min(out$terminals.results[,1])
    }else{
      out$best.terminal <- which.max(out$terminals.results[,1])
    }

    ##check for several minimums:
    if( length(out$best.terminal)>1 ){
      out$messages <- paste(out$messages,
        "- Several 'best' terminals, the first selected", sep="")
    }
    out$best.terminal <- out$best.terminal[1]
    out$specific.spec <- out$terminals[[ out$best.terminal ]] #the winner

    ##'prettify' out$specific.spec:
    if( length(out$specific.spec)==0 ){
      out$specific.spec <- NULL
    }else{
      out$specific.spec <- sort(out$specific.spec)
      names(out$specific.spec) <- colnames(x)[ out$specific.spec ]
    }

    ##'prettify' out$terminals.results and out$paths:
    if( gof.function$name=="infocrit" ){
      col.labels <- c(paste("info(", gofFunArg$method, ")", sep=""),
        "logl", "n", "k")
    }else{
      col.labels <- c("gof-value", "logl", "n", "k")
    }
    if( NCOL(out$terminals.results) != length(col.labels) ){
      col.labels <- c(col.labels[1], rep(NA,NCOL(out$terminals.results)-1))
    }
    colnames(out$terminals.results) <- col.labels
    rownames(out$terminals.results) <- row.labels
    if( length(out$paths)==0 ){ out$paths <- NULL }

  } #end if( !is.null(out$terminals.results) )


  ## OUTPUT ################################

  out$time.finished <- date()
  if(alarm){ alarm() }
  return(out)

} #close getsFun function


####################################################
##2 ARX FUNCTIONS
####################################################

##==================================================
##Estimate AR-X model with log-ARCH-X errors
arx <- function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
  zero.adj=0.1, vc.adj=TRUE,
  vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, user.estimator=NULL, user.diagnostics=NULL,
  tol=1e-07, LAPACK=FALSE, plot=NULL)
{
  ### ARGUMENTS: ###########

  vcov.type <- match.arg(vcov.type)

  ##regressand:
  y.name <- deparse(substitute(y))
  #if(is.ts(y)){ y <- as.zooreg(y) }
  if(is.zoo(y)){ y <- cbind(y) }else{ y <- as.zoo(cbind(y)) }
  #OLD: y <- as.zoo(cbind(y))
  y <- cbind(y)
  if(NCOL(y) > 1) stop("Dependent variable not 1-dimensional")
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y.n <- NROW(y)
  y.index <- index(y)
  y <- coredata(y)

  ##regressors:
  mX <- NULL
  mXnames <- NULL

  ##mean intercept:
  if(identical(as.numeric(mc),1)){
    mX <- cbind(rep(1,y.n))
    mXnames  <- "mconst"
  }

  ##ar terms:
  if(!is.null(ar) && !identical(as.numeric(ar),0) ){
    tmp <- NULL
    nas <- rep(NA, max(ar))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],y[1:c(y.n-i)]))
    }
    tmpfun <- sapply(ar,tmpfun)
    mX <- cbind(mX, tmp)
    mXnames <- c(mXnames, paste("ar", ar, sep=""))
  }

  ##ewma term:
  if(!is.null(ewma)){
    tmp <- do.call(eqwma, c(list(y),ewma) )
    mXnames <- c(mXnames, colnames(tmp))
    colnames(tmp) <- NULL
    mX <- cbind(mX, tmp)
  }

  ##adjust for NAs:
  tmp <- zoo(cbind(y,mX), order.by=y.index)
  tmp <- na.trim(tmp, sides="both", is.na="any")
  y <- tmp[,1]
  y.n <- NROW(y) #re-define y.n
  y.index <- index(y) #re-define y.index
  t1 <- y.index[1]
  t2 <- y.index[y.n]
  y <- coredata(y)
  if(!is.null(mX)){
    mX <- tmp[,2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL
  }

  ##mxreg:
  if(!is.null(mxreg) && !identical(as.numeric(mxreg),0) ){
    #if(is.ts(mxreg)){ mxreg <- as.zooreg(mxreg) }
    mxreg <- as.zoo(cbind(mxreg))
    mxreg.names <- colnames(mxreg)
    if(is.null(mxreg.names)){
      mxreg.names <- paste("mxreg", 1:NCOL(mxreg), sep="")
    }
    if(any(mxreg.names == "")){
      missing.colnames <- which(mxreg.names == "")
      for(i in 1:length(missing.colnames)){
        mxreg.names[i] <- paste("mxreg", i, sep="")
      }
    }
    #mxreg.names <- make.names(mxreg.names)
    mXnames <- c(mXnames, mxreg.names)
    mxreg <- window(mxreg, start=t1, end=t2)
    mxreg <- cbind(coredata(mxreg))
    mX <- cbind(mX, mxreg)

    ##re-adjust for NAs:
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    y <- tmp[,1]
    y.n <- NROW(y) #re-define y.n
    y.index <- index(y) #re-define y.index
    t1 <- y.index[1]
    t2 <- y.index[y.n]
    y <- coredata(y)
    mX <- tmp[,2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL

  } #end if(!is.null(mxreg))

  ##vxreg:
  if(!is.null(vxreg) && !identical(as.numeric(vxreg),0) ){
    vxreg <- as.zoo(cbind(vxreg))
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste("vxreg", 1:NCOL(vxreg), sep="")
    }
    if(any(vxreg.names == "")){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[i] <- paste("vxreg", i, sep="")
      }
    }
    vxreg <- window(vxreg, start=t1, end=t2)
#OLD:
#    vxreg <- cbind(coredata(vxreg))
    colnames(vxreg) <- NULL
  } #end if(!is.null(vxreg))

  ##determine qstat.options:
  if(is.null(qstat.options)){
    if(is.null(ar)){ar.lag <- 1}else{ar.lag <- max(ar)+1}
    if(is.null(arch)){arch.lag <- 1}else{arch.lag <- max(arch)+1}
    qstat.options <- c(ar.lag, arch.lag)
  }

  ##aux: info for getsm/getsv functions
  aux <- list()
  aux$y <- y
  aux$y.index <- y.index
  aux$y.name <- y.name
  aux$y.n <- y.n
  if(!is.null(mX)){
    colnames(mX) <- NULL
    aux$mX <- mX
    aux$mXnames <- mXnames
    aux$mXncol <- NCOL(mX)
  }
  aux$vc <- vc
  aux$zero.adj <- zero.adj
  aux$vc.adj <- vc.adj
  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$user.estimator <- user.estimator
  aux$user.diagnostics <- user.diagnostics
  aux$tol <- tol
  aux$LAPACK <- LAPACK

  ### INITIALISE ##########

  sysCall <- sys.call()
  #for the future: make the following objects part of the out-list?
  vcov.var <- NULL #make sure this object exists
  variance.results <- NULL #make sure this object exists

  #### for the future regarding user.estimator: check if
  #### user.estimator$spec is NULL, "mean", "variance" or "both"
  #### in order to determine what kind of estimator it is

  ##check if mean and log-garch spec:
  meanSpec <- !is.null(mX)
  noVarianceSpec <- if( vc==FALSE && is.null(arch)
    && is.null(asym) && is.null(log.ewma)
    && is.null(vxreg) ){ TRUE }else{ FALSE }

  #### MEAN ###############

  ##estimate:
  if( is.null(user.estimator) ){

    estMethod <- which(vcov.type==c("none", "none", "ordinary",
      "white", "newey-west"))
    out <- ols(y, mX, tol=tol, LAPACK=LAPACK,
      method=estMethod, user.fun=NULL, user.options=NULL)

    ##delete unneeded entries:
    #out$n <- NULL #this might have to be changed in order to enable gum.result in getsFun
    #out$k <- NULL ##this might have to be changed in order to enable gum.result in getsFun
    #out$df <- NULL: Do not delete!
    out$qr <- NULL
    out$rank <- NULL
    out$qraux <- NULL
    out$pivot <- NULL
    out$xtxinv <- NULL
    out$residuals2 <- NULL
    #out$rss <- NULL

    ##add entries if no variance spec:
    if(noVarianceSpec){
      out$var.fit <- rep(out$sigma2, aux$y.n)
      out$std.residuals <- out$residuals/sqrt(out$sigma2)
      aux$loge2.n <- aux$y.n #change to out$n?
    }

  }else{

    ##user-defined estimator:
    out <- do.call(user.estimator$name, list(y, mX),
      envir=.GlobalEnv)
    #delete?:
    if( is.null(out$vcov) && !is.null(out$vcov.mean) ){
      out$vcov <- out$vcov.mean
    }

  } #end if(..)else(..)

  ##make estimation results (a data frame):
  if(meanSpec){
    stderrs <- sqrt(diag(out$vcov))
    colnames(out$vcov) <- mXnames
    rownames(out$vcov) <- mXnames
    t.stat <- out$coefficients/stderrs
    p.val <- pt(abs(t.stat), out$df, lower.tail=FALSE)*2
    out$mean.results <- as.data.frame(cbind(out$coefficients,
      stderrs, t.stat, p.val))
    colnames(out$mean.results) <- c("coef", "std.error",
      "t-stat", "p-value")
    rownames(out$mean.results) <- mXnames
  } #end if(meanSpec)

  ##rename some stuff:
  if( is.null(user.estimator) ){
    outNames <- names(out)
    whereIs <- which(outNames=="vcov")
    if( length(whereIs) > 0 ){ names(out)[whereIs] <- "vcov.mean" }
    whereIs <- which(outNames=="fit")
    names(out)[whereIs] <- "mean.fit"
  }


  #### VARIANCE #############

  ##if log-arch spec:
  if( noVarianceSpec==FALSE ){

    ##regressand
    zero.where <- which(out$residuals==0)
    eabs <- abs(out$residuals)
    if(length(zero.where) > 0){
      eabs[zero.where] <- quantile(eabs[-zero.where], zero.adj)
    }
    loge2 <- log(eabs^2)

    ##regressor matrix:
    vX <- cbind(rep(1,y.n))
    vXnames <- "vconst"

    ##arch terms:
    if(!is.null(arch) && !identical(as.numeric(arch),0) ){
      tmp <- NULL
      nas <- rep(NA, max(arch))
      tmpfun <- function(i){
        tmp <<- cbind(tmp, c(nas[1:i],loge2[1:c(y.n-i)]))
      }
      tmpfun <- sapply(arch,tmpfun)
      vX <- cbind(vX, tmp)
      vXnames <- c(vXnames, paste("arch", arch, sep=""))
    }

    ##asym terms:
    if(!is.null(asym) && !identical(as.numeric(asym),0) ){
      tmp <- NULL
      nas <- rep(NA, max(asym))
      tmpfun <- function(i){
        tmp <<- cbind(tmp, c(nas[1:i],
          loge2[1:c(y.n-i)]*as.numeric(out$residuals[1:c(y.n-i)]<0)))
      }
      tmpfun <- sapply(asym,tmpfun)
      vX <- cbind(vX, tmp)
      vXnames <- c(vXnames, paste("asym", asym, sep=""))
    }

    ##log.ewma term:
    if(!is.null(log.ewma)){
      if(is.list(log.ewma)){
        log.ewma$lag <- 1
      }else{
        log.ewma <- list(length=log.ewma)
      }
      tmp <- do.call(leqwma, c(list(out$residuals),log.ewma) )
      vXnames <- c(vXnames, colnames(tmp))
      colnames(tmp) <- NULL
      vX <- cbind(vX, tmp)
    }

    ##adjust for NAs:
    tmp <- zoo(cbind(loge2,vX), order.by=y.index)
    tmp <- na.trim(tmp, sides="left", is.na="any")
    loge2 <- tmp[,1]
    loge2.n <- NROW(loge2)
    loge2.index <- index(loge2) #re-define y.index
    loge2 <- coredata(loge2)
    vX <- tmp[,2:NCOL(tmp)]
    vX <- coredata(vX)
    vX <- cbind(vX)
    colnames(vX) <- NULL

    ##vxreg:
    if( !is.null(vxreg) && !identical(as.numeric(vxreg),0) ){
      vxreg <- window(vxreg, start=loge2.index[1],
        end=loge2.index[loge2.n])
      vxreg <- cbind(vxreg)
      vX <- cbind(vX, coredata(vxreg))
      vXnames <- c(vXnames, vxreg.names)
      colnames(vxreg) <- vxreg.names

      ##re-adjust for NAs:
      tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
      tmp <- na.trim(tmp, sides="left", is.na="any")
      loge2 <- tmp[,1]
      loge2.n <- NROW(loge2)
      loge2.index <- index(loge2) #re-define index
      loge2 <- coredata(loge2)
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }

    ##aux: more info for getsm/getsv functions
    aux$loge2 <- loge2
    aux$loge2.n <- loge2.n
    aux$vX <- vX
    aux$vXnames <- vXnames
    aux$vXncol <- NCOL(vX)
    aux$arch <- arch
    aux$asym <- asym
    aux$log.ewma <- log.ewma
    aux$vxreg <- vxreg #note: NROW(vxreg)!=NROW(vX) is possible

    ##estimate, prepare results:
    estVar <- ols(loge2, vX, tol=tol, LAPACK=LAPACK,
      method=3)
    s.e. <- sqrt(as.vector(diag(estVar$vcov)))
    estVar$vcov <- as.matrix(estVar$vcov[-1,-1])
    colnames(estVar$vcov) <- vXnames[-1]
    rownames(estVar$vcov) <- vXnames[-1]
    Elnz2 <- -log(mean(exp(estVar$residuals)))
    t.stat <- estVar$coefficients/s.e.
    p.val <- pt(abs(t.stat), estVar$df, lower.tail=FALSE)*2
    if(vc.adj){
      t.stat[1] <- ((estVar$coefficients[1]-Elnz2)^2)/s.e.[1]^2
      p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
      estVar$coefficients[1] <- estVar$coefficients[1] - Elnz2
    }

    ##add/modify entries in out:
    out$n <- estVar$n
    out$vcov.var <- estVar$vcov
    out$var.fit <- exp(estVar$fit - Elnz2)
    out$ustar.residuals <- estVar$residuals
    out$std.residuals <- out$residuals[c(y.n-loge2.n+1):y.n]/sqrt(out$var.fit)
    out$logl <- -loge2.n*log(2*pi)/2 - sum(log(out$var.fit))/2 - sum(out$std.residuals^2)/2
    out$Elnz2 <- Elnz2
    out$variance.results <- as.data.frame(cbind(estVar$coefficients, s.e., t.stat, p.val))
    colnames(out$variance.results) <- c("coef", "std.error", "t-stat", "p-value")
    rownames(out$variance.results) <- vXnames

  } #end if( noVarianceSpec==FALSE )


  ### OUTPUT: ######################

  ##diagnostics:
  out$diagnostics <- diagnostics(out,
    ar.LjungB=c(qstat.options[1],0), arch.LjungB=c(qstat.options[2],0),
      normality.JarqueB=0, user.fun=user.diagnostics, verbose=TRUE)

  ##add NAs to variance series:
  if( noVarianceSpec==FALSE ){
    NAs2add <- rep(NA,y.n-loge2.n)
    out$var.fit <- c(NAs2add, out$var.fit)
    out$ustar.residuals <- c(NAs2add, out$ustar.residuals)
    out$std.residuals <- c(NAs2add, out$std.residuals)
  }

  ##add zoo-indices:
  out$mean.fit <- zoo(out$mean.fit, order.by=y.index)
  out$residuals <- zoo(out$residuals, order.by=y.index)
  if(!is.null(out$var.fit)){
    out$var.fit <- zoo(out$var.fit, order.by=y.index)
  }
  if(!is.null(out$ustar.residuals)){
    out$ustar.residuals <- zoo(out$ustar.residuals, order.by=y.index)
  }
  if(!is.null(out$std.residuals)){
    out$std.residuals <- zoo(out$std.residuals, order.by=y.index)
  }

  ##result:
  out <- c(list(call=sysCall, date=date(), aux=aux), out)
  class(out) <- "arx"

  ##plot:
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.arx(out) }

  ##return result:
  return(out)
} #close arx function

##==================================================
coef.arx <- function(object, spec=NULL, ...)
{
  ##spec argument:
  if(is.null(spec)){
    spec <- switch(as.character(object$call)[1],
      arx="both", getsm="mean", getsv="variance")
    #spec <- "both"
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  } #end if(..)else(..)

  ##mean results:
  result1 <- NULL
  if(spec=="mean" || spec=="both"){
    if(!is.null(object$mean.results)){
      result1 <- object$mean.results[,1]
      names(result1) <- rownames(object$mean.results)
    }
  } #end if(spec==..)

  ##variance results:
  result2 <- NULL
  if(spec=="variance" || spec=="both"){
    if(!is.null(object$variance.results)){
      result2 <- object$variance.results[,1]
      names(result2) <- rownames(object$variance.results)
      if(!is.null(object$Elnz2)){
        result2 <- c(result2,object$Elnz2)
        names(result2)[length(result2)] <- "Elnz2"
      }
    } #end if(..)else(..)
  } #end if(spec==..)

  result <- c(result1,result2)
  return(result)
} #end coef.arx

##==================================================
ES <- function(object, level=0.99, type=7, ...)
{
  ##check whether class is valid:
  classType <- class(object)
  if( !classType %in% c("arx", "gets", "isat") ){
    stop("object not of class 'arx', 'gets' nor 'isat'")
  }

  ##check the risk-levels:
  riskLevel <- 1-level
  if( any(riskLevel > 1) || any(riskLevel < 0) ){
    stop("risk-level(s) must be in the 0 to 1 interval")
  }

  ##fitted sd, standardised residuals, quantile:
  meanFit <- fitted(object, spec="mean")
  sdFit <- sqrt( fitted(object, spec="variance") )
  residsStd <- residuals(object, std=TRUE)
  qValue <- quantile(residsStd, probs=riskLevel, type=type,
    names=FALSE, na.rm=TRUE)
  colNames <- paste("ES", level, sep="")
  mExpShortF <- matrix(NA, length(sdFit), length(colNames))
  for(i in 1:length(colNames)){
    whereExceeds <- which( residsStd < qValue[i] )
    if( length(whereExceeds) == 0 ){
      stop("no standardised residual smaller than ", qValue[i])
    }else{
      ExpShortF <- mean( residsStd[whereExceeds] )
    }
    mExpShortF[,i] <- sdFit*ExpShortF
  }
  colnames(mExpShortF) <- colNames
  if(NCOL(mExpShortF)==1){ mExpShortF <- as.vector(mExpShortF) }
  mExpShortF <- zoo(mExpShortF, order.by=index(sdFit))
  mExpShortF <- meanFit + mExpShortF

  ##return
  return(-mExpShortF)

}  #close ES

##==================================================
fitted.arx <- function(object, spec=NULL, ...)
{
  ##spec argument:
  if(is.null(spec)){
    if(!is.null(object$mean.results)){
      spec <- "mean"
    }
    if(is.null(object$mean.results)
      && !is.null(object$variance.results) ){
      spec <- "variance"
    }
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  result <- NULL

  ##mean:
  if(spec=="mean"){
    result <- object$mean.fit
  }

  ##variance:
  if(spec=="variance"){
    result <- object$var.fit
  }

  ##both:
  if(spec=="both"){
    if(!is.null(object$mean.results)
      && !is.null(object$variance.results) ){
      result <- cbind(object$mean.fit, object$var.fit)
      colnames(result) <- c("yhat","sigma2hat")
    }
  }

  return(result)
} #end fitted.arx

##==================================================
logLik.arx <- function(object, ...)
{
  ## in the future: add a df.method argument with
  ## optional values "mean-coefficients" (default),
  ## "variance-coefficients" and "both"?
  
  result <- object$logl
  attr(result, "df") <- length(object$coefficients)
  attr(result, "nobs") <- object$n
  class(result) <- "logLik"
  return(result)

#OLD:
#  if( is.null(object$aux$user.estimator) ){
#
#    resids <- residuals.arx(object, std=FALSE)
#    sdhat <- sqrt(fitted.arx(object, spec="variance"))
#    sdhat <- na.trim(sdhat)
#    resids <- window(resids, start=index(sdhat)[1],
#      end=index(sdhat)[length(sdhat)])
#    result <- sum(dnorm(resids, sd=sdhat, log=TRUE))
#    attr(result, "df") <- length(coef.arx(object, spec="mean"))
#    attr(result, "nobs") <- length(sdhat)
#    class(result) <- "logLik"
#
#  }else{
#
#    result <- object$logl
#    if( !is.null(result) ){
#      attr(result, "df") <- length(object$coefficients)
#      if( "residuals" %in% names(object) ){
##OLD:
##      if( !is.null(objects$residuals) ){
#        attr(result, "nobs") <- length(na.trim(object$residuals))
#      }
#      class(result) <- "logLik"
#    }
#
#  }
#
#  return(result)

} #close logLik.arx

##==================================================
##plot results from arx
plot.arx <- function(x, spec=NULL, col=c("red","blue"),
  lty=c("solid","solid"), lwd=c(1,1), ...)
{
  ##check whether to plot:
  doPlot <- TRUE #default
  if( is.null(x$mean.results) && is.null(x$variance.results) ){
    doPlot <- FALSE
    warning("No estimated model, so no plot produced")
  }
  if(doPlot && !is.null(x$aux$user.estimator) ){
    doPlot <- FALSE
    warning("User defined estimation, so no plot produced")
  }

  ##proceed with plotting:
  if(doPlot){

    ##lwd argument:
    if(length(lwd)==1){
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd=rep(lwd,2)
    }else if (length(lwd)>2){
      print("lwd needs two arguments, but more provided. First two used.")
      lwd=lwd[1:2]
    }

    ##lty argument:
    if(length(lty)==1){
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty=rep(lty,2)
    }else if (length(lwd)>2){
      print("lty needs two arguments, but more provided. First two used.")
      lty=lty[1:2]
    }

    ##col argument:
    if(length(col)!=2){

      #####randomcol - returns random combination of colours of length 2
      randomcol <- function()
      {
        r.r <- runif(2)
        while(round(r.r[1],1)==round(r.r[2],1)) {#don't want colours too similar so add check
          r.r <- runif(2)
        }
        g.r <- runif(2)
        while(round(g.r[1],1)==round(g.r[2],1)) {#don't want colours too similar so add check
          g.r <- runif(2)
        }
        b.r <- runif(2)
        while(round(b.r[1],1)==round(b.r[2],1)) {#don't want colours too similar so add check
          b.r <- runif(2)
        }
        col<-rgb(runif(2),runif(2),runif(2))
        return(col)
      } #end randomcol

      ###clashcol function - returns clashing (opposite) combination of colours
      clashcol <- function()
      {
        #using http://forum.processing.org/one/topic/the-opposite-of-a-color.html formula for opposite colour
        r.1 <- runif(1)
        g.1 <- runif(1)
        b.1 <- runif(1)
        b.2 <- min(r.1,min(g.1,b.1)) + max(r.1,max(g.1,b.1))
        col <- rgb(c(r.1,b.2-r.1),c(g.1,b.2-g.1),c(b.1,b.2-b.1))
        return(col)
      } #end clashcol

      ##if random:
      if(col[1]=="random") {
        col <- randomcol()
      }else if(col[1]=="awful.clash") {
        col <- clashcol()
      }else{
        print("Wrong number of colours specified; using random set of colours instead.")
        col<-randomcol()
      }
    } #end col argument

    ##spec argument:
    ##(logically, this part should come before the col, lty and lwd arguments)
    if(is.null(spec)){
      if(!is.null(x$mean.results)){
        spec <- "mean"
      }
      if(is.null(x$mean.results)
         && !is.null(x$variance.results) ){
        spec <- "variance"
      }
      if(!is.null(x$mean.results)
         && !is.null(x$variance.results) ){
        spec <- "both"
      }
    }else{
      spec.type <- c("mean", "variance", "both")
      which.type <- charmatch(spec, spec.type)
      spec <- spec.type[which.type]
    }

    ##plot if spec is not NULL:
    if(!is.null(spec)){

      ##if variance modelled, plot square root of fitted variance and absolute residuals against time
      if(spec=="variance" || spec=="both"){
        vfitted <- sqrt(x$var.fit)
        vactual <- abs(x$residuals)
      }

      ##if mean modelled, plot fitted and actual values against time
      if(spec=="mean" || spec=="both"){
        mfitted <- x$mean.fit
        mactual <- zoo(x$aux$y, order.by=x$aux$y.index)
      }
      actual.name <- x$aux$y.name
      residsStd <- x$std.residuals

      ##do the plotting:
      ##----------------

      ##get current par-values:
      def.par <- par(no.readonly=TRUE)

      ##set new par values for plot
      if(spec=="both") {##if both mean and variance modelled, plot both
        par(mfrow=c(3,1))
      }else {##else just plot the one specified
        par(mfrow=c(2,1))
      }

      #set the plot margins:
      par(mar=c(2,2,0.5,0.5))

      ##plot the mean:
      if(spec=="mean" || spec=="both") {##plotting mean variables

        ##check whether ?? zoo object is regular, then plot:
        if(is.regular(mactual)) {
          plot(mactual, main = "",
             ylim=range(min(mactual,mfitted),max(mactual,mfitted)),
             type="l",ylab="",xlab="",col=col[2])
        } else {##if irregular, plot manually
          plot(as.Date(index(mactual)),coredata(mactual), main = "",
             ylim=range(min(mactual,mfitted),max(mactual,mfitted)),
             type="l",ylab="",xlab="",col=col[2])
        }

        ##check whether ?? zoo object is regular, then plot:
        if(is.regular(mfitted)) {
          lines(mfitted,col=col[1])
        } else {
          lines(as.Date(index(mfitted)),coredata(mfitted),col=col[1])
        }
        legend("topleft",lty=lty,lwd=lwd,ncol=2,col=col[c(2,1)],legend=c(actual.name,"fitted"),bty="n")

      } #close mean plotting

      ##plot the variance:
      if(spec=="variance" || spec=="both") {

        ##add comment?
        if(is.regular(vactual)) {
          plot(vactual, main = "",
             ylim=range(min(vactual,vfitted,na.rm=TRUE),max(vactual,vfitted,na.rm=TRUE)),
             type="l",ylab="",xlab="",col=col[2])
        } else {
          plot(as.Date(index(vactual)),coredata(vactual), main = "",
             ylim=range(min(vactual,vfitted,na.rm=TRUE),max(vactual,vfitted,na.rm=TRUE)),
             type="l",ylab="",xlab="",col=col[2])
        }

        ##add comment?
        if(is.regular(vfitted)) {
          lines(vfitted,col=col[1])
        } else {
          lines(as.Date(index(vfitted)),coredata(vfitted),col=col[1])
        }
        legend("topleft",lty=lty,lwd=lwd,ncol=2,col=col[c(2,1)],
          legend=c("abs(residuals)","fitted sd"),bty="n")

      } #close plotting variance parts

      ##if any standardised residuals:
      if(!is.null(residsStd)){
        if(is.regular(residsStd)) {
          plot(residsStd,type="h",col=col[1])
        } else {
          plot(as.Date(index(residsStd)),coredata(residsStd),type="h",col=col[1])
        }
        abline(0,0)
        legend("topleft",lty=1,col=col[1],legend=c("standardised residuals"),bty="n")
      }

      #return to old par-values:
      par(def.par)

    } #close if(!is.null(spec))

  } #close if(doPlot)

} #close plot.arx

##==================================================
## forecast up to n.ahead
predict.arx <- function(object, spec=NULL, n.ahead=12,
  newmxreg=NULL, newvxreg=NULL, newindex=NULL,
  n.sim=1000, innov=NULL, return=TRUE, plot=NULL,
  plot.options=list(), ...)
{

  ## contents:
  ## 0 initialise
  ## 1 if mean spec
  ## 2 if variance spec
  ## 3 if plot=TRUE
  ## 4 if return=TRUE

  ##n.ahead:
  if(n.ahead < 1){ stop("n.ahead must be 1 or greater") }

  ##spec:
  if(is.null(spec)){
    if(!is.null(object$mean.results)) spec <- "mean"
    if(is.null(object$mean.results)
       && !is.null(object$variance.results)) spec <- "variance"
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  } #end if(..)else(..)
  if(is.null(spec)){ stop("No estimated model") }

  ##newindex:
  yInSample <- zoo(object$aux$y, order.by=object$aux$y.index)

  if(!is.null(newindex)){
    yAsRegular <- FALSE
    if( n.ahead!=length(newindex) ){
      stop("length(newindex) must equal n.ahead")
    }
  }

  if( is.null(newindex) && is.regular(yInSample, strict=TRUE) ){
    endCycle <- cycle(yInSample)
    endCycle <- as.numeric(endCycle[length(endCycle)])
    endYear <- floor(as.numeric(object$aux$y.index[object$aux$y.n]))
    yFreq <- frequency(yInSample)
    yhataux <- rep(NA, n.ahead+1)
    yDeltat <- deltat(yInSample)
    if( yDeltat==1 && yFreq==1 ){
      yhataux <- zoo(yhataux,
        order.by=seq(endYear, endYear+n.ahead, by=1))
      yAsRegular <- FALSE
    }else{
      yhataux <- zooreg(yhataux, start=c(endYear, endCycle),
                frequency=yFreq)
      yAsRegular <- TRUE
    }
    yhataux <- yhataux[-1]
    newindex <- index(yhataux)
  }

  if(is.null(newindex)){ newindex <- 1:n.ahead }

  ##----------------
  ## 1 if mean spec
  ##----------------

  outMean <- NULL
  if(spec=="mean" || spec=="both"){
    coefs <- coef.arx(object, spec="mean")

    ##mc:
    if(!is.null(object$call$mc)){
      mconst <- as.numeric(coefs[1])
      mconstIndx <- 1
    }else{
      mconst <- 0
      mconstIndx <- 0
    }

    ##ar:
    arMax <- 0
    arIndx <- max(mconstIndx)
    arCoefs <- NULL ##J-dog change
    if(!is.null(object$call$ar)){
      arEval <- eval(object$call$ar)
      arIndx <- 1:length(arEval) + max(mconstIndx)
      arMax <- max(arEval)
      arCoefs <- rep(0,arMax)
      arCoefs[arEval] <- as.numeric(coefs[arIndx])
    }

    ##ewma:
    ewmaMax <- 0
    ewmaIndx <- max(arIndx)
    if(!is.null(object$call$ewma)) stop("Sorry, 'ewma' not implemented yet")
    #    if(!is.null(object$call$ewma)){
    #      modify arEval
    #      modify arIndx
    #      modify arMax
    #      modify arCoefs
    #    }

    ##backcast length:
    backcastMax <- max(arMax,ewmaMax)

    ##mxreg:
    mxreghat <- rep(0, n.ahead + backcastMax)

    if(!is.null(object$call$mxreg)){

      ##check newmxreg:
      if(is.null(newmxreg)){ stop("'newmxreg' is NULL") }
      if(NROW(newmxreg)!=n.ahead){ stop("NROW(newmxreg) must equal n.ahead") }
      #      if(NROW(newmxreg)!=n.ahead
      #        && NROW(newmxreg)!=1){ stop("NROW(newmxreg) must equal 1 or n.ahead") }
      #      if(NROW(newmxreg)!=n.ahead
      #        && NROW(newmxreg)==1){
      #          newmxreg <- coredata(rbind(as.zoo(newmxreg)))
      #          tmp <- matrix(NA,n.ahead,NCOL(newmxreg))
      #          tmp[1:NROW(tmp),] <- newmxreg
      #          newmxreg <- tmp
      #      }

      ##newmxreg:
      newmxreg <- coredata(cbind(as.zoo(newmxreg)))
      colnames(newmxreg) <- NULL

      ##mxreghat:
      mxregIndx <- c(max(ewmaIndx)+1):length(coefs) ##J-dog change

      mxreghat <-  newmxreg %*% as.numeric(coefs[mxregIndx])
      mxreghat <- c(rep(0,backcastMax),mxreghat)

    } ##end mxreg

    ##prepare prediction:
    yhat <- rep(NA, n.ahead + backcastMax)
    yhat.n <- length(yhat)
#OLD:
#    meanFit <- coredata(fitted.arx(object, spec="mean"))
    if(backcastMax>0) {
      ##actual y-values:
      yhat[1:backcastMax] <- tail(object$aux$y, n=backcastMax)
#        OR:
#        object$aux$y[c(length(object$aux$y)-backcastMax+1):length(object$aux$y)]
#        OLD (fitted y-values); this was erroneous!:
#        meanFit[c(length(meanFit)-backcastMax+1):length(meanFit)]
    }

    ##predict:
    for(i in c(backcastMax+1):yhat.n){
      yhat[i] <- mconst + sum(arCoefs*yhat[c(i-1):c(i-arMax)]) + mxreghat[i]
    }

    ##out:
    outMean <- yhat[c(yhat.n-n.ahead+1):yhat.n]
    #outMean <- as.zoo(outMean)

  } #end mean spec


  ##--------------------
  ## 2 if variance spec
  ##--------------------

  outVariance <- NULL
  if(spec=="variance" || spec=="both"){ # || !is.null(plot.options$errors.only)

    ##record coef estimates:
    coefs <- as.numeric(coef.arx(object, spec="variance"))
    Elnz2hat <- coefs[length(coefs)]
    coefs <- coefs[-length(coefs)]

    ##vc:
    vconst <- as.numeric(coefs[1])

    ##arch:
    archMax <- 0
    archIndx <- 1
    if(!is.null(object$call$arch)){
      archEval <- eval(object$call$arch)
      archIndx <- 1:length(archEval) + 1
      archMax <- max(archEval)
      archCoefs <- rep(0,archMax)
      archCoefs[archEval] <- as.numeric(coefs[archIndx])
    }

    ##asym:
    asymMax <- 0
    asymIndx <- max(archIndx)
    if(!is.null(object$call$asym)){
      asymEval <- eval(object$call$asym)
      asymIndx <- 1:length(asymEval) + max(archIndx)
      asymMax <- max(asymEval)
      asymCoefs <- rep(0,asymMax)
      asymCoefs[asymEval] <- as.numeric(coefs[asymIndx])
    }

    ##log.ewma:
    logewmaMax <- 0
    logewmaIndx <- max(asymIndx)
    if(!is.null(object$call$log.ewma)){
      logewmaEval <- eval(object$call$log.ewma)
      if(is.list(logewmaEval)){ logewmaEval <- logewmaEval$length }
      logewmaIndx <- 1:length(logewmaEval) + max(asymIndx)
      logewmaMax <- max(logewmaEval)
      logewmaCoefs <- as.numeric(coefs[logewmaIndx])
    }

    ##backcast length:
    backcastMax <- max(archMax,asymMax,logewmaMax)

    ##vxreg:
    vxreghat <- rep(0, n.ahead + backcastMax)
    if(!is.null(object$call$vxreg)){

      ##check newvxreg:
      if(is.null(newvxreg)){ stop("'newvxreg' is NULL") }
      if(NROW(newvxreg)!=n.ahead){ stop("NROW(newvxreg) must equal n.ahead") }

      ##newmxreg:
      newvxreg <- coredata(cbind(as.zoo(newvxreg)))
      colnames(newvxreg) <- NULL

      ##vxreghat:
      vxregIndx <- c(max(logewmaIndx)+1):length(coefs)
      vxreghat <-  newvxreg %*% coefs[vxregIndx]
      vxreghat <- c(rep(0,backcastMax),vxreghat)

    } #end vxreg

    ##prepare lnsd2:
    lnsd2hat <- rep(NA, n.ahead + backcastMax)
    lnsd2hat.n <- length(lnsd2hat)
    lnsd2Fit <- log(coredata(fitted.arx(object, spec="variance")))
    lnsd2hat[1:backcastMax] <- lnsd2Fit[c(length(lnsd2Fit)-backcastMax+1):length(lnsd2Fit)]
    mLnsd2Hat <- matrix(NA, lnsd2hat.n, n.sim)
    mLnsd2Hat[,1:NCOL(mLnsd2Hat)] <- lnsd2hat

    ##prepare lnz2:
    if(backcastMax>0){
      lnz2hat <- rep(NA, n.ahead + backcastMax)
      lnz2hat.n <- length(lnz2hat)
      lnz2Fit <- coredata(object$ustar.residuals + Elnz2hat)
      lnz2hat[1:backcastMax] <- lnz2Fit[c(length(lnz2Fit)-backcastMax+1):length(lnz2Fit)]
      mLnz2Hat <- matrix(NA, lnz2hat.n, n.sim)
      mLnz2Hat[,1:NCOL(mLnz2Hat)] <- lnz2hat
    }

    ##zhat:
    ##-----

    ##bootstrap innov (not user-provided):
    if(is.null(innov)){
      zhat <- coredata(na.trim(object$std.residuals))
      where.zeros <- which(zhat==0)
      if(length(where.zeros)>0){ zhat <- zhat[-where.zeros] }
      draws <- runif(n.ahead*n.sim, min=0.5+.Machine$double.eps,
                     max=length(zhat)+0.5+.Machine$double.eps)
      draws <- round(draws, digits=0)
      zhat <- zhat[draws]
    }

    ##user-provided innov:
    if(!is.null(innov)){
      if(length(innov)!=n.ahead*n.sim){ stop("length(innov) must equal n.ahead*n.sim") }
      if(any(innov==0)){ stop("innov cannot contain zeros") }
      zhat <- as.numeric(innov)
    }

    zhat <- matrix(zhat,n.ahead,n.sim)
    mZhat2 <- zhat^2
    mLnz2Hat[c(backcastMax+1):NROW(mLnz2Hat),] <- log(mZhat2)

    vEpsilon2 <- rep(NA, n.ahead+backcastMax)
    vEpsilon2[1:backcastMax] <- as.numeric(object$residuals[c(length(object$residuals)-backcastMax+1):length(object$residuals)]^2)
    mZhat2 <- rbind(matrix(NA,backcastMax,NCOL(mZhat2)),mZhat2)


    ##prepare asym:
    if(asymMax>0){
      zhatIneg <- rep(NA, n.ahead + backcastMax)
      zhatIneg.n <- length(zhatIneg)
      zhatFit <- coredata(object$std.residuals)
      zhatIneg[1:backcastMax] <- zhatFit[c(length(zhatFit)-backcastMax+1):length(zhatFit)]
      zhatIneg <- as.numeric(zhatIneg<0)
      mZhatIneg <- matrix(NA, zhatIneg.n, n.sim)
      mZhatIneg[,1:NCOL(mZhatIneg)] <- zhatIneg
      mZhatIneg[c(backcastMax+1):NROW(mZhatIneg),] <- matrix(as.numeric(zhat<0),NROW(zhat),NCOL(zhat))
    }

    ##prepare log.ewma:
    if(logewmaMax>0){
      mLogEwmaHat <- matrix(NA, n.ahead+backcastMax, length(logewmaCoefs))
      colnames(mLogEwmaHat) <- object$aux$vXnames[logewmaIndx]
      mLogEwmaHat[1:backcastMax,] <- object$aux$vX[c(NROW(object$aux$vX)-backcastMax+1):NROW(object$aux$vX),logewmaIndx]
      mLogEwmaHat <- as.matrix(mLogEwmaHat)
    }

    ##predict:
    archTerm <- 0
    lnz2Term <- 0
    asymTerm <- 0
    logewmaTerm <- 0
    for(j in 1:NCOL(mLnsd2Hat)){
      for(i in c(backcastMax+1):NROW(mLnsd2Hat)){
        if(archMax>0){
          archTerm <- sum( archCoefs*mLnsd2Hat[c(i-1):c(i-archMax),j] )
          lnz2Term <- sum( archCoefs*mLnz2Hat[c(i-1):c(i-archMax),j] )
        }
        if(asymMax>0){
          asymTermSd2 <- sum( asymCoefs*mLnsd2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
          asymTermLnz2 <- sum( asymCoefs*mLnz2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
          asymTerm <- asymTermSd2 + asymTermLnz2
        }
        if(logewmaMax>0){
          for(k in 1:NCOL(mLogEwmaHat)){
            mLogEwmaHat[i,k] <- log( mean(vEpsilon2[c(i-logewmaEval[k]):c(i-1)]) )
          }
          logewmaTerm <- sum( coefs[logewmaIndx] * mLogEwmaHat[i,] )
        }
        mLnsd2Hat[i,j] <- vconst + archTerm + lnz2Term + asymTerm + logewmaTerm + vxreghat[i]
        vEpsilon2[i] <- exp(mLnsd2Hat[i,j])*mZhat2[i,j]
      } ##end for(i)
    } ##end for(j)

    ##out:
    outVariance <- mLnsd2Hat[c(lnsd2hat.n-n.ahead+1):lnsd2hat.n,]
    outVariance <- rowMeans(rbind(exp(outVariance)))

  } #end variance spec


  ##--------------------
  ## 3 out object
  ##--------------------

  ##out:
  if(spec=="mean"){ out <- outMean }
  if(spec=="variance"){ out <- outVariance }
  if(spec=="both"){
    out <- cbind(outMean,outVariance)
    colnames(out) <- c("mean","variance")
  }
  if(yAsRegular){
    startYear <- floor(as.numeric(index(yhataux))[1])
    startCycle <- cycle(yhataux)[1]
    out <- zooreg(out, start=c(startYear, startCycle),
      frequency=yFreq)
  }else{
    out <- zoo(out, order.by=newindex)
  }


  ##----------------
  ## 4 if plot=TRUE
  ##----------------

  ##determine plot-argument:
  plotArg <- plot
  if( is.null(plotArg) ){
    plotArg <- getOption("plot")
    if( is.null(plotArg) ){ plotArg <- FALSE }
  }

  ##if plot=TRUE:
  if(plotArg){

    ##how many observations to retain?:
    if(is.null(plot.options$keep)) {
      plot.options$keep=12L
    }

    ##include retained fitted?
    if(is.null(plot.options$fitted)) {
      plot.options$fitted <- FALSE
    }

#    ##WHAT DOES THIS DO???:
#    ##errors only?
#    if(!is.null(plot.options$errors.only)) {
##      plot.options$errors.only <- TRUE
#      if((spec=="mean")){
#        print("Cannot plot estimated error bars when mean specified")
#      }
#    }

    ##location of legend:
    if(is.null(plot.options$legend.loc)) {
      plot.options$legend.loc <- "topleft"
    }

    ##xlab:
    xlabArg <- ifelse(is.null(newindex), "Step ahead", "")

    ##make dataForPlot:
    dataForPlot <- matrix(NA, NROW(out), 6)
    colnames(dataForPlot) <- c("MeanActual", "MeanFitted",
      "MeanPrediction", "VarianceActual", "VarianceFitted",
      "VariancePrediction")
    if(!is.null(plot.options$newmactual)){
      dataForPlot[1:length(plot.options$newmactual),"MeanActual"] <-
        plot.options$newmactual
    }
    if(!is.null(outMean)){ dataForPlot[,"MeanPrediction"] <- outMean }
    if(!is.null(plot.options$newvactual)){
      dataForPlot[1:length(plot.options$newvactual),"VarianceActual"] <-
        plot.options$newvactual
    }
    if(is.null(outVariance)){
      dataForPlot[,"VariancePrediction"] <- sigma.arx(object)^2
    }else{
      dataForPlot[,"VariancePrediction"] <- outVariance
    }
    if(plot.options$keep > 0){
      retainedData <- matrix(NA, plot.options$keep, NCOL(dataForPlot))
      colnames(retainedData) <- c("MeanActual", "MeanFitted",
        "MeanPrediction", "VarianceActual", "VarianceFitted",
        "VariancePrediction")
      retainedData[,"MeanActual"] <-
        tail(coredata(yInSample), n=plot.options$keep)
      retainedData[,"MeanFitted"] <-
        tail(coredata(object$mean.fit), n=plot.options$keep)
      retainedData[,"VarianceActual"] <-
        tail(coredata(object$residuals^2), n=plot.options$keep)
      retainedData[,"VarianceFitted"] <-
        tail(coredata(object$var.fit), n=plot.options$keep)
      dataForPlot <- rbind(retainedData, dataForPlot)
      tmpIndx <- c(tail(index(yInSample), n=plot.options$keep),
        index(out))
      dataForPlot <- zoo(dataForPlot, order.by=tmpIndx)
    }else{
      dataForPlot <- zoo(dataForPlot, order.by=index(out))
    }

    ##get current par-values:
    def.par <- par(no.readonly=TRUE)

    ##prepare plotting:
    if(xlabArg=="") { #plot margins?
      par(mar = c(2,3,0.5,0.5) + 0.1) #b,l,t,r
    }else{
      par(mar = c(3,3,0.5,0.5) + 0.1) #b,l,t,r
    }
    plotTypeForecast <- ifelse(n.ahead==1, "p", "l")
    plotTypeRetained <- ifelse(plot.options$keep > 1, "l", "p")

    ##if spec="mean" or "both":
    ##-------------------------

    if( spec=="mean" || spec=="both" ){

      ##y-axis (limits):
      ylimVals <- range( na.trim(dataForPlot[,"MeanActual"]),
        na.trim(dataForPlot[,"MeanFitted"]),
        na.trim(dataForPlot[,"MeanPrediction"]
          + 2*sqrt(dataForPlot[,"VariancePrediction"])),
        na.trim(dataForPlot[,"MeanPrediction"]
          - 2*sqrt(dataForPlot[,"VariancePrediction"])) )

      ##plot the mean-predictions:
      plot.zoo(dataForPlot[,"MeanPrediction"], xlab="", ylab="",
        col="red", type=plotTypeForecast, ylim=ylimVals)

      ##add 2 x SEs:
      SE <- sqrt(dataForPlot[,"VariancePrediction"])
      lines(dataForPlot[,"MeanPrediction"] + 2*SE,
        col="darkgreen", lty=2, type=plotTypeForecast)
      lines(dataForPlot[,"MeanPrediction"] - 2*SE,
        col="darkgreen", lty=2, type=plotTypeForecast)

      ##add actuals (pre- and post-prediction):
      if(plot.options$keep > 0){
        lines(dataForPlot[,"MeanActual"], col="blue",
          type=plotTypeRetained)
      }

      ##add fitted (pre-prediction):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"MeanFitted"], col="red",
          type=plotTypeRetained)
      }

      ##add text closer to plot than xlab or ylab would do
      if(xlabArg!="") {
        mtext(xlabArg, side=1, line=2)
      }
      mtext("Mean", side=2, line=2)

      ##add legend:
      legend(plot.options$legend.loc, col=c("red","darkgreen","blue"),
        lty=c(1,2,1),
        legend=c("Forecast","+/- 2 x SE of regression","Actual"),
        bty="n")

    } #close if spec="mean" or "both"


    ##if spec="variance":
    ##-------------------

    if(spec=="variance"){

      ##y-axis (limits):
      ylimVals <- range( na.trim(dataForPlot[,"VarianceActual"]),
        na.trim(dataForPlot[,"VarianceFitted"]),
        na.trim(dataForPlot[,"VariancePrediction"]) )

      ##plot the variance-predictions:
      plot.zoo(dataForPlot[,"VariancePrediction"], xlab="", ylab="",
        col="red", type=plotTypeForecast, ylim=ylimVals)

      ##add fitted (pre-prediction):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"VarianceFitted"], col="red",
          type=plotTypeRetained)
      }

      ##add actuals (pre- and post-prediction):
      if(plot.options$keep > 0){
        lines(dataForPlot[,"VarianceActual"], col="blue",
          type=plotTypeRetained)
      }

      ##add text closer to plot than xlab or ylab would do
      if(xlabArg!="") {
        mtext(xlabArg, side=1, line=2)
      }
      mtext("Variance", side=2, line=2)

      ##add legend:
      legend(plot.options$legend.loc, col=c("red","blue"),
        lty=c(1,1), legend=c("Forecast (conditional variance)", "Actual (residuals^2)"), bty="n")

    } #close if spec="variance"


    ##add vertical green line between last observation
    ##and first forecast:
    if(plot.options$keep > 0){
      abline( v=( as.numeric(index(dataForPlot))[plot.options$keep] + as.numeric(index(dataForPlot))[plot.options$keep+1] )/2,
        col="darkgreen", lty=3)
    }

    ##return to old par-values:
    par(def.par)

  } #end if(plot)


  ##------------------
  ## 5 if return=TRUE
  ##------------------

  if(return){ return(out) }

} ##close predict.arx

##==================================================
## print estimation result
print.arx <- function(x, signif.stars=FALSE, ...)
{
  ##check if mean and variance have been fitted:
  xNames <- names(x)
  meanResults <- ifelse("mean.results" %in% xNames, TRUE, FALSE)
  varianceResults <- ifelse("variance.results" %in% xNames, TRUE, FALSE)

  ##header - first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(meanResults || varianceResults){
    estType <- ifelse(is.null(x$aux$user.estimator),
      "Ordinary Least Squares (OLS)", "User defined")
    cat("Dependent var.:", x$aux$y.name, "\n")
    cat("Method:", estType, "\n")
  }

  ##header - if mean results:
  if(meanResults){
    if(is.null(x$aux$user.estimator)){
      cat("Variance-Covariance:", switch(x$aux$vcov.type,
        ordinary = "Ordinary", white = "White (1980)",
        "newey-west" = "Newey and West (1987)"), "\n")
    }
    if("residuals" %in% xNames){
      cat("No. of observations (mean eq.):",
        length(na.trim(x$residuals)), "\n")
    }
  }

  ##header - if variance results:
  if( varianceResults && "resids.std" %in% xNames ){
    cat("No. of observations (variance eq.):",
      length(na.trim(x$std.residuals)), "\n")
  }

  ##header - sample info:
  if( "residuals" %in% xNames ){
    indexTrimmed <- index(na.trim(x$residuals))
    isRegular <- is.regular(x$residuals, strict=TRUE)
    isCyclical <- frequency(x$residuals) > 1
    if(isRegular && isCyclical){
      cycleTrimmed <- cycle(na.trim(x$residuals))
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
  } #end if( "residuals" %in% xNames )

  ##print mean results:
  if(meanResults){
    cat("\n")
    cat("Mean equation:\n")
    cat("\n")
    printCoefmat(x$mean.results, signif.stars=signif.stars)
  }

  ##print variance results:
  if(varianceResults){
    cat("\n")
    cat("Log-variance equation:\n")
    cat("\n")
    printCoefmat(x$variance.results, signif.stars=signif.stars)
  }

  ##print if no results:
  if( !meanResults && !varianceResults ){
    cat("\n")
    cat("   No estimation results\n")
  }

  ##goodness-of-fit:
  if( !"gof" %in% xNames && is.null(x$aux$user.estimator) ){
    gof <- matrix(NA, 3, 1)
    rownames(gof) <- c("SE of regression", "R-squared",
      paste("Log-lik.(n=", length(na.trim(x$std.residuals)), ")", sep=""))
    colnames(gof) <- ""
    gof[1,1] <- sigma.arx(x)
    gof[2,1] <- rsquared(x)
    gof[3,1] <- as.numeric(logLik.arx(x))
    x$gof <- gof
  }

  ##print diagnostics and fit:
  if( !is.null(x$diagnostics) ) {
    cat("\n")
    cat("Diagnostics and fit:\n")
    cat("\n")
    printCoefmat(x$diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
    if( !is.null(x$gof) ){printCoefmat(x$gof, digits=6, signif.stars=FALSE) }
  }

} #end print.arx

##==================================================
recursive <- function(object, spec=c("mean","variance"),
  std.errors=TRUE, from=40, tol=1e-07, LAPACK=FALSE,
  plot=NULL, return=TRUE)
{
  ##check if user-defined estimator:
  if( !is.null(object$aux$user.estimator) ){
    stop("Not available for user-defined estimators")
  }

  ##which specification:
  specType <- match.arg(spec)
#Change to?:
#  specType <- c("mean", "variance")
#  whichType <- charmatch(spec, specType)
#  specType <- specType[whichType]

  ##if mean-specification:
  if(specType=="mean"){
    if(is.null(object$mean.results)){
      stop("No mean-equation")
    }
    vY <- object$aux$y
    yNrow <- NROW(vY)
    mX <- object$aux$mX
    mXncol <- object$aux$mXncol
    mXnames <- object$aux$mXnames
    mainlab <- "Recursive estimates: mean equation"
  }

  ##if variance-specification:
  if(specType=="variance"){
    if(is.null(object$variance.results)){
      stop("No variance-equation")
    }
    vY <- object$aux$loge2
    yNrow <- NROW(vY)
    mX <- object$aux$vX
    mXncol <- object$aux$vXncol
    mXnames <- object$aux$vXnames
    mainlab <- "Recursive estimates: log-variance equation"
  }

  ##determine ols method:
  if(specType=="mean"){
    if(std.errors){
      if(object$aux$vcov.type=="ordinary"){ olsMethod=3 }
      if(object$aux$vcov.type=="white"){ olsMethod=4 }
      if(object$aux$vcov.type=="newey-west"){ olsMethod=5 }
    }else{
      olsMethod=1
    }
  } #close if(mean)

  if(specType=="variance"){
    if(std.errors){
      olsMethod=3
    }else{
      olsMethod=2
    }
  } #close if(variance)

  ##initialise:
  colnames(mX) <- mXnames
  recursiveEstimates <- matrix(NA, yNrow, mXncol)
  if(specType=="variance"){
    recursiveEstimatesElnz2 <- rep(NA, yNrow)
  }
  colnames(recursiveEstimates) <- mXnames
  if(std.errors){
    recursiveStdErrs <- recursiveEstimates
  }
  startIndx <- max(mXncol, min(from, yNrow))
  compute.at <- seq.int(from=yNrow, to=startIndx, by=-1)

  ##recursion:
  for(i in 1:length(compute.at)){

    ##estimate:
    vY <- vY[1:compute.at[i]]
    mXnames <- colnames(mX)
    NCOLmX <- NCOL(mX)
    mX <- dropvar(as.matrix(mX[1:compute.at[i], ]), tol=tol,
      LAPACK=LAPACK, silent=TRUE)
    if(NCOLmX==1){ colnames(mX) <- mXnames }
    tmpEst <- ols(vY, mX, tol=tol, LAPACK=LAPACK,
      method=olsMethod)
    recursiveEstimates[compute.at[i],colnames(mX)] <- tmpEst$coefficients
    if(std.errors){
      recursiveStdErrs[compute.at[i],colnames(mX)] <- sqrt(diag(tmpEst$vcov))
    }

    ##if variance-specification:
    if(specType=="variance"){
      Elnz2est <- -log(mean(exp(tmpEst$residuals)))
      recursiveEstimates[compute.at[i], "vconst"] <- recursiveEstimates[compute.at[i], "vconst"] - Elnz2est
      recursiveEstimatesElnz2[compute.at[i]] <- Elnz2est
    }

  } #close for loop

  ##rename std.errors columns:
  if(std.errors){
    colnames(recursiveStdErrs) <- paste(colnames(recursiveStdErrs),
      "SE", sep="")
  }

  ##set vconstSE to NA:
  if(std.errors==TRUE && specType=="variance"){
    recursiveStdErrs[,1] <- NA
  }

  ##handle zoo-index:
  naDiff <- object$aux$y.n - yNrow
  if(naDiff==0){
    zooIndx <- object$aux$y.index
  }else{
    zooIndx <- object$aux$y.index[-c(1:naDiff)]
  }
  recursiveEstimates <- zoo(recursiveEstimates,
    order.by=zooIndx)
  if(is.regular(recursiveEstimates, strict=TRUE)){ recursiveEstimates <- as.zooreg(recursiveEstimates) }
  if(std.errors){
    recursiveStdErrs <- zoo(recursiveStdErrs,
      order.by=zooIndx)
    if(is.regular(recursiveStdErrs, strict=TRUE)){ recursiveStdErrs <- as.zooreg(recursiveStdErrs) }
  }

  ##if return=TRUE:
  if(return){
    out <- list()
    out$estimates <- recursiveEstimates
    if(std.errors){
      out$standard.errors <- recursiveStdErrs
    }
  }

  ##plot:
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){
    recursiveEstimates <- na.trim(recursiveEstimates,
      is.na="all")
#    recursiveEstimates <- zoo(coredata(recursiveEstimates),
#      order.by=index(recursiveEstimates))
    plot(recursiveEstimates, main=mainlab, xlab="", col="blue")
  }

  ##out:
  if(return){
    return(out)
  }
} #close recursive

##==================================================
residuals.arx <- function(object, std=FALSE, ...)
{
  ##determine spec:
  if(is.null(std)){
    std <- switch(as.character(object$call)[1],
      arx=FALSE, getsm=FALSE, getsv=TRUE)
  }

  if(std){
    result <- object$std.residuals
  }else{
    result <- object$residuals
  }
  return(result)
} #end residuals.arx

##==================================================
## SE of regression
sigma.arx <- function(object, ...)
{
  residsTrimmed <- na.trim(object$residuals)
  RSS <- sum(residsTrimmed^2)
  nobs <- length(residsTrimmed)
  DFs <- length(coef.arx(object, spec="mean"))
  return( sqrt(RSS/(nobs-DFs)) )
} #close sigma.arx

##==================================================
## R-squared
rsquared <- function(object, ...)
{
  classOK <- class(object) %in% c("arx", "gets", "isat")
  if(!classOK){ stop("Object not of class 'arx', 'gets' nor 'isat'") }
  if( class(object) == "gets" ){
    specType <- switch(as.character(object$call)[1],
      getsm="mean", getsv="variance")
  }
  if( class(object) == "gets" && specType=="variance" ){
    Rsquared <- NA
  }else{
    TSS <- sum( (object$aux$y - mean(object$aux$y))^2 )
    residsTrimmed <- na.trim(object$residuals)
    RSS <- sum(residsTrimmed^2)
    Rsquared <- 1 - RSS/TSS
    ##to do: adjusted R-squared
  }
  return(Rsquared)
} #close rsquared function

##==================================================
## summarise output
summary.arx <- function(object, ...)
{
  object.name <- deparse(substitute(object))
  xnames <- cbind(names(object))
  xrows <- nrow(xnames)
  colnames(xnames) <- ""
  rownames(xnames) <- rep(" ",xrows)
  cat("\n")
  cat("Items in list '", object.name, "':", sep="")
  cat("\n")
  print(xnames, quote=FALSE)
  cat("\n")
} #end summary.arx

##==================================================
VaR <- function(object, level=0.99, type=7, ...)
{
  ##check whether class is valid:
  classType <- class(object)
  if( !classType %in% c("arx", "gets", "isat") ){
    stop("object not of class 'arx', 'gets' nor 'isat'")
  }

  ##check the risk-levels:
  riskLevel <- 1-level
  if( any(riskLevel > 1) || any(riskLevel < 0) ){
    stop("risk-level(s) must be in the 0 to 1 interval")
  }

  ##fitted mean and sd, standardised residuals, quantile:
  meanFit <- fitted(object, spec="mean")
  sdFit <- sqrt( fitted(object, spec="variance") )
  residsStd <- residuals(object, std=TRUE)
  qValue <- quantile(residsStd, probs=riskLevel, type=type,
    names=FALSE, na.rm=TRUE)
  if( length(riskLevel)==1 ){
    VaR <- meanFit + sdFit*qValue
  }else{
    colNames <- paste("VaR", level, sep="")
    VaR <- matrix(NA, length(sdFit), length(colNames))
    for(i in 1:length(colNames)){
      VaR[,i] <- meanFit + sdFit*qValue[i]
    }
    colnames(VaR) <- colNames
    VaR <- zoo(VaR, order.by=index(sdFit))
  }

  ##return
  return(-VaR)
} #close VaR

##==================================================
vcov.arx <- function(object, spec=NULL, ...)
{

  ##spec argument:
  if(is.null(spec)){
    if(!is.null(object$mean.results)){
      spec <- "mean"
    }
    if(is.null(object$mean.results)
      && !is.null(object$variance.results) ){
      spec <- "variance"
    }
  }else{
    spec.type <- c("mean", "variance")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  ##mean:
  if(spec=="mean"){
    result <- object$vcov.mean
  }

  ##variance:
  if(spec=="variance"){
    result <- object$vcov.var
  }

#  ##check and change if 0 x 0?:
#  if(all(dim(result)==0)){ result <- NULL }

  ##check if user-specified estimator:
  if( !is.null(object$aux$user.estimator) && is.null(result) ){
    result <- object$vcov
  }
  
  return(result)
} #end vcov.arx

####################################################
##3 GETS FUNCTIONS
####################################################

##==================================================
## Multi-path GETS modelling of mean specification
getsm <- function(object, t.pval=0.05, wald.pval=t.pval,
  vcov.type=NULL, do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc", "aic", "hq"),
  keep=NULL, include.gum=FALSE, include.empty=FALSE,
  max.paths=NULL, max.regs=NULL, zero.adj=NULL, vc.adj=NULL,
  verbose=TRUE, print.searchinfo=TRUE, estimate.specific=TRUE,
  plot=NULL, alarm=FALSE)
{
  ### ARGUMENTS: ###########

  info.method <- match.arg(info.method)

  ##determine no. of paths:
  if( !is.null(max.paths) && max.paths < 1){
    stop("'max.paths' cannot be smaller than 1")
  }

  ##determine vcov.type:
  if(is.null(vcov.type)){
    vcov.type <- object$aux$vcov.type
  }else{
    vcovTypes <- c("ordinary", "white", "newey-west")
    which.type <- charmatch(vcov.type, vcovTypes)
    vcov.type <- vcovTypes[which.type]
  }

  ##determine ar and arch lags:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])

  ## max.regs, tol, LAPACK:
  if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }
  tol <- object$aux$tol
  LAPACK <- object$aux$LAPACK

  ##check if mean equation:
  if( is.null(object$aux$mX) ){ stop("Mean equation empty") }

  ##check if variance equation:
  ##in the future, check object$aux$vX instead?:
  if( is.null(object$variance.results )){
    var.spec.chk <- FALSE
  }else{
    var.spec.chk <- TRUE
  }

  ## check if estimator user-defined:
  if( !is.null(object$aux$user.estimator) ){
    if(estimate.specific){
      estimate.specific <- FALSE
      message("  'estimate.specific' set to FALSE")
    }
    if( is.null(plot) || identical(plot,TRUE) ){
      plot <- FALSE
      message("  'plot' set to FALSE")
    }
  }

  ### INITIALISE ##########

  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  messages <- NULL
  spec <- list() #terminal specifications
  spec.results <- NULL

  ##deletable, non-deletable regressors, re-organise:
  keep.n <- length(keep)
  gum <- 1:object$aux$mXncol
  delete <- setdiff(gum, keep) #integer(0) if empty
  delete.n <- length(delete)
  if(delete.n > 0){mXdel <- cbind(object$aux$mX[,delete])}else{mXdel <- NULL}
  if(is.null(keep)){mXndel <- NULL}else{mXndel <- cbind(object$aux$mX[,keep])}
  mXadj <- cbind(mXdel,mXndel)

  ## add keep column to results:
  tmp <- rep(0,object$aux$mXncol)
  if(!is.null(keep)){ tmp[keep] <- 1 }
  tmpdf <- cbind(tmp, object$mean.results)
  tmp <- 1:object$aux$mXncol
  tmpdf <- cbind(tmp, tmpdf)
  colnames(tmpdf)[1:2] <- c("reg.no", "keep")
  out$gum.mean <- tmpdf
  out$gum.variance <- object$variance.results
  out$gum.diagnostics <- object$diagnostics

  ## GUM: ##################

  ##estimate GUM:
  if( is.null(object$aux$user.estimator) ){

    ##usual ols:
    estMethod <- which( vcov.type==c("none", "none", "ordinary",
      "white", "newey-west") )
    est <- ols(object$aux$y, mXadj, tol=object$aux$tol,
      LAPACK=object$aux$LAPACK, method=estMethod,
      user.fun=NULL, user.options=NULL)
    est$std.residuals <- coredata(na.trim(object$std.residuals))
    est$logl <- object$logl
    if( !is.null(object$aux$loge2.n) ){
      est$n <- object$aux$loge2.n
    }

  }else{

    ##user-defined estimator:
    est <- do.call(object$aux$user.estimator$name,
      list(object$aux$y, mXadj), envir=.GlobalEnv)
    #delete?:
    if( is.null(est$vcov) && !is.null(est$vcov.mean) ){
      est$vcov <- est$vcov.mean
    }

  } #end if( is.null(..) )

  ##diagnostics:
  if( !is.null(est$residuals) ){
    gum.chk <- diagnostics(est, ar.LjungB=ar.LjungB,
      arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
      verbose=FALSE, user.fun=user.diagnostics)
  }else{
    gum.chk <- TRUE
  }

  ## if GUM passes diagnostic checks:
  if(gum.chk){

    spec[[1]] <- spec.gum <- gum #add gum to list of terminal specs

    ## vcov, t-stats, p-vals:
    stderrs <- sqrt(diag(est$vcov))
    t.stat <- est$coefficients/stderrs
    p.val <- pt(abs(t.stat), est$df, lower.tail=FALSE)*2

    ##specification results
    info.results <- info.criterion(est$logl, est$n, est$k,
      method=info.method)
    spec.results <- rbind( c(info.results$value, est$logl,
      info.results$n, info.results$k) )
    col.labels <- c(paste("info(", info.method, ")", sep=""),
      "logl", "n", "k")
    row.labels <- c("spec 1 (gum):")

    ##record data for tests against gum:
    gum.regs <- c(delete, keep)
    gum.coefs <- object$mean.results[gum.regs,1]
    gum.varcovmat <- est$vcov #OLD: varcovmat

  }else{
    messages <- paste(messages,
      "- MGUM does not pass one or more diagnostic checks", sep="")
  }

  ## FUTURE: ADD 1-CUT MODEL #########

  ## EMPTY MODEL: ################

  if( gum.chk && delete.n>0 && include.empty ){

    ## DO NOT do pet in order to enable reality check!

    ## estimate model:
    if( is.null(object$aux$user.estimator) ){
      est <- ols(object$aux$y, mXndel, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, method=estMethod)
    }else{
      est <- do.call(object$aux$user.estimator$name,
        list(object$aux$y, mXndel), envir=.GlobalEnv)
      #delete?:
      if( is.null(est$vcov) && !is.null(est$vcov.mean) ){
        est$vcov <- est$vcov.mean
      }
    } #end if( is.null(..) )

    #add var.fit, modify resids, resids.std, logl and n:
    if(var.spec.chk){
      residsAdj <- zoo(est$residuals, order.by=object$aux$y.index)
      est.var <- arx(residsAdj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
        zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
        tol=object$aux$tol, LAPACK=object$aux$LAPACK, plot=FALSE)
      est$std.residuals <- coredata(na.trim(est.var$std.residuals))
      est$var.fit <- coredata(na.trim(est.var$var.fit))
      est$logl <- est.var$logl
      est$residuals <- est$residuals[c(object$aux$y.n-object$aux$loge2.n+1):object$aux$y.n]
      est$n <- length(est$residuals)
    }

    diagnostics.chk <- diagnostics(est, ar.LjungB=ar.LjungB,
      arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
      verbose=FALSE, user.fun=user.diagnostics)

    ##if empty model passes diagnostic checks:
    if(diagnostics.chk){

      ##add empty to spec:
      spec[[length(spec)+1]] <- if(is.null(keep)){0}else{keep}

      ##specification results
      info.results <- info.criterion(est$logl, est$n, est$k,
        method = info.method)
      spec.results <- rbind(spec.results,
        c(info.results$value, est$logl, info.results$n,
        info.results$k))
      row.labels <- c(row.labels,
        paste("spec ", length(spec), " (empty):", sep=""))

    }else{
        messages <- paste(messages,
          "- Empty mean model does not pass one or more diagnostic checks",
          sep="")
    } #end if(empty passes diagnostics==TRUE){..}else{..}

  } #end if(include empty model==TRUE)

## MULTI-PATH SEARCH: #################

#future:
#regs <- list()
#regs.info <- list()

#the following should probably be moved to just before
#"single path search", that is, just beneath "prepare
#single-path search":
#regs.current.path <- list() #the regressions of the current path
#regs.info.current.path <- list() #the regression info of current path
#regression info: list(which.path=??, where.in.path=??)

insig.regs <- NULL
paths <- list()
if( gum.chk && delete.n>1 ){

  ## no. of paths:
  insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
  if( !is.null(max.paths) ){
    if(max.paths < length(insig.regs)){
      pvalRanksInv <- rank( 1-p.val[insig.regs] )
      insig.regs <- delete[ pvalRanksInv <= max.paths ]
    }
  }
  n.paths <- length(insig.regs)

  ## if paths = 0:
  if(n.paths == 0){
    messages <- paste(messages,
      "- All regressors significant in GUM mean equation", sep="")
  }

  ## if paths > 0:
  if(n.paths > 0){

    if(print.searchinfo){
      message(n.paths, " path(s) to search")
      message("Searching: ", appendLF=FALSE)
    }

    ## paths:
    for(i in 1:n.paths){

      ## print searchinfo:
      if(print.searchinfo){
        newLine <- ifelse(i==n.paths, TRUE, FALSE)
        message(i, " ", appendLF=newLine)
      }

      ## prepare single-path search:
      path <- insig.regs[i]
      delete.adj <- setdiff(delete, insig.regs[i])
      keep.adj <- as.numeric(keep)

      ## single-path search of path i:
      for(j in 1:max.regs){

        ## matrices:
        mXdell <- if(length(delete.adj)==0){NULL}else{object$aux$mX[,delete.adj]}
        mXndell <- if(is.null(keep.adj)){NULL}else{object$aux$mX[,keep.adj]}
        mXadj <- cbind(mXdell,mXndell)
        mXadj.k <- NCOL(mXadj)

        ## estimate model:
        if( is.null(object$aux$user.estimator) ){
          est <- ols(object$aux$y, mXadj, tol=object$aux$tol,
            LAPACK=object$aux$LAPACK, method=estMethod,
            user.fun=NULL, user.options=NULL)
        }else{
          ##user-defined estimator:
          est <- do.call(object$aux$user.estimator$name,
            list(object$aux$y, mXadj), envir=.GlobalEnv)
          #delete?:
          if( is.null(est$vcov) && !is.null(est$vcov.mean) ){
            est$vcov <- est$vcov.mean
          }
        } #end if(is.null(..))else(..)

        #add var.fit, modify resids, resids.std, logl and n:
        if(var.spec.chk){
          residsAdj <- zoo(est$residuals, order.by=object$aux$y.index)
          est.var <- arx(residsAdj, vc=object$aux$vc,
            arch=object$aux$arch, asym=object$aux$asym,
            log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
            zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
            tol=object$aux$tol, LAPACK=object$aux$LAPACK, plot=FALSE)
          est$std.residuals <- coredata(na.trim(est.var$std.residuals))
          est$var.fit <- coredata(na.trim(est.var$var.fit))
          est$logl <- est.var$logl
          est$residuals <- est$residuals[c(object$aux$y.n-object$aux$loge2.n+1):object$aux$y.n]
          est$n <- length(est$residuals)
        }

        ##diagnostics:
        diagnostics.chk <- diagnostics(est, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          verbose=FALSE, user.fun=user.diagnostics)

        ## if diagnostics fails (FALSE),
        ## then move path[length(path)] to keep.adj
        if(!diagnostics.chk){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*c(-1))
          next #next j
        }

        #if diagnostics are ok (TRUE):
        if(diagnostics.chk){

          ## stop if no more deletable regressors:
          if( length(delete.adj)==0 ){
            spec.adj <- keep.adj
            break
          } #end if(length(..)==0)

          ## vcov, t-stats, p-vals:
#for the future:
#          if( is.null(est$vcov) ){
#            est$vcov <- vcovFun(est, method="ordinary")
#          }
          stderrs <- sqrt(diag(est$vcov))
          t.stat <- est$coefficients/stderrs
          p.val <- pt(abs(t.stat), est$df, lower.tail=FALSE)*2

# for new version: here is where I need to identify the coefficient,
# among those not in keep, with highest p-value. sketch?:
# gumRegs # regressors in gum, e.g. 1:NCOL(x)
# keep # the regressors to be kept, e.g. 1:2
# keepAdj # keep + those put in keep
# regsAdj # regressors left
# pvalRanks <- rank(p.val)
# pvalRanksAdj <- pvalRanks[-keepAdj]
# pvalMin <- which.min(pvalRanksAdj)
# reg.no <- keepAdj[ which(pvalRanks == pvalMin) ]

          ## if any p-value > t.pval:
          if(sum(p.val[1:c(length(delete.adj))] > t.pval) > 0){

            reg.no <- which.max(p.val[1:I(length(delete.adj))])

            ## do pet test (i.e. wald-test):
            if(do.pet){
              deleted <- setdiff(delete, delete.adj[-reg.no])
              n.deleted <- length(deleted)
#for the new version:
#              mR <- matrix(0, length(gum.regs), length(gum.regs))
#              diag(mR) <- 1
#              mR <- rbind(mR[deleted,])
              mR <- NULL #initiate restriction matrix
              for(k in 1:length(gum.regs)){
              #for(k in 1:object$aux$mXncol){
                if(gum.regs[k] %in% deleted){
                  mR <- rbind(mR, c(rep(0,I(k-1)), 1, rep(0, I(object$aux$mXncol-k) )))
                } #end if(gum.regs[k}..)
              } #end for(k in ..)

              mRestq <- mR%*%cbind(gum.coefs)
              wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=object$aux$tol)%*%mRestq
              pet.chk <- as.logical(wald.pval < pchisq(wald.stat, n.deleted, lower.tail = FALSE))
            }else{
              pet.chk <- TRUE
            } #end if(do.pet)else..

            ## delete regressor if(pet.chk), else move to keep:
            if(pet.chk){
              path <- union(path, delete.adj[reg.no])
              delete.adj <- delete.adj[-reg.no]
            }else{
              path <- union(path, delete.adj[reg.no]*I(-1))
              keep.adj <- union(delete.adj[reg.no], keep.adj)
              delete.adj <- delete.adj[-reg.no]
            } #end if(pet.chk)else{..}

          }else{
            spec.adj <- union(delete.adj, keep.adj)
            break
          } #end if..else.. any p-value > t.pval

        } #end if diagnostics are ok

      } #### end single-path search: for(j in..

      #it is probably at this point that I should introduce
      #the 'bookkeeping' with respect to regs, regs.info,
      #regs.current.path and regs.info.current.path

      #add path to the paths list:
      paths[[length(paths)+1]] <- path

      #check if spec.adj (terminal) is already in spec:
      if(length(spec.adj)==0){spec.adj <- 0} #check if completely empty
      for(l in 1:length(spec)){
        chk.spec <- setequal(spec.adj, spec[[l]])
        if(chk.spec==TRUE){break} #stop for(l in..)
      }

      #if spec.adj not in spec (among terminals):
      if(chk.spec==FALSE){

        #add spec.adj to spec:
        spec[[length(spec)+1]] <- spec.adj

        #specification results
        if(spec.adj[1]==0){
          n.spec.adj <- 0
        }else{
          n.spec.adj <- length(spec.adj)
        }
#        if( is.null(object$aux$user.estimator) ){
#          est$logl <- -resids.adj.n*log(2*pi)/2 - sum(log(est$var.fit))/2 - sum(est$std.residuals^2)/2
#        }
        info.results <- info.criterion(est$logl, est$n,
          n.spec.adj, method=info.method)

        #add terminal to spec.results:
        spec.results <- rbind(spec.results,
          c(info.results$value, est$logl, info.results$n,
          info.results$k))
        row.labels <- c(row.labels, paste("spec ", length(spec), ":", sep=""))

      } #end if(chk.spec==FALSE)

    } #end multi-path search: for(i in 1:n.paths) loop

  } #end if paths > 0
} #end if( gum.chk!=0 && delete.n>1 )

  ## FIND THE BEST MODEL: ########################

  #future?: if include.gum=FALSE (default), then first check if
  #spec results is empty, then add gum to it if it is

  if(!is.null(spec.results)){

    J <- 1:NROW(spec.results)
    models <- cbind(J, spec.results)
    colnames(models) <- NULL

    #find best model and check for several minimums:
    if(include.gum){
      min.value <- min(models[,2])
      where <- which(min.value==models[,2])
    }else{
      if(length(spec)==1){
        where <- 1
      }else{
        min.value <- min(models[-1,2])
        where <- which(min.value==models[-1,2]) + 1
      } #end if(length(spec)==1)
    } #end if(include.gum)..
    if(length(where)>1){
      messages <- paste(messages,
        "- Several terminal specifications attain the minimum information criterion",
        sep="")
      }
    best.spec <- spec[[where[1]]] #winner

  } #end if(!is.null(spec.results))

  ## OUTPUT ################################

  out$keep <- keep
  #out$insigs.in.gum <- insig.regs #don't think I need this one...

  ##if no search has been undertaken:
  if(is.null(spec.results)){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }

  ##if search has been undertaken:
  if(!is.null(spec.results)){

    ##terminals results:
    if(length(paths)==0){
      out$paths <- NULL
    }else{ out$paths <- paths }
    out$terminals <- spec
    colnames(spec.results) <- col.labels
    where.empty <- which(spec.results[,"k"]==0)
    if(include.empty==FALSE && length(where.empty) > 0){
      row.labels[where.empty] <- paste("spec ", where.empty,
        " (empty):", sep="")
    }
    rownames(spec.results) <- row.labels
    out$terminals.results <- spec.results

    if(!estimate.specific){
      if(best.spec==0 || is.na(best.spec) || length(best.spec)==0 ){
        out$specific.spec <- NULL
      }else{
        specific <- sort(best.spec)
        names(specific) <- object$aux$mXnames[specific]
        out$specific.spec <- specific
      }
    } #end if(!estimate.specific)

    if(estimate.specific){

      ##prepare for estimation:
      yadj <- zoo(cbind(object$aux$y),
        order.by=object$aux$y.index)
      colnames(yadj) <- object$aux$y.name
      specific <- sort(best.spec)
      if(specific[1]==0){
        mXadj <- NULL
      }else{
        mXadj <- cbind(object$aux$mX[,specific])
        colnames(mXadj) <- object$aux$mXnames[specific]
        mXadj <- zoo(mXadj, order.by=object$aux$y.index)
      }
      if(is.null(object$aux$vxreg)){
        vxregAdj <- NULL
      }else{
        vxregAdj <- zoo(object$aux$vxreg,
          order.by=object$aux$y.index)
      }
      if(is.null(ar.LjungB)){
        ar.LjungB <- object$aux$qstat.options[1]
      }
      if(is.null(arch.LjungB)){
        arch.LjungB <- object$aux$qstat.options[2]
      }

      ##estimate specific model:
      est <- arx(yadj, mxreg=mXadj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=vxregAdj,
        zero.adj=object$aux$zero.adj,
        vc.adj=object$aux$vc.adj, vcov.type=vcov.type,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        user.diagnostics=user.diagnostics, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, plot=FALSE)

      ##rename various stuff:
      est$call <- est$date <- NULL
      where.diagnostics <- which(names(est)=="diagnostics")
      if(length(where.diagnostics)>0){
        names(est)[where.diagnostics] <- "specific.diagnostics"
      }
      est$aux$y.name <- object$aux$y.name
      est <- unclass(est)
      names(specific) <- colnames(mXadj)
      out$specific.spec <- specific
      out <- c(out,est)
      ##solution to issue raised by J-dog?:
#      if(specific[1]==0){
#        out$aux$mX <- out$aux$mXnames <- NULL
#      }

    } #end if(estimate.specific)

  } #end if(!is.null(spec.results))

  if(!is.null(messages)){ out$messages <- messages }
  if( is.null(object$aux$user.estimator) ){
    out$aux$y <- object$aux$y
    out$aux$y.index <- object$aux$y.index
    out$aux$y.n <- object$aux$y.n
    out$aux$y.name <- object$aux$y.name
    out$aux$mXnames.gum <- object$aux$mXnames
    out$aux$call.gum <- object$call
    if(is.null(out$aux$vcov.type)){ out$aux$vcov.type <- vcov.type }
  }
  out <- c(list(date=date(), gets.type="getsm"), out)
  out$time.finished <- date()
  class(out) <- "gets"

  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.gets(out) }
  return(out)
} #close getsm function

##==================================================
## Multi-path GETS modelling of log-variance
getsv <- function(object, t.pval=0.05, wald.pval=t.pval,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025),
  normality.JarqueB=NULL, user.diagnostics=NULL,
  info.method=c("sc", "aic", "hq"), keep=c(1), include.gum=FALSE,
  include.1cut=FALSE, include.empty=FALSE, max.paths=NULL,
  turbo=FALSE, max.regs=NULL, zero.adj=NULL, vc.adj=NULL,
  print.searchinfo=TRUE, estimate.specific=TRUE, plot=NULL, alarm=FALSE)
{
  ### ARGUMENTS ###########

  info.method <- match.arg(info.method)
  vc=TRUE #obligatory
  vcov.type <- "ordinary" #obligatory
  tol <- object$aux$tol
  LAPACK <- object$aux$LAPACK

  ##zoo and NA related:
  e <- object$residuals #should not contain NAs
  e.index <- index(e) #use object$aux$y.index instead?
  e <- coredata(e)
  e.n <- length(e) #use object$aux$y.n instead?
  eadj <- e[c(e.n-object$aux$loge2.n+1):e.n] #Note: log(eadj^2)=loge2
  eadj.n <- length(eadj)
  eadj.index <- e.index[c(e.n-object$aux$loge2.n+1):e.n]

  ##diagnostics options, max.regs:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
  if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }

  ##zero-handling:
  if(is.null(zero.adj)){ zero.adj <- object$aux$zero.adj }
  if(is.null(vc.adj)){ vc.adj <- object$aux$vc.adj }

  ##tolerancy for non-invertibility of matrices:
  if(is.null(tol)){ tol <- object$aux$tol }
  if(is.null(LAPACK)){ LAPACK <- object$aux$LAPACK }


  ### INITIALISE ##########

  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  loge2 <- object$aux$loge2
  mX <- object$aux$vX
  colnames(mX) <- object$aux$vXnames
  if( !(1 %in% keep) ){
    keep <- union(1,keep)
    warning("Regressor 1 set to 'keep'")
  }

  ##add gum results and diagnostics to out:
  out$gum.mean <- object$mean.results
  tmp <- matrix(0, NROW(object$variance.results), 2)
  colnames(tmp) <- c("reg.no.", "keep")
  tmp[,1] <- 1:NROW(tmp) #fill reg.no. column
  tmp[keep,2] <- 1 #fill keep column
  out$gum.variance <- cbind(tmp, object$variance.results)
  out$gum.diagnostics <- object$diagnostics


  ### DO MULTI-PATH GETS ##########

  ##do the gets:
  est <- getsFun(loge2, mX, user.estimator=list(name="ols",
    untransformed.residuals=eadj, tol=tol, LAPACK=LAPACK, method=6),
    gum.result=NULL, t.pval=t.pval, wald.pval=wald.pval, do.pet=do.pet,
    ar.LjungB=ar.LjungB, arch.LjungB=arch.LjungB,
    normality.JarqueB=normality.JarqueB, user.diagnostics=user.diagnostics,
    gof.function=list(name="infocrit", method=info.method),
    gof.method="min", keep=keep, include.gum=include.gum,
    include.1cut=include.1cut, include.empty=include.empty,
    max.paths=max.paths, turbo=turbo, max.regs=max.regs,
    print.searchinfo=print.searchinfo, alarm=alarm)
  est$time.started <- NULL
  est$time.finished <- NULL
  est$call <- NULL
  out <- c(out, est)

  ## if no search has been undertaken:
  if(is.null(est$terminals.results)){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }


  ### ESTIMATE SPECIFIC ################

  if(estimate.specific){

    ## prepare estimation:
    e <- zoo(cbind(eadj), order.by=eadj.index)
    colnames(e) <- "e"
    specificadj <- setdiff(out$specific.spec, 1)
    if(length(specificadj)==0){
      vXadj <- NULL
    }else{
      vXadj <- cbind(object$aux$vX[,specificadj])
      colnames(vXadj) <- object$aux$vXnames[specificadj]
      vXadj <- zoo(vXadj, order.by=eadj.index)
    }
    if(is.null(ar.LjungB)){
      ar.LjungB <- object$aux$qstat.options[1]
    }
    if(is.null(arch.LjungB)){
      arch.LjungB <- object$aux$qstat.options[2]
    }

    ## estimate model:
    est <- arx(e, vc=TRUE, vxreg=vXadj,
      zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
      qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
      user.diagnostics=user.diagnostics, tol=object$aux$tol,
      LAPACK=object$aux$LAPACK, plot=FALSE)

    ## delete, rename and change various stuff:
    est$call <- est$date <- NULL
    where.diagnostics <- which(names(est)=="diagnostics")
    if(length(where.diagnostics)>0){
      names(est)[where.diagnostics] <- "specific.diagnostics"
    }
    est$mean.fit <- object$mean.fit[ index(object$mean.fit) %in% eadj.index ]
    #est$mean.fit <- object$mean.fit[ eadj.index ] #should work, but doesn't!
    est$aux$vxreg <- est$aux$vxreg.index <- NULL
    est$aux$y.name <- "e"

    ## finalise:
    est <- unclass(est)
    out <- c(out,est)

  } #end if(estimate.specific)

  ### OUTPUT ########

  out$aux$vXnames.gum <- object$aux$vXnames
  out$aux$call.gum <- object$call
  if(is.null(out$aux$vcov.type)){ out$aux$vcov.type <- vcov.type }
  out <- c(list(date=date(), gets.type="getsv"), out)
  out$time.finished <- date()
  class(out) <- "gets"

  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.gets(out) }
  return(out)
} #close getsv function


##==================================================
coef.gets <- function(object, spec=NULL, ...)
{
  coef.arx(object, spec=spec)
} #end coef.gets

##==================================================
## fitted values
fitted.gets <- function(object, spec=NULL, ...)
{
  fitted.arx(object, spec=spec)
} #end fitted.gets

##==================================================
logLik.gets <- function(object, ...)
{
  logLik.arx(object)
} #end logLik.gets

##==================================================
## extract paths
paths <- function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$paths)
  }else{
    stop("The object does not belong to the 'gets' nor 'isat' class")
  }
} #end paths

##==================================================
## plot gets object
plot.gets <- function(x, spec=NULL, col=c("red","blue"),
  lty=c("solid","solid"), lwd=c(1,1), ...)
{
  plot.arx(x, spec=spec, col=col, lty=lty, lwd=lwd)
} #close plot.gets


##==================================================
## forecast up to n.ahead
predict.gets <- function(object, spec=NULL, n.ahead=12,
  newmxreg=NULL, newvxreg=NULL, newindex=NULL,
  n.sim=1000, innov=NULL, return=TRUE, plot=NULL,
  plot.options=list(), ...)
{
  ##create new object to add stuff to in order to use predict.arx()
  object.new <- object

  ##spec:
  if(is.null(spec)){
    if(object$gets.type=="getsm"){ spec <- "mean" }
    if(object$gets.type=="getsv"){ spec <- "variance" }
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  } #end if(..)else(..)
  if( spec %in% c("variance","both") ){ stop("variance-forecasts with predict.gets not available yet") }


  ##-----------------------------------
  ## arguments mean-equation:
  ##-----------------------------------

  if("mX" %in% names(object$aux)) {
    ##what dynamics specified in gum?
    gum.ar <- eval(object$aux$call.gum$ar)
    ##what dynamics remain in spec?
    spec.ar <- as.numeric(gsub("ar(\\d+)","\\1",object$aux$mXnames[grep("^ar\\d+$",object$aux$mXnames)]))

    if(NROW(spec.ar)==0) {
      object.new$call$ar <- NULL
    } else { ##check that dynamics in specific are subset of those in gum
      object.new$call$ar <- spec.ar[spec.ar %in% gum.ar]
    }

    ##need to check for constant
    if(!is.null(object$aux$call.gum$mc)) {
      ##but the constant needs to have been retained...
      if ( "mconst" %in% object$aux$mXnames ) { object.new$call$mc <- TRUE }
    }

    if(is.null(object$aux$call.gum$mc) && "mconst" %in% object$aux$mXnames) {
      ##ensure that when we call predict.arx we set mc=TRUE
      object.new$call$mc <- !is.null(object$aux$call.gum$mc)
      ##need to also get rid of mconst in mX
      #object.new$mX <- object$aux$mX[,-grep("mconst",object$aux$mXnames)]
    }

    ##mxreg needs to be null if only thing in it is dynamics (i.e. ar variables)
    if(NROW(grep("^ar\\d+$",object$aux$mXnames))==object$aux$mXncol) {
      object.new$call$mxreg <- NULL
    } else {
      object.new$call$mxreg <- object$aux$call.gum$mxreg
    } ##if there is non-arch variables in mXreg, user needs to specify newmxreg

  } else {

    object.new$call$mc <- NULL
    object.new$call$ar <- NULL
    ##to do: log.ewma
    object.new$call$mxreg <- NULL

  }

  ##-----------------------------------
  ## arguments variance-equation:
  ##-----------------------------------

#  if("vX" %in% names(object$aux)) {
#
#    ##what dynamics specified in gum?
#    gum.arch <- eval(object$aux$call.gum$arch)
#    ##what dynamics remain in spec?
#    spec.arch <- as.numeric(gsub("arch(\\d+)","\\1",object$aux$vXnames[grep("^arch\\d+$",object$aux$vXnames)]))
#    if(NROW(spec.arch)==0) {
#      object.new$call$arch <- NULL
#    } else { ##check that dynamics in specific are subset of those in gum
#      object.new$call$arch <- spec.arch[spec.arch %in% gum.arch]
#    }
#
#    ##which asyms specified in gum?
#    gum.asym <- eval(object$aux$call.gum$asym)
#    ##which asyms remain in spec?
#    spec.asym <- as.numeric(gsub("asym(\\d+)","\\1",object$aux$vXnames[grep("^asym\\d+$",object$aux$vXnames)]))
#    if(NROW(spec.asym)==0) {
#      object.new$call$asym <- NULL
#    } else { ##check that dynamics in specific are subset of those in gum
#      object.new$call$asym <- spec.asym[spec.asym %in% gum.asym]
#    }
#
#    ##which asyms specified in gum?
#    gum.log.ewma <- eval(object$aux$call.gum$log.ewma)
#    if(!is.null(gum.log.ewma)){ stop("'log.ewma' not available yet for 'predict.gets' ") }
#
#    ##vxreg needs to be null if only thing in it is dynamics (i.e. arch variables)
#    if(NROW(grep("^arch\\d+$",object$aux$vXnames))==object$aux$vXncol) {
#      object.new$call$vxreg <- NULL
#    } else {
#      object.new$call$vxreg <- object$aux$call.gum$vxreg
#    } ##if there is non-arch variables in vXreg, user needs to specify newvxreg
#
#  } else {

    object.new$call$arch <- NULL
    object.new$call$asym <- NULL
    object.new$call$log.ewma <- NULL
    object.new$call$vxreg <- NULL

#  }

  ##-----------------------------------
  ## pass on arguments to predict.arx:
  ##-----------------------------------

  out <- predict.arx(object.new, spec=spec, n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=newvxreg, newindex=newindex,
    n.sim=n.sim, innov=innov, return=return, plot=plot,
    plot.options=plot.options)

  ##-------------------
  ## return forecasts:
  ##-------------------

  if(return){ return(out) }

} #close predict.gets

##==================================================
## print gets results
print.gets <- function(x, ...)
{
  ##determine spec:
  specType <- switch(as.character(x$call)[1],
    getsm="mean", getsv="variance")

  ##header - first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(specType=="mean"){
    cat("Dependent var.:", x$aux$y.name, "\n")
  }
  cat("Method: Ordinary Least Squares (OLS)\n")

  ##header - if mean:
  if(specType=="mean"){
    cat("Variance-Covariance:", switch(x$aux$vcov.type,
      ordinary = "Ordinary", white = "White (1980)",
      "newey-west" = "Newey and West (1987)"), "\n")
    if(!is.null(x$aux$y.n)){
      cat("No. of observations (mean eq.):", x$aux$y.n, "\n") }
  }

  ##header - if variance:
  if(specType=="variance"){
    if(!is.null(x$aux$loge2.n)){
      cat("No. of observations (variance eq.):",
        x$aux$loge2.n, "\n") }
  }

  ##header - sample info:
  if(!is.null(x$residuals)){
    indexTrimmed <- index(na.trim(x$residuals))
    isRegular <- is.regular(x$residuals, strict=TRUE)
    isCyclical <- frequency(x$residuals) > 1
    if(isRegular && isCyclical){
      cycleTrimmed <- cycle(na.trim(x$residuals))
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
  } #end if(!is.null..)

  ##gum:
  if( specType=="mean" && !is.null(x$gum.mean) ){
    cat("\n")
    cat("GUM mean equation:\n")
    cat("\n")
    printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2),
      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }
  if( !is.null(x$gum.variance) ){
    cat("\n")
    cat("GUM log-variance equation:\n")
    cat("\n")
    if(specType=="mean"){
      printCoefmat(x$gum.variance, signif.stars=FALSE)
    }
    if(specType=="variance"){
      printCoefmat(x$gum.variance, dig.tst=0, tst.ind=c(1,2),
      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
    }
  }
  if(!is.null(x$gum.diagnostics)){
    cat("\n")
    cat("Diagnostics:\n")
    cat("\n")
    printCoefmat(x$gum.diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
  }

  ##paths:
  cat("\n")
  cat("Paths searched: \n")
  cat("\n")
  if(is.null(x$paths)){
    print(NULL)
  }else{
    for(i in 1:length(x$paths)){
      cat("path",i,":",x$paths[[i]],"\n")
    }
  } #end if(is.null(x$paths))

  ##terminal models and results:
  if(!is.null(x$terminals)){
    cat("\n")
    cat("Terminal models: \n")
    if(!is.null(x$terminals)){
      cat("\n")
      for(i in 1:length(x$terminals)){
        cat("spec",i,":",x$terminals[[i]],"\n")
      }
    }
  }
  if(!is.null(x$terminals.results)){
    cat("\n")
    printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
      signif.stars=FALSE)
  }
  
  ##specific model:
  if(specType=="mean" && !is.null(x$specific.spec)){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if(!is.null(x$mean.results)){
      printCoefmat(x$mean.results, signif.stars=FALSE)
    }
    if(x$specific.spec[1]==0){
      cat("empty\n")
    }
##in the future: use estimate.specific=FALSE more directly?
    if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
      cat("Not estimated, since estimate.specific=FALSE\n")
    }
  }
  if(!is.null(x$variance.results)){
    cat("\n")
    cat("SPECIFIC log-variance equation:\n")
    cat("\n")
    printCoefmat(x$variance.results, signif.stars=FALSE)
#    printCoefmat(x$variance.results, dig.tst=0, tst.ind=c(1,2),
#      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }

  ##diagnostics and fit:
  if(!is.null(x$specific.diagnostics)){

    #fit-measures:
    mGOFnames <- "SE of regression"
    mGOF <- sigma.gets(x) #OLD: sqrt( RSS/(nobs-DFs) )
    if(specType == "mean"){
      mGOFnames <- c(mGOFnames, "R-squared")
      mGOF <- rbind(mGOF, rsquared(x))
    }
    mGOFnames <- c(mGOFnames,
      paste("Log-lik.(n=", length(na.trim(x$std.residuals)), ")", sep=""))
    mGOF <- rbind(mGOF, as.numeric(logLik.arx(x)))
    rownames(mGOF) <- mGOFnames
    colnames(mGOF) <- ""

#OLD:
#    #fit-measures:
#    mGOF <- matrix(NA, 3, 1)
#    rownames(mGOF) <- c("SE of regression", "R-squared",
#      paste("Log-lik.(n=", length(na.trim(x$std.residuals)), ")", sep=""))
#    colnames(mGOF) <- ""
#    mGOF[1,1] <- sigma.gets(x) #OLD: sqrt( RSS/(nobs-DFs) )
#    mGOF[2,1] <- rsquared(x) #OLD: x$specific.diagnostics[4,1]
#    mGOF[3,1] <- as.numeric(logLik.arx(x))

    cat("\n")
    cat("Diagnostics and fit:\n")
    cat("\n")
    printCoefmat(x$specific.diagnostics, dig.tst=0, tst.ind=2,
      signif.stars=FALSE)
    printCoefmat(mGOF, digits=6, signif.stars=FALSE)

  }

  ##messages:
  if(!is.null(x$messages)){
    message("\n", appendLF=FALSE)
    message(x$messages)
  }

} #close print.gets

##==================================================
## extract residuals of specific model
residuals.gets <- function(object, std=NULL, ...)
{
  residuals.arx(object, std=std)
} #end residuals.gets

##==================================================
## SE of regression
sigma.gets <- function(object, ...)
{
  sigma.arx(object)
} #close sigma.gets

##==================================================
## summarise output
summary.gets <- function(object, ...)
{
  object.name <- deparse(substitute(object))
  xnames <- cbind(names(object))
  xrows <- nrow(xnames)
  colnames(xnames) <- ""
  rownames(xnames) <- rep(" ",xrows)
  cat("\n")
  cat("Items in list '", object.name, "':", sep="")
  cat("\n")
  print(xnames, quote=FALSE)
  cat("\n")
} #end summary.gets

##==================================================
## extract terminal models
terminals <- function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$terminals)
  }else{
    cat("The object does not belong to the 'gets' nor 'isat' class\n")
  }
} #end terminals

##==================================================
vcov.gets <- function(object, spec=NULL,  ...)
{
  vcov.arx(object, spec=spec)
} #end vcov.gets


####################################################
## 6 ADDITIONAL CONVENIENCE FUNCTIONS
####################################################

##==================================================
##make periodicity (e.g. seasonal) dummies for
##regular time series
periodicdummies <- function(x, values=1)
{
  ##prepare:
  if(!is.regular(x, strict=TRUE)) stop("Vector/matrix not strictly regular")
  iFreq <- frequency(x)
  if(iFreq==1) stop("Frequency must be greater than 1")
  if(!is.zoo(x)){ x <- as.zooreg(x) }
  vCycle <- as.numeric(cycle(x))

  ##values argument:
  if(length(values)==1){ values <- rep(1,iFreq) }
  if(length(values)!=iFreq) stop("length(values) must be 1 or equal to frequency")

  ##make dummies:
  mDums <- matrix(0,NROW(x),iFreq)
  colnames(mDums) <- paste("dum", 1:iFreq, sep="")
  for(i in 1:NCOL(mDums)){
    whereIs <- which(vCycle==i)
    mDums[whereIs,i] <- values[i]
  }

  ##out:
  mDums <- zoo(mDums, order.by=index(x), frequency=iFreq)
  return(mDums)
} #close periodicdummies


##==================================================
## export to EViews
eviews <- function(object, file=NULL, print=TRUE,
  return=FALSE)
{
  out <- list()
  out$object.name <- deparse(substitute(object))

  ##index, data, names:
  out$index <- object$aux$y.index
  out$data <- cbind(object$aux$y, object$aux$mX)
  out$data <- as.data.frame(out$data)
  out$data <- cbind(as.character(out$index), out$data)
  out$names <- c("index", object$aux$y.name, object$aux$mXnames)
  where.mconst <- which(out$names=="mconst")
  if(length(where.mconst) > 0){ out$names[where.mconst] <- "c" }
  colnames(out$data) <- out$names

  ##equation command:
  tmp <- paste(out$names[-1], collapse=" ")
  vcov.type <- NULL
  if(object$aux$vcov.type=="white"){
    vcov.type <- "(cov=white)"
  }
  if(object$aux$vcov.type=="newey-west"){
    vcov.type <- "(cov=hac)"
  }
  out$equation <- paste("equation ", out$object.name,
    ".ls", vcov.type, " ", tmp, sep="")

  ##if print=TRUE and is.null(file):
  if(print && is.null(file)){

    ##EViews code to estimate the model:
    message("EViews code to estimate the model:\n")
    message("  ", out$equation, "\n")

    ##R code to export the data:
    message("R code (example) to export the data of the model:\n")
    message(paste("  eviews(", out$object.name, ", file='C:/Users/myname/Documents/getsdata.csv')\n", sep=""))

  } #close if(print)

  ##if save data:
  if(!is.null(file)){
    write.csv(out$data, file, row.names=FALSE)
    ##if print=TRUE:
    if(print){
      message("Data saved in:\n")
      message("  ", file, "\n", sep="")
      message("EViews code to estimate the model:\n")
      message(" ", out$equation, "\n")
    }
  } #end if(!is.null(file))

  ##out:
  if(return){ return(out) }

} #close eviews

###==================================================
### export to Stata
stata <- function(object, file=NULL, print=TRUE,
  return=FALSE)
{
  out <- list()
  out$object.name <- deparse(substitute(object))

  ##index, data, names:
  out$index <- object$aux$y.index
  out$data <- cbind(object$aux$y, object$aux$mX)
  out$data <- as.data.frame(out$data)
  out$data <- cbind(as.character(out$index), out$data)
  out$names <- gsub("[.]","",tolower(c("index", object$aux$y.name, object$aux$mXnames)))
  where.mconst <- which(out$names=="mconst")
  if(length(where.mconst) > 0){
    out$data <- out$data[-where.mconst]
    out$names <- out$names[-where.mconst]
    noConstant <- FALSE
  }else{
    noConstant <- TRUE
  }
  colnames(out$data) <- out$names

  ##Stata code to estimate the model:
  outNames <- out$names
  outNames[1] <- "regress"
  out$regress <- paste(outNames, collapse=" ")
  if( noConstant==TRUE || object$aux$vcov.type!="ordinary" ){

    cmdOptions <- NULL
    if(noConstant){ cmdOptions <- c(cmdOptions, "noconstant") }
    if(object$aux$vcov.type!="ordinary"){
      cmdOptions <- c(cmdOptions, "vce(robust)")
    }
    cmdOptions <- paste(cmdOptions, collapse=" ")
    out$regress <- paste(out$regress, ",", cmdOptions, collapse="")

  }

  ##if print=TRUE and is.null(file):
  if(print && is.null(file)){

    ##Stata code to estimate the model:
    message("Stata code to estimate the model:\n")
    message(" ", out$regress, "\n")

    ##R code to export the data:
    message("R code (example) to export the data of the model:\n")
    message(paste("  stata(", out$object.name, ", file='C:/Users/myname/Documents/getsdata.csv')\n", sep=""))

  } #close if(print && is.null(file))

  ##if save data:
  if(!is.null(file)){
    write.csv(out$data, file, row.names=FALSE)
    ##if print=TRUE:
    if(print){
      message("Data saved in:\n")
      message("  ", file, "\n", sep="")
      message("Stata code to estimate the model:\n")
      message(" ", out$regress, "\n")
    }
  } #end if(!is.null(file))

  ##out:
  if(return){ return(out) }

} #close stata()

##==================================================
##generate latex-code (equation form):
printtex <- function(x, fitted.name=NULL, xreg.names=NULL,
  digits=3, intercept=TRUE, gof=TRUE, diagnostics=TRUE)
{

  ##is class(x)="arx"/"gets"/"isat"?:
  ##---------------------------------
  
  xName <- deparse(substitute(x))
  xClass <- class(x)
  if( xClass %in% c("arx","gets","isat") ){
    yName <- ifelse(is.null(fitted.name), x$aux$y.name, fitted.name)
  }else{
    yName <- ifelse(is.null(fitted.name), "y", fitted.name)
    message(paste0("\n '", xName, "'", " is not of class 'arx', ",
      "'gets' or 'isat', LaTeX code may contain errors:\n"))
  }
  yName <- paste0("\\widehat{", yName, "}")

  ##equation:
  ##---------

  coefs <- coef(x)
  if(is.null(xreg.names)){
    coefsNames <- names(coefs)
  }else{
    coefsNames <- xreg.names
    if( length(coefs) != length(xreg.names) ){
      message(paste0("\n length of 'xreg.names' does not match",
        " length of 'coef(x)'\n"))
    }
  }
  intercept <- as.numeric(intercept)
  if( intercept > 0 ){ coefsNames[ intercept ] <- "" }
  coefs <- as.numeric(coefs)
  stderrs <- as.numeric(sqrt(diag(vcov(x))))

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

  txtAddEq <- ifelse(gof+diagnostics>0, " \\\\[2mm]", "")
  eqtxt <- paste0("  ", yName, " &=& ", eqtxt, "", txtAddEq, " \n")

  ##goodness of fit:
  ##----------------

  goftxt <- NULL
  if(gof){
    goftxt <- "   &&"
    iT <- ""
    if(xClass %in% c("arx","gets","isat") ){
      goftxt <- paste(goftxt, " R^2=",
        format(round(rsquared(x), digits=digits), nsmall=digits),
        " \\qquad \\widehat{\\sigma}=",
        format(round(sigma(x), digits=digits), nsmall=digits),
        sep="")
      iT <- x$aux$y.n
    }
    goftxt <- paste(goftxt, " \\qquad LogL=",
      format(round(as.numeric(logLik(x)), digits=digits), nsmall=digits),
        "\\qquad T = ", iT, " \\nonumber \\\\ \n", sep="")
  }
  
  ##diagnostics:
  ##------------

  diagtxt <- NULL
  if(xClass=="arx" && diagnostics==TRUE){
    dfDiags <- diagnostics(x,
      ar.LjungB=c(ar.LjungB=x$aux$qstat.options[1],1),
      arch.LjungB=c(ar.LjungB=x$aux$qstat.options[2],1),
      normality.JarqueB=TRUE, verbose=TRUE)
    diagtxt <- paste("  ", " && \\underset{[p-val]}{ AR(",
      x$aux$qstat.options[1], ") }:", " \\underset{[",
      format(round(dfDiags[1,3], digits=digits), nsmall=digits), "]}{",
      format(round(dfDiags[1,1], digits=digits), nsmall=digits), "}",
      "\\qquad \\underset{[p-val]}{ ARCH(",
      x$aux$qstat.options[2], ")}:", "\\underset{[",
      format(round(dfDiags[2,3], digits=digits), nsmall=digits), "]}{",
      format(round(dfDiags[2,1], digits=digits), nsmall=digits), "}",
      "\\qquad \\underset{[p-val]}{ Normality }:", "\\underset{[",
      format(round(dfDiags[3,3], digits=digits), nsmall=digits), "]}{",
      format(round(dfDiags[3,1], digits=digits), nsmall=digits), "}",
      " \\nonumber \n", sep="")
  }
  
  ##print code:
  ##-----------

  cat("\\begin{eqnarray}\n")
  cat(eqtxt)
  cat(goftxt)
  cat(diagtxt)
  cat("\\end{eqnarray}\n")

}   #close printtex

