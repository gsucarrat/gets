####################################################
## This file contains the base-source of the gets
## package.
##
## Current version: 0.25
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
## create S3 generics/methods:
## - gets
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
## regressorsMean
## regressorsVariance
## getsFun
## blocksFun
##
####################################################
## 3 ARX FUNCTIONS
####################################################
##
## arx
## coef.arx         #extraction functions
## ES               #(some are S3 methods)
## fitted.arx
## gets.arx
## logLik.arx
## plot.arx
## predict.arx
## print.arx
## recursive
## residuals.arx
## sigma.arx
## rsquared
## summary.arx
## toLatex.arx
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
## toLatex.gets
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
## as.lm      #convert arx/gets/isat object to 'lm' object
##
####################################################


####################################################
## 1 INITIATE
####################################################
                                                       
##==================================================
##create S3 generic/method 'gets':
gets <- function(x, ...){ UseMethod("gets") }
###test the method:
#set.seed(123); y <- arima.sim(list(ar=0.4, ma=0.1), 100)
#mod01 <- arima(y, order=c(1,0,1))
#gets.Arima <- function(x){ print("It works - cool!") }
#gets(mod01)

##for the future?:
#gets.default <- function(x, ...){ print("It works!") }

##for the future?:
##==================================================
##create S3 method 'isat':
#isat <- function(x, ...){ UseMethod("isat") }
###test the method:
#set.seed(123); y <- arima.sim(list(ar=0.4, ma=0.1), 100)
#mod01 <- arima(y, order=c(1,0,1))
#isat.Arima <- function(x){ print("Cool!") }
#isat(mod01)

##for the future?:
#isat.default <- today's isat function


####################################################
##2 BASE FUNCTIONS
####################################################

##==================================================
##diagnostics checking
diagnostics <- function(x, ar.LjungB=c(1,0.025), arch.LjungB=c(1,0.025),
  normality.JarqueB=NULL, verbose=TRUE, user.fun=NULL, ...)
{
  ##initiate:
  ##---------
  if( is.null(x$std.residuals) ){
    zhat <- x$residuals
  }else{
    zhat <- x$std.residuals
  }
  diagnosticsGood <- TRUE

  ##test for autocorrelation:
  ##-------------------------
  if( !is.null(ar.LjungB) ){
    ar.LjungBox <- Box.test(zhat, lag=ar.LjungB[1], type="L")
    if( ar.LjungBox$p.value <= ar.LjungB[2] ){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  }

  ##test for arch:
  ##--------------
  if( diagnosticsGood && !is.null(arch.LjungB) ){
    zhat2 <- zhat^2
    arch.LjungBox <- Box.test(zhat2, lag=arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){
      diagnosticsGood <- FALSE
      diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  }


  ##test for normality:
  ##-------------------
  if( diagnosticsGood && !is.null(normality.JarqueB) ){
    zhatadj <- coredata(na.trim(zhat))
    n <- length(zhatadj)
    avgzhat <- mean(zhatadj) #do I really need this?
    zhat.avgzhat <- zhatadj-avgzhat #do I really need this?
    zhat.avgzhat2 <- zhat.avgzhat^2
    K <- n*sum(zhat.avgzhat^4)/(sum(zhat.avgzhat2)^2)
    S <- (sum(zhat.avgzhat^3)/n)/(sum(zhat.avgzhat2)/n)^(3/2)
    JB <- (n/6)*(S^2 + 0.25*((K-3)^2))
    JBpval <- pchisq(JB, df = 2, lower.tail=FALSE)
    if(JBpval <= normality.JarqueB){
        diagnosticsGood <- FALSE
        diagnosticsGood <- as.logical(max(diagnosticsGood,verbose))
    }
  }

  ##user-defined test(s):
  ##---------------------
  if( diagnosticsGood && !is.null(user.fun) ){
    ##make user.fun argument:
    userFunArg <- user.fun
    userFunArg$name <- NULL
    userFunArg$envir <- NULL
    userFunArg$pval <- NULL
    if( length(userFunArg)==0 ){ userFunArg <- NULL }
    ##'do' user diagnostics:
    if( is.null(user.fun$envir) ){
      userVals <- do.call(user.fun$name, c(list(x=x),userFunArg))
    }else{
      userVals <- do.call(user.fun$name, c(list(x=x),userFunArg),
        envir=user.fun$envir)
    }
    userVals <- rbind(userVals)
    if( !is.null(user.fun$pval) ){
      userFunPval <- as.numeric(userVals[,3])
      if( any(userFunPval <= user.fun$pval) ){
        diagnosticsGood <- FALSE
      }
    }
  } #end if( user.fun )

  ##result:
  ##-------

  ##if(!verbose): return logical only
  if( !verbose ){ result <- diagnosticsGood }

  ##if(verbose): return diagnostics table
  if( verbose ){
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
eqwma <- function(x, length=5, k=1, p=1, abs=FALSE, log=FALSE,
  as.vector=FALSE, lag=NULL, start=NULL)
{
  ##deprecated arguments:
  if(!is.null(start)){ stop("'start' has been deprecated") }
  if(!is.null(lag)){ stop("'lag' has been deprecated, use 'k' instead") }

  ##zoo related:
  if(is.zoo(x)){
    isZoo <- TRUE
    xindex <- index(x)
    x <- coredata(x)
  }else{
    isZoo <- FALSE
  }

  ##is x a matrix?
  if(NCOL(x)>1){
    stop("'x' must be one-dimensional")
  }else{
    x <- as.vector(x)
  }

  ##na's, power p, absolute value:
  xn <- length(x)
  x <- na.trim(x, sides="left")
  xnadj <- length(x)
  if(xn > xnadj){nachk <- TRUE}else{nachk <- FALSE}
  if(abs){xabs <- abs(x)}else{xabs <- x}
  if(p!=1){xabsp <- xabs^p}else{xabsp <- xabs}

  ##compute eqwma:
  result <- NULL
  for(i in 1:length(length)){
    movavg <- rollmean(xabsp, length[i], fill=NA, align="r")
    result <- cbind(result,movavg)
  }

  ##lag?:
  if(k>0){
    result <- cbind(result[ 1:c(NROW(result)-k), ]) #lag
    result <- rbind( matrix(NA,k,NCOL(result)), result) #add NAs
  }
  
  #apply log?:
  if(log){
    result <- log(result)
    colnames(result) <- paste0("logEqWMA(", length, ")")
  }else{
    colnames(result) <- paste0("EqWMA(", length, ")")
  }
  if(nachk){ result <- rbind( matrix(NA,c(xn-xnadj),NCOL(result)), result) }
  if(as.vector && NCOL(result)==1 ){ result <- as.vector(result) }
  if(isZoo){ result <- zoo(result, order.by=xindex) }

  return(result)
} #close eqwma

##==================================================
##Compute information criterion:
infocrit <- function(x, method=c("sc", "aic", "aicc", "hq"))
{
  method <- match.arg(method)
  ##Schwarch criterion:
  if(method == "sc") infoVal <- -2*x$logl/x$n + x$k*log(x$n)/x$n
  ##Akaike criterion:
  if(method == "aic") infoVal <- -2*x$logl/x$n + 2*x$k/x$n
  ##AICC of Hurvich and Tsai (1989):
  if(method == "aicc") infoVal <- -2*x$logl/x$n + 2*x$k/(x$n-x$k-1)
  ##Hannan-Quinn criterion:
  if(method == "hq") infoVal <- -2*x$logl/x$n + 2*x$k*log(log(x$n))/x$n
  ##return:
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
##make log(EqWMA) terms
leqwma <- function(x, length=5, k=1, p=2, as.vector=FALSE,
  lag=NULL, start=NULL)
{
  eqwma(x, length=length, k=k, p=p, log=TRUE, abs=TRUE,
    as.vector=as.vector, lag=lag, start=start)
} #close leqwma

##==================================================
##OLS estimation using the QR decomposition
ols <- function(y, x, untransformed.residuals=NULL, tol=1e-07,
  LAPACK=FALSE, method=3, variance.spec=NULL, ...)
{

  ##for the future:
  ## - new argument: options=NULL (default), to control how the
  ## Newey and West (1987) coefficient-covariance is computed,
  ## amongst other
  ## - rename ols to estFun? Split estFun into two functions,
  ## estFun and vcovFun?

  ##user-specified:
  ##---------------
  if(method==0){
    stop("method = 0 has been deprecated")
  }

  ##fastest, usually only for estimates:
  ##------------------------------------
  if(method==1){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK)
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
  }

  ##second fastest (slightly more output):
  ##--------------------------------------
  if(method==2){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$residuals <- y - out$fit
  }

  ##ordinary vcov:
  ##--------------
  if(method==3){

    ##mean specification:
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr) #(x'x)^-1
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

  } #close method=3

  ##White (1980) vcov:
  ##------------------
  if(method==4){

    ##mean specification:
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr) #(x'x)^-1
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

    ##variance specification:
    
  }

  ##Newey and West(1987) vcov:
  ##--------------------------
  if(method==5){

    ##mean specification:
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k>0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr) #(x'x)^-1
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

    ##variance specification:
    
  }

  ##log-variance w/ordinary vcov (note: y = log(e^2)):
  ##--------------------------------------------------
  if(method==6){

    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr) #(x'x)^-1
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

  ##if variance specification:
  ##--------------------------
  if( !is.null(variance.spec) ){

    if(method==6){ stop("not compatible with method=6") }
    if( !is.null(variance.spec$vxreg) ){
      if( length(y)!=NROW(variance.spec$vxreg) ){
        stop("length(y) != NROW(vxreg)")
      }
      variance.spec$vxreg <- coredata(variance.spec$vxreg)
    }
    e <- out$residuals
    variance.spec <- c(list(e=e), variance.spec)
    variance.spec$return.regressand <- TRUE #some protection
    variance.spec$return.as.zoo <- FALSE
    variance.spec$na.trim <- TRUE #some protection
    variance.spec$na.omit <- FALSE #--||--
    tmp <- do.call("regressorsVariance", variance.spec)
    loge2 <- tmp[,1]
    vX <- cbind(tmp[,-1])
    e <- e[c(length(e)-length(loge2)+1):length(e)]
    estVar <- ols(loge2, vX, untransformed.residuals=e, tol=tol,
      LAPACK=LAPACK, method=6)
    out$regressorsVariance <- tmp
    out$var.coefficients <- estVar$coefficients
    out$Elnz2 <- estVar$Elnz2
    out$vcov.var <- estVar$vcov
    NAs2add <- rep(NA, length(y)-length(loge2))
    out$var.fit <- c(NAs2add, estVar$var.fit)
    out$std.residuals <- c(NAs2add, estVar$std.residuals)
    out$ustar.residuals <- c(NAs2add, estVar$residuals)
    out$logl <- estVar$logl

  }

  ##return result:
  ##--------------
  return(out)

} #close ols() function

##==================================================
##Create the mean regressors of an arx model:
regressorsMean <- function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  return.regressand=TRUE, return.as.zoo=TRUE, na.trim=TRUE,
  na.omit=FALSE)
{

  ##regressand:
  y.name <- deparse(substitute(y))
  if(is.zoo(y)){ y <- cbind(y) }else{ y <- as.zoo(cbind(y)) }
  if(NCOL(y) > 1) stop("Dependent variable not 1-dimensional")
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y.n <- NROW(y)
  y.index <- index(y)
  y <- coredata(y)
  t1 <- y.index[1]
  t2 <- y.index[y.n]

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
    ewma$as.vector <- FALSE #force result to be a matrix
    tmp <- do.call(eqwma, c(list(y),ewma) )
    mXnames <- c(mXnames, colnames(tmp))
    colnames(tmp) <- NULL
    mX <- cbind(mX, tmp)
  }

  ##trim for NAs:
  if(na.trim){
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    y.n <- NROW(tmp) #re-define y.n
    y.index <- index(tmp) #re-define y.index
    t1 <- y.index[1] #re-define t1
    t2 <- y.index[y.n] #re-define t2
    y <- coredata(tmp[,1]) #re-define y
    if(!is.null(mX)){ #re-define mX
      mX <- tmp[,2:NCOL(tmp)]
      mX <- coredata(mX)
      mX <- cbind(mX)
      colnames(mX) <- NULL
    }
  }
  
  ##mxreg:
  if( !is.null(mxreg) ){
#OLD (until version 0.23):
#  if( !is.null(mxreg) && !identical(as.numeric(mxreg),0) ){
    mxreg <- as.zoo(cbind(mxreg))
    mxreg.names <- colnames(mxreg)
    if(is.null(mxreg.names)){
      mxreg.names <- paste("mxreg", 1:NCOL(mxreg), sep="")
    }
    if(any(mxreg.names == "")){
      missing.colnames <- which(mxreg.names == "")
      for(i in 1:length(missing.colnames)){
        mxreg.names[missing.colnames[i]] <- paste0("mxreg", i)
      }
    }
    #mxreg.names <- make.names(mxreg.names)
    mXnames <- c(mXnames, mxreg.names)
    mxreg <- window(mxreg, start=t1, end=t2)
    mxreg <- cbind(coredata(mxreg))
    mX <- cbind(mX, mxreg)

    ##re-trim for NAs:
    if(na.trim){
      tmp <- zoo(cbind(y,mX), order.by=y.index)
      tmp <- na.trim(tmp, sides="both", is.na="any")
      y.n <- NROW(tmp) #re-define y.n
      y.index <- index(tmp) #re-define y.index
      t1 <- y.index[1] #re-define t1
      t2 <- y.index[y.n] #re-define t2
      y <- coredata(tmp[,1])
      mX <- tmp[,2:NCOL(tmp)]
      mX <- coredata(mX)
      mX <- cbind(mX)
      colnames(mX) <- NULL
    }

  } #end if(!is.null(mxreg))

  ##remove rows with NAs:
  if(na.omit){
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.omit(tmp)
    y.n <- NROW(tmp) #re-define y.n
    y.index <- index(tmp) #re-define y.index
    t1 <- y.index[1] #re-define t1
    t2 <- y.index[y.n] #re-define t2
    y <- coredata(tmp[,1]) #re-define y
    if(!is.null(mX)){ #re-define mX
      mX <- tmp[,2:NCOL(tmp)]
      mX <- coredata(mX)
      mX <- cbind(mX)
      colnames(mX) <- NULL
    }
  }

  ### OUTPUT: ######################

  if(return.regressand){
    result <- cbind(y, mX)
    colnames(result) <- c(y.name, mXnames)
  }else{
    result <- mX
    if(!is.null(result)){ colnames(result) <- mXnames }
  }
  if(return.as.zoo && !is.null(result) ){ result <- zoo(result, order.by=y.index) }
  return(result)
  
} #close regressorsMean()

##==================================================
##Create the variance regressors of an arch-x model:
regressorsVariance <- function(e, vc=TRUE, arch=NULL, asym=NULL,
  log.ewma=NULL, vxreg=NULL, zero.adj=0.1, vc.adj=TRUE,
  return.regressand=TRUE, return.as.zoo=TRUE, na.trim=TRUE,
  na.omit=FALSE)
{

  ##regressand:
  if(is.zoo(e)){ e <- cbind(e) }else{ e <- as.zoo(cbind(e)) }
  if(NCOL(e) > 1) stop("Dependent variable not 1-dimensional")
  e.n <- NROW(e)
  loge2.index <- index(e)
  e <- coredata(e)
  t1 <- loge2.index[1]
  t2 <- loge2.index[e.n]
  zero.where <- which(e==0)
  eabs <- abs(e)
  if( length(zero.where)>0 ){
    eabs[zero.where] <- quantile(eabs[-zero.where], zero.adj, na.rm=TRUE)
  }
  loge2 <- log(eabs^2)

  ##create regressor matrix:
  vX <- NULL
  vXnames <- NULL

  ##variance intercept:
  if( identical(as.numeric(vc),1) ){
    vX <- cbind(rep(1,e.n))
    vXnames <- "vconst"
  }

  ##arch terms:
  if(!is.null(arch) && !identical(as.numeric(arch),0) ){
    tmp <- NULL
    nas <- rep(NA, max(arch))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],loge2[1:c(e.n-i)]))
    }
    tmpfun <- sapply(arch,tmpfun)
    vX <- cbind(vX, tmp)
    vXnames <- c(vXnames, paste0("arch", arch))
  }

  ##asym terms:
  if(!is.null(asym) && !identical(as.numeric(asym),0) ){
    tmp <- NULL
    nas <- rep(NA, max(asym))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],
        loge2[1:c(e.n-i)]*as.numeric(e[1:c(e.n-i)]<0)))
    }
    tmpfun <- sapply(asym,tmpfun)
    vX <- cbind(vX, tmp)
    vXnames <- c(vXnames, paste0("asym", asym))
  }

  ##log.ewma term:
  if(!is.null(log.ewma)){
    if(is.list(log.ewma)){
      log.ewma$k <- 1
    }else{
      log.ewma <- list(length=log.ewma)
    }
    tmp <- do.call(leqwma, c(list(e),log.ewma) )
    vXnames <- c(vXnames, colnames(tmp))
    colnames(tmp) <- NULL
    vX <- cbind(vX, tmp)
  }

  ##trim for NAs:
  if(na.trim){
    tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    loge2.n <- NROW(tmp)
    loge2.index <- index(tmp) #re-define index
    t1 <- loge2.index[1] #re-define t1
    t2 <- loge2.index[loge2.n] #re-define t2
    loge2 <- tmp[,1]
    loge2 <- coredata(loge2)
    if(!is.null(vX)){ #re-define vX
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }
  }

  ##vxreg:
  if(!is.null(vxreg)){
    vxreg <- as.zoo(cbind(vxreg))
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste0("vxreg", 1:NCOL(vxreg))
    }
    if(any(vxreg.names == "")){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[missing.colnames[i]] <- paste0("vxreg", i)
      }
    }
    vXnames <- c(vXnames, vxreg.names)
    vxreg <- window(vxreg, start=t1, end=t2)
    vxreg <- cbind(coredata(vxreg))
    vX <- cbind(vX, vxreg)
    colnames(vxreg) <- NULL

    ##re-trim for NAs:
    if(na.trim){
      tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
      tmp <- na.trim(tmp, sides="both", is.na="any")
      loge2.n <- NROW(tmp)
      loge2.index <- index(tmp) #re-define index
      t1 <- loge2.index[1] #re-define t1
      t2 <- loge2.index[loge2.n] #re-define t2
      loge2 <- coredata(tmp[,1])
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }

  } #end if(!is.null(vxreg))

  ##remove rows with NAs:
  if(na.omit){
    tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
    tmp <- na.omit(tmp)
    loge2.n <- NROW(tmp) #re-define
    loge2.index <- index(tmp) #re-define
    t1 <- loge2.index[1] #re-define t1
    t2 <- loge2.index[loge2.n] #re-define t2
    loge2 <- coredata(tmp[,1]) #re-define
    if(!is.null(vX)){ #re-define vX
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }
  }

  ### OUTPUT: ######################

  if(return.regressand){
    result <- cbind(loge2, vX)
    colnames(result) <- c("loge2", vXnames)
  }else{
    result <- vX
    if(!is.null(result)){ colnames(result) <- vXnames }
  }
  if(return.as.zoo && !is.null(result) ){ result <- zoo(result, order.by=loge2.index) }
  return(result)

} #close regressorsVariance()

##==================================================
##do gets fast and with full flexibility (for advanced users)
getsFun <- function(y, x, untransformed.residuals=NULL,
  user.estimator=list(name="ols"), gum.result=NULL, t.pval=0.05,
  wald.pval=t.pval, do.pet=TRUE, ar.LjungB=NULL, arch.LjungB=NULL,
  normality.JarqueB=NULL, user.diagnostics=NULL,
  gof.function=list(name="infocrit"), gof.method=c("min","max"),
  keep=NULL, include.gum=FALSE, include.1cut=FALSE,
  include.empty=FALSE, max.paths=NULL, turbo=FALSE, tol=1e-07,
  LAPACK=FALSE, max.regs=NULL, print.searchinfo=TRUE, alarm=FALSE)
{
  ## DO NOT:
  ## - introduce a check of the type NROW(y)==NCOL(x), since this will
  ##   invalidate situations where the x's contain coefficients rather
  ##   than regressors (e.g. when models are non-linear in parameters)
  ## TO DO:
  ## - introduce check for is.vector(y)==TRUE?
  ## - introduce check for is.matrix(x)==TRUE?
  ## - let out$specific.spec be equal to the GUM in the case where
  ##   all regressors are significant in the GUM?
  ## - if gof.function is not default, e.g. adjusted R-squared, then
  ##   it seems the value of logl is added to terminals.results
  ##   unnecessarily. Look into?
  ## - turbo: replace length(regsDeleteList) with regsDeleteList.n?
  ## - turbo: redefine regsFun function (careful!: setequal is delicate)
  ## - envir argument in user.estimator: change default behaviour?

  ## contents:
  ## 1 arguments
  ## 2 initialise
  ## 3 gum
  ## 4 1-cut model
  ## 5 empty model
  ## 6 multi-path search
  ## 7 find the best model
  ## 8 output


  ##-----------------------
  ## 1 arguments
  ##-----------------------

  gof.method <- match.arg(gof.method)

  ##y, x, make auxiliary list:
  if( is.null(x) || NCOL(x)==0 ){ stop("GUM regressor matrix is empty") }
  x <- cbind(x) #ensure x is a matrix
  aux <- list()
  aux$y.n <- NROW(y)
  aux$xNCOL <- NCOL(x)

  ##make user-estimator argument:
  userEstArg <- user.estimator
  userEstArg$name <- NULL
  userEstArg$envir <- NULL
  if( length(userEstArg)==0 ){ userEstArg <- NULL }

  ##make gof.function argument:
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


  ##-----------------------
  ## 2 initialise
  ##-----------------------

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

  ##if all regressors in keep, add gum to terminals:
  if( delete.n==0 && include.gum==FALSE ){
    include.gum <- TRUE
    out$messages <- paste(out$messages,
      "- All regressors in 'keep', GUM added to terminals", sep="")
  }

  ##-----------------------
  ## 3 gum
  ##-----------------------

  ##estimate GUM:
  if( is.null(gum.result) ){
    if( is.null(user.estimator$envir) ){
      est <- do.call(user.estimator$name, c(list(y,x), userEstArg))
    }else{
      est <- do.call(user.estimator$name,
        c(list(y,x), userEstArg), envir=user.estimator$envir)
    }
    out$no.of.estimations <- out$no.of.estimations + 1
  }else{ est <- gum.result }

  ##do diagnostics:
  if( doDiagnostics ){
    gumDiagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
      arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
      verbose=FALSE, user.fun=user.diagnostics)
  }else{ gumDiagnosticsOK <- TRUE }

  ## if GUM passes diagnostic checks:
  if( gumDiagnosticsOK ){

    ##record data for Wald-tests (pets) against gum:
    gum.regs <- gum
    gum.coefs <- est$coefficients
    gum.varcovmat <- est$vcov

    ##compute stderrs, t-stats, p-vals:
    stderrs <- sqrt(diag(est$vcov))
    gum.tstat <- est$coefficients/stderrs
    gum.pval <- pt(abs(gum.tstat), est$df, lower.tail=FALSE)*2

    ##these two lines are repeated later under 1-cut, and - if
    ##max.paths < n.paths - adjusted in the multi-path search:
    insig.regs <- setdiff( which(gum.pval > t.pval), keep)
    n.paths <- length(insig.regs)

    ##if all regressors significant, ensure gum is a terminal:
    if( n.paths==0 ){ include.gum <- TRUE }

    if( include.gum ){

      out$terminals[[1]]  <- gum #add gum to list of terminal specs

      ##specification results
      if( is.null(gof.function$envir) ){
        gofValue <- do.call(gof.function$name, c(list(est),gofFunArg))
      }else{
        gofValue <- do.call(gof.function$name,
          c(list(est),gofFunArg), envir=gof.function$envir)
      }
      out$terminals.results <- rbind(out$terminals.results,
        c(gofValue, est$logl, est$n, est$k))
      row.labels <- c(row.labels, "spec 1 (gum):")

    } #end if(include.gum)

  }else{
    out$messages <- paste(out$messages,
      "- GUM does not pass one or more diagnostic checks", sep="")
  }


  ##-----------------------
  ## 4 1-cut model
  ##-----------------------

  if( gumDiagnosticsOK && delete.n>0 && include.1cut ){

    ##all non-keep regressors significant:
    if( n.paths==0 ){
      out$messages <- paste(out$messages,
        "- 1-CUT not included (all non-keep regressors are significant)",
        sep="")
    }

    ##one or more non-keep regressor insignificant:
    if( n.paths>0 ){

      ##estimate 1cut:
      mXadj <- cbind(x[,-insig.regs])
      if( is.null(user.estimator$envir) ){
        est <- do.call(user.estimator$name, c(list(y,mXadj), userEstArg))
      }else{
        est <- do.call(user.estimator$name, c(list(y,mXadj), userEstArg),
          envir=user.estimator$envir)
      }
      out$no.of.estimations <- out$no.of.estimations + 1

      ##do diagnostics:
      if( doDiagnostics ){
        diagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          verbose=FALSE, user.fun=user.diagnostics)
      }else{ diagnosticsOK <- TRUE }

      ## if 1cut passes diagnostic checks:
      if( diagnosticsOK ){

        ## do pet (i.e. wald-test):
        if( do.pet ){
          mR <- rbind(aux$mR[insig.regs,])
          mRestq <- mR %*% cbind(gum.coefs)
          wald.stat <- t(mRestq)%*%qr.solve(mR%*%gum.varcovmat%*%t(mR), tol=tol) %*% mRestq
          petOK <- as.logical(wald.pval < pchisq(wald.stat, n.paths, lower.tail = FALSE))
        }else{ petOK <- TRUE }

        ##add 1-cut to terminals?:
        if( petOK ){

          #add 1cut to list of terminal specs:
          spec.1cut <- setdiff(gum, insig.regs)
          out$terminals[[ length(out$terminals)+1 ]] <- spec.1cut

          ##specification results
          if( is.null(gof.function$envir) ){
            gofValue <- do.call(gof.function$name, c(list(est),
              gofFunArg))
          }else{
            gofValue <- do.call(gof.function$name, c(list(est),
              gofFunArg), envir=gof.function$envir)
          }
          out$terminals.results <- rbind(out$terminals.results,
            c(gofValue, est$logl, est$n, est$k))
          row.labels <- c(row.labels,
            paste("spec ", length(out$terminals), " (1-cut):", sep=""))

        } #end if(petOK)

      } ##end if(diagnosticsOK)

    } ###end if(n.paths > 0)

  } ####end if(1-cut model)


  ##-----------------------
  ## 5 empty model
  ##-----------------------

  if( gumDiagnosticsOK && delete.n>0 && include.empty ){

    ##Here: do NOT do pet in order to enable reality check!

    ##check if empty = 1-cut:
    if( include.1cut && exists("spec.1cut") ){
      emptyEqualTo1cut <- identical(keep, spec.1cut)
    }else{ emptyEqualTo1cut <- FALSE }

    ##empty equal to 1cut?:
    if( emptyEqualTo1cut ){
        out$messages <- paste0(out$messages,
          "- The empty model is equal to the 1-cut model")
    }else{

      ## estimate model:
      mXadj <- cbind(x[,keep])
      if( is.null(user.estimator$envir) ){
        est <- do.call(user.estimator$name, c(list(y,mXadj), userEstArg))
      }else{
        est <- do.call(user.estimator$name, c(list(y,mXadj), userEstArg),
          envir=user.estimator$envir)
      }
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

        ##specification results
        if( is.null(gof.function$envir) ){
          gofValue <- do.call(gof.function$name, c(list(est),
            gofFunArg))
        }else{
          gofValue <- do.call(gof.function$name, c(list(est),
            gofFunArg), envir=gof.function$envir)
        }
        out$terminals.results <- rbind(out$terminals.results,
          c(gofValue, est$logl, est$n, est$k))
        row.labels <- c(row.labels,
          paste("spec ", length(out$terminals), " (empty):", sep=""))

      }else{

          out$messages <- paste0(out$messages,
            "- Empty model not included (it does not pass one or more diagnostics)")

      } #end if(empty passes diagnostics==TRUE){..}else{..}

    } ##end if( emptyEqualTo1cut )else(...)

  } ###end if(include empty model==TRUE)


  ##-----------------------
  ## 6 multi-path search
  ##-----------------------

#OLD (already set under 'gum'):
#insig.regs <- NULL
pathsTerminals <- list()
if( gumDiagnosticsOK && delete.n>0 ){

  ##number of paths:
#OLD:
#  insig.regs <- setdiff( which(gum.pval > t.pval), keep)
  if( !is.null(max.paths) ){
    if(max.paths < length(insig.regs)){
      pvalRanksInv <- rank( 1-gum.pval[insig.regs] )
      insig.regs <- insig.regs[ pvalRanksInv <= max.paths ]
    }
  }
  n.paths <- length(insig.regs) #re-define n.paths

  ## if paths = 0:
  if(n.paths == 0){
    out$messages <- paste(out$messages,
      "- All non-keep regressors significant in GUM", sep="")
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
        if( is.null(user.estimator$envir) ){
          est <- do.call(user.estimator$name, c(list(y,mXadj),
            userEstArg))
        }else{
          est <- do.call(user.estimator$name, c(list(y,mXadj),
            userEstArg), envir=user.estimator$envir)
        }
        out$no.of.estimations <- out$no.of.estimations + 1

        ##do diagnostics:
        if(doDiagnostics){
          diagnosticsOK <- diagnostics(est, ar.LjungB=ar.LjungB,
            arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
            verbose=FALSE, user.fun=user.diagnostics)
        }else{ diagnosticsOK <- TRUE }

        ## move regressor to keep.adj?:
        if( !diagnosticsOK ){
          path.n <- length(path)
          keep.adj <- union(path[path.n], keep.adj)
          path <- union(path, path[path.n]*c(-1))
          next #next j
        }

        ## if empty model passes diagnostic checks:
        if( diagnosticsOK ){

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
        if( is.null(gof.function$envir) ){
          gofValue <- do.call(gof.function$name, c(list(est),
            gofFunArg))
        }else{
          gofValue <- do.call(gof.function$name, c(list(est),
            gofFunArg), envir=gof.function$envir)
        }
        out$terminals.results <- rbind(out$terminals.results,
          c(gofValue, est$logl, est$n, est$k))
        row.labels <- c(row.labels, paste("spec ", length(out$terminals), ":", sep=""))

      } #end if(chk.spec==FALSE)

    } ##end multi-path search: for(i in 1:n.paths) loop

  } ###end if paths > 0

} #####end if( gumDiagnosticsOK && delete.n>0 )


  ##-----------------------
  ## 7 find the best model
  ##-----------------------

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


  ##-----------------------
  ## 8 output
  ##-----------------------

  out$time.finished <- date()
  if(alarm){ alarm() }
  return(out)

} #close getsFun function


##==================================================
##do block-based gets with full flexibility (for advanced users)
blocksFun <- function(y, x, untransformed.residuals=NULL,
  blocks=NULL, no.of.blocks=NULL, max.block.size=30,
  ratio.threshold=0.8, gets.of.union=TRUE, force.invertibility=FALSE,
  user.estimator=list(name="ols"), t.pval=0.001, wald.pval=t.pval,
  do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL,
  normality.JarqueB=NULL, user.diagnostics=NULL,
  gof.function=list(name="infocrit"), gof.method=c("min","max"),
  keep=NULL, include.gum=FALSE, include.1cut=FALSE,
  include.empty=FALSE, max.paths=NULL, turbo=FALSE,
  parallel.options=NULL, tol=1e-07, LAPACK=FALSE, max.regs=NULL,
  print.searchinfo=TRUE, alarm=FALSE)
{
  ## contents:
  ## 1 initiate
  ## 2 x and blocks arguments
  ## 3 loop on x matrices
  ## 4 make return object
    
  ##-------------------------------
  ## 1 initiate
  ##-------------------------------

  gof.method <- match.arg(gof.method)

  ##make result list, add to list:
  result <- list()
  result$call <- sys.call()
  result$time.started <- date()
  result$time.finished <- NA
  result$messages <- NULL

  ##parallel.options argument:
  if(!is.null(parallel.options)){

    ##if(numeric):
    if(is.numeric(parallel.options)){
      clusterSpec <- parallel.options
      OScores <- detectCores()
      if(parallel.options > OScores){
        stop("parallel.options > number of cores/threads")
      }
    }

    ##varlist for clusterExport:
    if(is.list(parallel.options)){
      clusterVarlist <- parallel.options$varlist
    }else{
      clusterVarlist <- NULL
    }
    clusterVarlist <- c(clusterVarlist,
      "dropvar", "getsFun", "ols", "infocrit", "diagnostics")
    if(!is.null(user.diagnostics)){
      clusterVarlist <- c(clusterVarlist, user.diagnostics$name)
    }
    if(!is.null(user.estimator)){
      clusterVarlist <- c(clusterVarlist, user.estimator$name)
    }
    if(!is.null(gof.function)){
      clusterVarlist <- c(clusterVarlist, gof.function$name)
    }

    #for the future?: add memory.limit()/memory.size() = max cores check?

  } #end if(!is.null(parallel.options))


  ##-------------------------------
  ## 2 x argument
  ##-------------------------------
  
  ##x is a matrix:
  if( is.matrix(x) ){
    xMatrixName <- deparse(substitute(x))
    x <- list(x=x) #convert to list
    names(x) <- xMatrixName
  }

  ##x is a list of matrices:
  if( is.list(x) ){

    ##ensure matrices are named:
    xMatrixNames <- paste0("X", 1:length(x))
    if( is.null(names(x)) ){
      names(x) <- xMatrixNames
    }else{
      for(i in 1:length(x)){
        if( names(x)[i] %in% c("", NA) ){
          names(x)[i] <- xMatrixNames[i]
        }
      } #close for..loop
    }

    ##handle colnames:
    for(i in 1:length(x)){

      xColNames <- colnames(x[[i]])
      if( is.null(xColNames) ){
        xColNames <- paste0("X", i, ".xreg", 1:NCOL(x[[i]]))
      }
      if( any(xColNames=="") ){
        missing.colnames <- which(xColNames == "")
        for(j in 1:length(missing.colnames)){
          #fixed by Jonas: 
          xColNames[ missing.colnames[j] ] <-
            paste0("X", i, ".xreg", missing.colnames[j])
        }
      }
      xColNames <- make.unique(xColNames)
      colnames(x[[i]]) <- xColNames

    } #end for(..) loop
    
    ##do NOT check that colnames() are unique across matrices,
    ##so that matrices can contain the same regressors
    
  } #end if( is.list(x) )

  ##-------------------------------
  ## 3 blocks and keep arguments
  ##-------------------------------

  ##check blocks:
  if( is.list(blocks) ){
    if( length(x)!=length(blocks) ){
      stop("No. of matrices unequal to length(blocks)")
    }
    blocks.is.list <- TRUE
  }else{
    blocks.is.list <- FALSE
    blocks <- list()
  }

  ##keep is vector:
  if( !is.null(keep) && !is.list(keep) && is.vector(keep) ){
    keeptmp <- keep
    keep <- list()
    keep[[1]] <- keeptmp
    if( length(x)>1 ){
      for(i in 2:length(x)){ keep[[i]] <- integer(0) }
    }
  }

  ##check keep argument, name keep items:
  if( is.list(keep) ){
  
    ##check keep argument:
    if( length(x)!=length(keep) ){
      stop("Length(keep) unequal to no. of matrices in 'x'")
    }
    
    ##name keep items:
    for(i in 1:length(keep)){
      if( length(keep[[i]])>0 ){
        names(keep[[i]]) <- colnames(x[[i]])[ keep[[i]] ]
      } 
    }

    ##name entries in keep list:
    names(keep) <- names(x)
        
  } #end if( is.list(keep) )
 
  
  ##-------------------------------
  ## 4 loop on x matrices
  ##-------------------------------

  ##create list w/union of retained regressors from
  ##each x matrix:
  xUnionOfModels <- list() 
  
  ##loop on x:
  for(i in 1:length(x)){

    ##add entry i to list:
    xUnionOfModels[[i]] <- integer(0)
    
    ##blocks:
    if( !blocks.is.list ){

      y.n <- NROW(y)
      ncol.adj <- NCOL(x[[i]])

      ##determine no. of blocks:
      if( is.null(no.of.blocks) ){
        blockratio.value <- ncol.adj/(ratio.threshold*ncol.adj)
        blocksize.value <-
          ncol.adj/min(y.n*ratio.threshold, max.block.size)
        no.of.blocks <- max(2,blockratio.value,blocksize.value)
        no.of.blocks <- ceiling(no.of.blocks)
        no.of.blocks <- min(ncol.adj, no.of.blocks) #ensure blocks < NCOL
      }

      ##make partitions:
      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for(j in 1:no.of.blocks){
        if( blocksize*j <= ncol.adj ){
          partitions.t2[j] <- blocksize*j
        }
      }      
      ##check if last block contains last regressor:
      if( partitions.t2[ length(partitions.t2) ] < ncol.adj ){
        partitions.t2 <- c(partitions.t2, ncol.adj)
      }
      blocksadj <- length(partitions.t2)
      partitions.t1 <- partitions.t2 + 1
      partitions.t1 <- c(1, partitions.t1[ -blocksadj ])

      ##finalise:
      tmp <- list()
      for(j in 1:blocksadj){
        tmp[[j]] <- partitions.t1[j]:partitions.t2[j]
      }
      blocks[[i]] <- tmp

    } #end if(!blocks.is.list)

#    ##name the blocks:
#    names(blocks[[i]]) <- paste0("block", 1:length(blocks[[i]]))

    ##add keep entries to blocks:
    if( length(keep[[i]])>0 ){
      for(j in 1:length(blocks[[i]])){
        blocks[[i]][[j]] <- sort(union(keep[[i]], blocks[[i]][[j]]))
      }
    }
            
    ##xkeep argument:
    if( length(keep[[i]])==0 ){
      xkeep <- NULL
    }else{ xkeep <- names(keep[[i]]) }
    
    ##make blocks function for lapply/parLapply:
    XblocksFun <- function(j, i, x, blocks,
      parallel.options, y, untransformed.residuals, user.estimator,
      t.pval, wald.pval, do.pet, ar.LjungB, arch.LjungB,
      normality.JarqueB, user.diagnostics, gof.function, gof.method,
      xkeep, include.gum, include.1cut, include.empty, max.paths,
      turbo, force.invertibility, tol, LAPACK, max.regs,
      print.searchinfo){

      ##check if block contains 1 regressor:
      if( length(blocks[[i]][[j]])==1 ){
        tmp <- colnames(x[[i]])[ blocks[[i]][[j]] ]
        mX <- cbind(x[[i]][, blocks[[i]][[j]] ])
        colnames(mX) <- tmp
      }else{
        mX <- cbind(x[[i]][, blocks[[i]][[j]] ])
      }

      ##apply dropvar:
      if( force.invertibility ){
        mX <- dropvar(mX, tol=tol, LAPACK=LAPACK, silent=TRUE)
      }

      ##set xkeep argument:
      if( !is.null(xkeep) ){
        xkeep <- which( colnames(mX) %in% xkeep )
      }
      
      ##print info:
      if( is.null(parallel.options) ){
        if(print.searchinfo){
          message("\n", appendLF=FALSE)
          message(names(x)[i],
            " block ", j, " of ", length(blocks[[i]]), ":",
            appendLF=TRUE)
        }
      }

      ##do gets inside XblocksFun:
      getsx <- getsFun(y, mX,
        untransformed.residuals=untransformed.residuals,
        user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, gof.function=gof.function,
        gof.method=gof.method, keep=xkeep, include.gum=include.gum,
        include.1cut=include.1cut, include.empty=include.empty,
        max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
        max.regs=max.regs, print.searchinfo=print.searchinfo,
        alarm=FALSE)
      if( !is.null(getsx$messages) ){ message(getsx$messages) }

      if( is.null(getsx$specific.spec) ){
        xSpecificmodels <- NULL
      }else{
        xSpecificmodels <- names(getsx$specific.spec)
      }

      ##return
      return(xSpecificmodels)
      ##for the future?: return(getsx) - careful! this would change
      ##the subsequent code substantially

    } #close XblocksFun

    ##call XblocksFun/do gets on each block: NO parallel computing
    if( is.null(parallel.options) ){
      xSpecificmodels <- lapply(1:length(blocks[[i]]),
        XblocksFun, i, x, blocks, parallel.options,
        y, untransformed.residuals, user.estimator, t.pval, wald.pval,
        do.pet, ar.LjungB, arch.LjungB, normality.JarqueB,
        user.diagnostics, gof.function, gof.method, xkeep, include.gum,
        include.1cut, include.empty, max.paths, turbo,
        force.invertibility, tol, LAPACK, max.regs, print.searchinfo)
    }
      
    ##call XblocksFun/do gets on each block: WITH parallel computing
    if( !is.null(parallel.options) ){

      ##print info:
      if(print.searchinfo){
        message("\n", appendLF=FALSE)
        message("Preparing parallel computing...",
          appendLF=TRUE)
        message(names(x)[i],
          " blocks to search in parallel: ", length(blocks[[i]]),
          appendLF=TRUE)
        message("Searching...", appendLF=TRUE)
      }

      blocksClust <- makeCluster(clusterSpec, outfile="") #make cluster
      clusterExport(blocksClust, clusterVarlist,
        envir=.GlobalEnv) #idea for the future?: envir=clusterEnvir
      xSpecificmodels <- parLapply(blocksClust,
        1:length(blocks[[i]]), XblocksFun, i, x,
        blocks, parallel.options, y, untransformed.residuals,
        user.estimator, t.pval, wald.pval, do.pet, ar.LjungB,
        arch.LjungB, normality.JarqueB, user.diagnostics, gof.function,
        gof.method, xkeep, include.gum, include.1cut, include.empty,
        max.paths, turbo, force.invertibility, tol, LAPACK, max.regs,
        print.searchinfo)
      stopCluster(blocksClust)

    } #end if( parallel computing )

    ##union of retained variables:
    ##------------------------------------

    ##union of retained variables (names):
    xNames <- NULL
    if( length(xSpecificmodels)>0 ){
      #which variables retained?:
      for(j in 1:length(xSpecificmodels)){
        #check if non-empty:
        if( !is.null(xSpecificmodels[[j]]) ){
          xNames <- union(xNames, xSpecificmodels[[j]])
        }
      }
    }

    ##NOT do gets of union:
    if( gets.of.union==FALSE ){ xUnionOfModels[[i]] <- xNames }
    
    ##DO gets of union:
    if( gets.of.union==TRUE && length(xNames)>0 ){

      if( print.searchinfo ){
        message("\n", appendLF=FALSE)
        message("GETS of union of retained ",
          names(x)[i], " variables... ",
          appendLF=TRUE)
        message("\n", appendLF=FALSE)
      }
  
      ##build regressor matrix:
      mX <- cbind(x[[i]][,xNames])
      colnames(mX) <- xNames
      if( force.invertibility ){
        mX <- dropvar(mX, tol=tol, LAPACK=LAPACK, silent=TRUE)
      }
      
      ##build xkeep:
      if( !is.null(keep) ){
        xkeep <- NULL
        for(j in 1:length(keep)){
          xkeep <- union(xkeep, names(keep[[j]]))
        }
        xkeep <- which( colnames(mX) %in% xkeep)
      }
      
      ##do gets of union:
      getsx <- getsFun(y, mX,
        untransformed.residuals=untransformed.residuals,
        user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, gof.function=gof.function,
        gof.method=gof.method, keep=xkeep, include.gum=include.gum,
        include.1cut=include.1cut, include.empty=include.empty,
        max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
        max.regs=max.regs, print.searchinfo=print.searchinfo,
        alarm=FALSE)
      if( !is.null(names(getsx$specific.spec)) ){
        xUnionOfModels[[i]] <- names(getsx$specific.spec)
      }

    } #end if( do gets )

  } #end for(i) loop (on x matrices)

  ##add names:
  names(blocks) <- names(x)
  names(xUnionOfModels) <- names(x)


  ##-------------------------------
  ## 5 make return object:
  ##-------------------------------

  result$y <- y
  result$x <- list()
  for(i in 1:length(x)){ result$x[[i]] <- colnames(x[[i]]) }
  names(result$x) <- names(x)  
  result$blocks <- blocks
  result$keep <- keep
  result$specific.spec <- xUnionOfModels
  result$time.finished <- date()
  if(alarm){ alarm() }
  return(result)
  
} #close blocksFun()


####################################################
##2 ARX FUNCTIONS
####################################################

##==================================================
##Estimate AR-X model with log-ARCH-X errors
arx <- function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
  zero.adj=0.1, vc.adj=TRUE,
  vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
  user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, plot=NULL)
{
  ### ARGUMENTS: ###########

  vcov.type <- match.arg(vcov.type)

  ##regressand, regressors:
  tmp <- regressorsMean(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg,
    return.regressand=TRUE, return.as.zoo=TRUE,
    na.trim=TRUE,
    na.omit=FALSE)
  
  ##aux: auxiliary list, also used by getsm/getsv
  aux <- list()
  aux$y <- coredata(tmp[,1])
  aux$y.n <- length(aux$y)
  aux$y.name <- colnames(tmp)[1]
  aux$y.index <- index(tmp)
  if( NCOL(tmp)>1 ){
    aux$mX <- cbind(coredata(tmp[,-1]))
    aux$mXnames <- colnames(tmp)[-1]
    colnames(aux$mX) <- NULL
    aux$mXncol <- NCOL(aux$mX)
  }

  ##modify vxreg:
  if( !is.null(vxreg) ){
    ##time:
    vxreg <- as.zoo(cbind(vxreg))
    vxreg <- window(vxreg, start=aux$y.index[1],
      end=aux$y.index[length(aux$y.index)])
    ##colnames:
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste0("vxreg", 1:NCOL(vxreg))
    }
    if( any(vxreg.names == "") ){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[missing.colnames[i]] <- paste0("vxreg", i)
      }
    }
    colnames(vxreg) <- vxreg.names
    ##add to aux (necessary for getsm/getsv):
    aux$vxreg <- vxreg #note: NROW(vxreg)!=NROW(vX) is possible
  }

  ##determine qstat.options:
  if(is.null(qstat.options)){
    if(is.null(ar)){ar.lag <- 1}else{ar.lag <- max(ar)+1}
    if(is.null(arch)){arch.lag <- 1}else{arch.lag <- max(arch)+1}
    qstat.options <- c(ar.lag, arch.lag)
  }
  
  ##info for getsm/getsv functions
  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$user.estimator <- user.estimator
  aux$user.diagnostics <- user.diagnostics
  aux$tol <- tol
  aux$LAPACK <- LAPACK

  ### INITIALISE ##########

  sysCall <- sys.call()
  #for the future: make sure the following objects are part of the out-list?
  vcov.var <- NULL #make sure this object exists
  variance.results <- NULL #make sure this object exists

  #### for the future regarding user.estimator: check if
  #### user.estimator$spec is NULL, "mean", "variance" or "both"
  #### in order to determine what kind of estimator it is

  ##check if mean and log-garch spec:
  meanSpec <- !is.null(aux$mX)
  varianceSpec <- if( vc==FALSE && is.null(arch)
    && is.null(asym) && is.null(log.ewma)
    && is.null(vxreg) ){ FALSE }else{ TRUE }

  #### DEFAULT ESTIMATOR ###############

  if( is.null(user.estimator) ){

    ##estimate:
    estMethod <- which(vcov.type==c("none", "none", "ordinary",
      "white", "newey-west"))
    varianceSpecArg <- NULL
    if( varianceSpec ){
      ##note: vc must be TRUE
      varianceSpecArg <- list(vc=TRUE, arch=arch, asym=asym,
        log.ewma=log.ewma, vxreg=vxreg)
    }
    out <- ols(aux$y, aux$mX, tol=tol, LAPACK=LAPACK, method=estMethod,
      variance.spec=varianceSpecArg)

    ##delete some unneeded entries:
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

    ##re-organise stuff related to mean spec:
    colnames(out$vcov) <- aux$mXnames
    rownames(out$vcov) <- aux$mXnames
    outNames <- names(out)
    whereIs <- which(outNames=="vcov")
    if( length(whereIs) > 0 ){ names(out)[whereIs] <- "vcov.mean" }
    whereIs <- which(outNames=="fit")
    names(out)[whereIs] <- "mean.fit"

    ##if no variance spec:
    if( varianceSpec==FALSE ){
      out$var.fit <- rep(out$sigma2, aux$y.n)
      out$std.residuals <- out$residuals/sqrt(out$sigma2)
      aux$loge2.n <- aux$y.n #same as out$n, change?
      aux$vc <- FALSE #needed for specific in getsm()
    }

    ##if variance spec:
    if( varianceSpec ){

      ##aux: info for getsm and getsv:
      aux$vc <- TRUE #obligatory if varianceSpec
      aux$zero.adj <- zero.adj
      aux$vc.adj <- vc.adj
      aux$loge2 <- out$regressorsVariance[,1]
      aux$loge2.n <- length(aux$loge2)
      aux$vX <- cbind(out$regressorsVariance[,-1])
      aux$vXnames <- colnames(out$regressorsVariance)[-1]
      colnames(aux$vX) <- NULL
      aux$vXncol <- NCOL(aux$vX)
      aux$arch <- arch
      aux$asym <- asym
      aux$log.ewma <- log.ewma
      out$regressorsVariance <- NULL #delete, not needed anymore

      ##re-organise stuff related to variance spec:
      s.e. <- sqrt(as.vector(diag(out$vcov.var)))
      tmpdf <- aux$loge2.n - length(out$var.coefficients)
      tmpvcov <- as.matrix(out$vcov.var[-1,-1])
      colnames(tmpvcov) <- aux$vXnames[-1]
      rownames(tmpvcov) <- aux$vXnames[-1]
      t.stat <- out$var.coefficients/s.e.
      p.val <- pt(abs(t.stat), tmpdf, lower.tail=FALSE)*2
      t.stat[1] <- ((out$var.coefficients[1]-out$Elnz2)^2)/s.e.[1]^2
      p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
      out$var.coefficients[1] <- out$var.coefficients[1] - out$Elnz2
      out$n <- aux$loge2.n
      out$vcov.var <- tmpvcov
      out$variance.results <-
        as.data.frame(cbind(out$var.coefficients, s.e., t.stat, p.val))
      colnames(out$variance.results) <- c("coef", "std.error", "t-stat", "p-value")
      rownames(out$variance.results) <- aux$vXnames
      out$var.coefficients <- NULL

    } #close if( varianceSpec )
        
  } #close if( is.null(user.estimator) )
  

  #### USER-DEFINED ESTIMATOR ###############

  if( !is.null(user.estimator) ){

    ##make user-estimator argument:
    if( is.null(user.estimator$envir) ){ user.estimator$envir <- .GlobalEnv }
    userEstArg <- user.estimator
    userEstArg$name <- NULL
    userEstArg$envir <- NULL
    if( length(userEstArg)==0 ){ userEstArg <- NULL }

    ##user-defined estimator:
    if( is.null(user.estimator$envir) ){
      out <- do.call(user.estimator$name,
        c(list(aux$y,aux$mX), userEstArg))
    }else{
      out <- do.call(user.estimator$name,
        c(list(aux$y,aux$mX), userEstArg), envir=user.estimator$envir)
#        c(list(y=aux$y,x=aux$mX), userEstArg), envir=user.estimator$envir)
    }
    
#delete?:
    ##just in case...:
    if( is.null(out$vcov) && !is.null(out$vcov.mean) ){
      out$vcov <- out$vcov.mean
    }

  } #end if( user.estimator )


  ### OUTPUT: ######################

  ##mean estimation results (a data frame):
  if( meanSpec ){
    if( !is.null(out$vcov) ){
      coefvar <- out$vcov
    }else{
      coefvar <- out$vcov.mean
    }
    stderrs <- sqrt(diag(coefvar))
    t.stat <- out$coefficients/stderrs
    p.val <- pt(abs(t.stat), out$df, lower.tail=FALSE)*2
    out$mean.results <- as.data.frame(cbind(out$coefficients,
      stderrs, t.stat, p.val))
    colnames(out$mean.results) <- c("coef", "std.error",
      "t-stat", "p-value")
    rownames(out$mean.results) <- aux$mXnames
  } #end if(meanSpec)

  ##diagnostics:
  if( any( names(out) %in% c("residuals", "std.residuals") ) ){
    ar.LjungBarg <- c(qstat.options[1],0)
    arch.LjungBarg <- c(qstat.options[2],0)
    if( identical(normality.JarqueB,FALSE) ){
      normality.JarqueBarg <- NULL
    }else{
      normality.JarqueBarg <- as.numeric(normality.JarqueB)
    }
  }else{
    ar.LjungBarg <- NULL
    arch.LjungBarg <- NULL
    normality.JarqueBarg <- NULL
  }
  out$diagnostics <- diagnostics(out,
    ar.LjungB=ar.LjungBarg, arch.LjungB=arch.LjungBarg,
    normality.JarqueB=normality.JarqueBarg,
    user.fun=user.diagnostics, verbose=TRUE)
    
  ##add zoo-indices:
  if(!is.null(out$mean.fit)){
    out$mean.fit <- zoo(out$mean.fit, order.by=aux$y.index)
  }
  if(!is.null(out$residuals)){
    out$residuals <- zoo(out$residuals, order.by=aux$y.index)
  }
  if(!is.null(out$var.fit)){
    out$var.fit <- zoo(out$var.fit, order.by=aux$y.index)
  }
  if(!is.null(out$ustar.residuals)){
    out$ustar.residuals <- zoo(out$ustar.residuals, order.by=aux$y.index)
  }
  if(!is.null(out$std.residuals)){
    out$std.residuals <- zoo(out$std.residuals, order.by=aux$y.index)
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

} #close arx() function

##==================================================
coef.arx <- function(object, spec=NULL, ...)
{
  ##spec argument:
  if(is.null(spec)){
#OLD:
#    spec <- switch(as.character(object$call)[1],
#      arx="both", getsm="mean", getsv="variance")
    spec <- "both"
  }else{
    specType <- c("mean", "variance", "both")
    whichType <- charmatch(spec, specType)
    spec <- specType[ whichType ]
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
  if( !classType %in% c("arx", "gets") ){
    stop("object not of class 'arx' or 'gets'")
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
gets.arx <- function(x, spec=NULL, ...)
{
  ##determine spec:
  if(is.null(spec)){
    if( !is.null(x$mean.results) ){ spec <- "mean" }
    if( is.null(x$mean.results)
      && !is.null(x$variance.results) ){ spec <- "variance" }
  }else{
    specType <- c("mean", "variance")
    whichType <- charmatch(spec, specType)
    spec <- specType[ whichType ]  
  }
  
  ##do the gets modelling:
  if( spec=="mean" ){
    result <- getsm(x, ...)
  }else{
    result <- getsv(x, ...)
  }

  ##return result:
  return(result)
    
} #close gets.arx()

##==================================================
logLik.arx <- function(object, ...)
{
  ## in the future: add a df.method argument with
  ## optional values "mean-coefficients" (default),
  ## "variance-coefficients" and "both"?
  
  result <- object$logl
  if(is.null(result)){
    result <- numeric(0)
    warning("'object$logl' is NULL")
  }else{
    attr(result, "df") <- length(object$coefficients)
    attr(result, "nobs") <- object$n
  }
  class(result) <- "logLik"
  return(result)

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
    message("No estimated model, so no plot produced")
  }
  if(doPlot && !is.null(x$aux$user.estimator) ){
    doPlot <- FALSE
    message("User defined estimation, so no plot produced")
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
  n.sim=5000, innov=NULL, probs=NULL, ci.levels=NULL,
  quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
  plot.options=list(), ...)
{

  ## contents:
  ## 0 initialise
  ## 1 simulate innov
  ## 2 variance predictions
  ## 3 mean predictions
  ## 4 probs (quantiles)
  ## 5 newindex
  ## 6 if plot=TRUE
  ## 7 if return=TRUE

  ##-----------------------
  ## 0 initialise
  ##-----------------------

  ##name of object:
  objectName <- deparse(substitute(object))

  ##check n.ahead:
  if(n.ahead < 1){ stop("n.ahead must be 1 or greater") }

  ##determine spec argument:
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
  
  ##what needs to be predicted?
  predictMean <- switch(spec, "mean"=TRUE, "both"=TRUE,
    "variance"=FALSE)
  predictVariance <- switch(spec, "mean"=FALSE, "both"=TRUE,
    "variance"=TRUE)
  
  ##is there a mean spec?
  coefs <- as.numeric(coef.arx(object, spec="mean"))
  if(length(coefs)>0){ specMean <- TRUE }else{ specMean <- FALSE }

  ##is there a variance spec?
  coefs <- as.numeric(coef.arx(object, spec="variance"))
  if(length(coefs)>0){ specVar <- TRUE }else{ specVar <- FALSE }

  ##determine plot-argument:
  plotArg <- plot
  if( is.null(plotArg) ){
    plotArg <- getOption("plot")
    if( is.null(plotArg) ){ plotArg <- FALSE }
  }

  ##probs argument:
  if( !is.null(probs) ){
    if( any(probs <= 0) || any(probs >= 1) ){
      stop("the values of 'probs' must be between 0 and 1")
    }
    probs <- union(probs,probs) #ensure values are distinct/not repeated
    probs <- probs[order(probs, decreasing=FALSE)] #re-order to increasing
  }
  probsArg <- probs

  ##ci.levels argument:
  if( is.null(ci.levels) && plotArg==TRUE ){
    if( class(object)=="isat" ){ ##if isat:
      ciLevelsArg <- c(0.68,0.95)
    }else{ ##if not isat:
      ciLevelsArg <- c(0.5,0.9)    
    }
  }else{
    ciLevelsArg <- ci.levels
  }
  if( !is.null(ciLevelsArg) ){
    if( any(ciLevelsArg <= 0) || any(ciLevelsArg >= 1) ){
      stop("'ci.levels' must be between 0 and 1")
    }
    ciLevelsArg <- union(ciLevelsArg, ciLevelsArg) #ensure levels are distinct/not repeated
    ciLevelsArg <- ciLevelsArg[order(ciLevelsArg, decreasing=FALSE)] #ensure levels are increasing
    ciLower <- (1-ciLevelsArg)/2
    ciLower <- ciLower[order(ciLower, decreasing=FALSE)]    
    ciUpper <- ciLevelsArg + (1-ciLevelsArg)/2
    ciUpper <- ciUpper[order(ciUpper, decreasing=TRUE)]    
    probsArg <- union(probsArg, c(ciLower,ciUpper) ) #add to probsArg
    probsArg <- probsArg[order(probsArg, decreasing=FALSE)] #ascending
  }
  
  ##are simulations of innov needed?
  doSimulations <- FALSE
  if( !is.null(probsArg) ){ doSimulations <- predictVariance <- TRUE }
  if( predictVariance && specVar ){ doSimulations <- TRUE }


  ##-----------------------
  ## 1 simulate innov
  ##-----------------------

  ##simulate:
  mZhat <- NULL
  if(doSimulations){

    ##bootstrap (innov not user-provided):
    if(is.null(innov)){
      zhat <- coredata(na.trim(object$std.residuals))
      if(specVar){
        where.zeros <- which(zhat==0)
        if(length(where.zeros)>0){ zhat <- zhat[-where.zeros] }
      }
      draws <- runif(n.ahead*n.sim, min=0.5+.Machine$double.eps,
                     max=length(zhat)+0.5+.Machine$double.eps)
      draws <- round(draws, digits=0)
      zhat <- zhat[draws]
    }
    
    ##user-provided innov:
    if(!is.null(innov)){
      if(length(innov)!=n.ahead*n.sim){ stop("length(innov) must equal n.ahead*n.sim") }
      if(specVar){
        if(any(innov==0)){ stop("'innov' cannot contain zeros") }
      }
      zhat <- as.numeric(innov)
    }

    ##matrix of innovations:
    mZhat <- matrix(zhat,n.ahead,n.sim)
    colnames(mZhat) <- paste0("mZhat.", seq(1,n.sim))
  
  } #end if(doSimulations)
  
  
  ##-----------------------
  ## 2 variance predictions
  ##-----------------------

  sd2hat <- NULL #variance predictions
  mEpsilon <- NULL #matrix of simulated errors
  if(predictVariance){
    
    ##there is no variance specification:
    if(specVar==FALSE){
      sigmahat <- sigma.arx(object)
      sd2hat <- rep(sigmahat^2, n.ahead)
    }

    ##there is a variance specification:
    if(specVar==TRUE){

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
      if(backcastMax>0){
        lnsd2hat[1:backcastMax] <- 
          lnsd2Fit[c(length(lnsd2Fit)-backcastMax+1):length(lnsd2Fit)]
      }
      mLnsd2Hat <- matrix(NA, lnsd2hat.n, n.sim) #matrix of lnsd2 predictions
      mLnsd2Hat[,1:NCOL(mLnsd2Hat)] <- lnsd2hat  #fill with backcast values
  
      ##prepare lnz2:
      lnz2hat <- rep(NA, n.ahead + backcastMax)
      lnz2hat.n <- length(lnz2hat)
      lnz2Fit <- coredata(object$ustar.residuals) + Elnz2hat
      if(backcastMax>0){
        lnz2hat[1:backcastMax] <- lnz2Fit[c(length(lnz2Fit)-backcastMax+1):length(lnz2Fit)]
      }
      mLnz2Hat <- matrix(NA, lnz2hat.n, n.sim)
      mLnz2Hat[,1:NCOL(mLnz2Hat)] <- lnz2hat
      mZhat2 <- mZhat^2 #mZhat from section 1
      mLnz2Hat[c(backcastMax+1):NROW(mLnz2Hat),] <- log(mZhat2)
      vEpsilon2 <- rep(NA, n.ahead+backcastMax) #needed for log(ewma) term(s)
      if(backcastMax>0){
        vEpsilon2[1:backcastMax] <- as.numeric(object$residuals[c(length(object$residuals)-backcastMax+1):length(object$residuals)]^2)
        mZhat2 <- rbind(matrix(NA,backcastMax,NCOL(mZhat2)),mZhat2)
      }
   
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
      mSd2Hat <- exp( mLnsd2Hat[c(lnsd2hat.n-n.ahead+1):lnsd2hat.n,] )
      if(n.ahead==1){ mSd2Hat <- rbind(mSd2Hat) } #rbind() needed when n.ahead=1
      sd2hat <- as.vector(rowMeans(mSd2Hat))
  
    } #end if(specVar==TRUE)

    ##matrix of errors:
    if(!is.null(mZhat)){
      mEpsilon <- sqrt(sd2hat)*mZhat
      colnames(mEpsilon) <- paste0("mEpsilon.", seq(1,n.sim))
    }

  } #end if(predictVariance)


  ##-----------------------
  ## 3 mean predictions
  ##-----------------------

  yhat <- NULL #mean predictions
  mY <- NULL #matrix of simulated y's
  if(predictMean){
    
    ##there is no mean specification:
    if(specMean==FALSE){
      yhat <- rep(0, n.ahead)
      if(!is.null(mEpsilon)){
        mY <- yhat + mEpsilon
        colnames(mY) <- paste0("mY.", seq(1,n.sim))
      }
    }

    ##there is a mean specification:
    if(specMean==TRUE){

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
      if( !is.null(object$call$ewma) ){
        ewmaEval <- eval(object$call$ewma)
        if(is.list(ewmaEval)){ ewmaEval <- ewmaEval$length }
        ewmaIndx <- 1:length(ewmaEval) + max(arIndx)
        ewmaMax <- max(ewmaEval)
        ewmaCoefs <- as.numeric(coefs[ewmaIndx])
      }
  
      ##backcast length:
      backcastMax <- max(arMax,ewmaMax)

      ##mxreg:
      mxreghat <- rep(0, n.ahead + backcastMax)
      if(!is.null(object$call$mxreg)){

        ##check newmxreg:
        if(is.null(newmxreg)){ stop("'newmxreg' is NULL") }
        if(NROW(newmxreg)!=n.ahead){ stop("NROW(newmxreg) must equal n.ahead") }
  
        ##newmxreg:
        newmxreg <- coredata(cbind(as.zoo(newmxreg)))
        colnames(newmxreg) <- NULL
  
        ##mxreghat:
        mxregIndx <- c(max(ewmaIndx)+1):length(coefs)
        mxreghat <-  newmxreg %*% as.numeric(coefs[mxregIndx])
        mxreghat <- c(rep(0,backcastMax),mxreghat)
  
      } ##end mxreg
  
      ##prepare prediction:
      yhat <- rep(NA, n.ahead + backcastMax)
      yhat.n <- length(yhat)
      if(backcastMax>0) {
        ##actual y-values:
        yhat[1:backcastMax] <- tail(object$aux$y, n=backcastMax)
      }

      ##prepare ewma:
      if(ewmaMax>0){
        mEwmaHat <- matrix(NA, n.ahead+backcastMax, length(ewmaCoefs))
        colnames(mEwmaHat) <- object$aux$mXnames[ewmaIndx]
        mEwmaHat[1:backcastMax,] <- object$aux$mX[c(NROW(object$aux$mX)-backcastMax+1):NROW(object$aux$mX),ewmaIndx]
        mEwmaHat <- as.matrix(mEwmaHat)
      }

      ##predict yhat:
      arTerm <- 0
      ewmaTerm <- 0
      for(i in c(backcastMax+1):yhat.n){
        if( arMax>0 ){ arTerm <- sum(arCoefs*yhat[c(i-1):c(i-arMax)]) }
        if( ewmaMax>0 ){
          for(k in 1:NCOL(mEwmaHat)){
            mEwmaHat[i,k] <- mean( yhat[c(i-ewmaEval[k]):c(i-1)] )
          }
          ewmaTerm <- sum( coefs[ewmaIndx] * mEwmaHat[i,] )
        }
        yhat[i] <- mconst + arTerm + ewmaTerm + mxreghat[i]
      } #end loop
  
      ##out:
      yhat <- yhat[c(yhat.n-n.ahead+1):yhat.n]
  
      ##simulate mY?:
      if( !is.null(mEpsilon) ){
        
        ##loop on j:
        mY <- matrix(NA, NROW(mEpsilon), NCOL(mEpsilon))
        for(j in 1:NCOL(mEpsilon)){

          ##prepare prediction no. j:
          yhatadj <- rep(NA, n.ahead + backcastMax)
          if(backcastMax>0) {
            ##actual y-values:
            yhatadj[1:backcastMax] <- tail(object$aux$y, n=backcastMax)
          }

          ##prepare ewma:
          if(ewmaMax>0){
            mEwmaHat <- matrix(NA, n.ahead+backcastMax, length(ewmaCoefs))
            colnames(mEwmaHat) <- object$aux$mXnames[ewmaIndx]
            mEwmaHat[1:backcastMax,] <- object$aux$mX[c(NROW(object$aux$mX)-backcastMax+1):NROW(object$aux$mX),ewmaIndx]
            mEwmaHat <- as.matrix(mEwmaHat)
          }

          ##predict yhatadj:
          arTerm <- 0
          ewmaTerm <- 0
          for(i in c(backcastMax+1):yhat.n){
            if( arMax>0 ){ arTerm <- sum(arCoefs*yhatadj[c(i-1):c(i-arMax)]) }
            if( ewmaMax>0 ){
              for(k in 1:NCOL(mEwmaHat)){
                mEwmaHat[i,k] <- mean( yhatadj[c(i-ewmaEval[k]):c(i-1)] )
              }
              ewmaTerm <- sum( coefs[ewmaIndx] * mEwmaHat[i,] )
            }
            yhatadj[i] <- mconst + arTerm + ewmaTerm +
              mxreghat[i] + mEpsilon[c(i-backcastMax),j]
          } #end loop
  
          ##store the simulation of yhatadj:
          mY[,j] <- yhatadj[c(yhat.n-n.ahead+1):yhat.n]

        } #end for(j)

        ##add colnames:
        colnames(mY) <- paste0("mY.", seq(1,n.sim))

      } #end if( not null(mEpsilon) )
    
    } #end if(specMean==TRUE)

  } #end if(predictMean)


  ##-----------------------
  ## 4 probs (quantiles)
  ##-----------------------

  ##mean:
  mMeanQs <- NULL
  if( predictMean && !is.null(probsArg) ){
    mMeanQs <- matrix(NA, n.ahead, length(probsArg))
    for(i in 1:NROW(mY)){
      mMeanQs[i,] <- quantile(mY[i,], probs=probsArg, type=quantile.type)
    }
    colnames(mMeanQs) <- paste0(probsArg)
  }
   
  ##variance:
  mVarianceQs <- NULL
  if( predictVariance && specVar && !is.null(probsArg) ){
    mVarianceQs <- matrix(NA, n.ahead, length(probsArg))
    for(i in 1:NROW(mSd2Hat)){
      mVarianceQs[i,] <- quantile(mSd2Hat[i,], probs=probsArg, type=quantile.type)
    }
    colnames(mVarianceQs) <- paste0(probsArg)
  }


  ##-----------------------
  ## 5 newindex
  ##-----------------------
  
  ##in-sample:
  yInSample <- zoo(object$aux$y, order.by=object$aux$y.index)

  #newindex user-provided:
  if( !is.null(newindex) ){
    yAsRegular <- FALSE
    if( n.ahead!=length(newindex) ){
      stop("length(newindex) must equal 'n.ahead'")
    }
    newindexInSample <- any( newindex %in% object$aux$y.index )
  }else{ newindexInSample <- FALSE }

  #in-sample index regular:
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

  ##neither user-provided nor regular:
  if( is.null(newindex) ){ newindex <- 1:n.ahead }

  ##add index to results:
  if(!is.null(mZhat)){ mZhat <- zoo(mZhat, order.by=newindex) }
  if(!is.null(sd2hat)){ sd2hat <- zoo(sd2hat, order.by=newindex) }
  if(!is.null(mEpsilon)){ mEpsilon <- zoo(mEpsilon, order.by=newindex) }
  if(!is.null(yhat)){ yhat <- zoo(yhat, order.by=newindex) }
  if(!is.null(mY)){ mY <- zoo(mY, order.by=newindex) }
  if(!is.null(mMeanQs)){ mMeanQs <- zoo(mMeanQs, order.by=newindex) }
  if(!is.null(mVarianceQs)){ mVarianceQs <- zoo(mVarianceQs, order.by=newindex) }

  
  ##-----------------------
  ## 6 if plot=TRUE
  ##-----------------------

  ##check special case:
  if( plotArg && spec=="variance" && specVar==FALSE ){
    message("Set 'vc = TRUE' to enable a plot of the variance predictions")
    plotArg <- FALSE #change argument
  }

  ##check another special case (plot.zoo does not work if
  ##the index is not unique):
  if( plotArg && newindexInSample==TRUE ){
    message("'newindex' not entirely out-of-sample, so no plot produced")
    plotArg <- FALSE #change argument
  }
  
  ##idea?: the special case where the out-of-sample index
  ##is not of the same type as the in-sample index. for 
  ##regular zoo-series, this can possibly be checked by checking
  ##whether the numeric delta is the same both in-sample and
  ##out-of-sample.
  
  ##way to check for this is 
  ##plot?:
  if( plotArg ){

    ##some of the plot.options:
    ##-------------------------

    ##remember: both ylab and hlines argument can be specified,
    ##even though they do not appear below
                
    ##how many in-sample observations to include in plot?:
    if( is.null(plot.options$keep) ){ plot.options$keep <- 12L }
    if( plot.options$keep < 1 ){
      plot.options$keep <- 1L
      message("'plot.options$keep' changed to 1")
    }
    
    ##if "main" argument:
    if( is.null(plot.options$main) ){
      parMarVals <- c(2.1,3.1,0.6,0.6) #bottom,left,top,right
    }else{
      parMarVals <- c(2.1,3.1,1.5,0.6) #bottom,left,top,right
    }
        
    ##linetype (solid=1, dashed=2, etc.). Order: lty=c(Forecast,Actual)
    if( is.null(plot.options$lty )){ plot.options$lty <- c(1,1) }

    ##linewidth. Order: lwd=c(Forecast,Actual)
    if( is.null(plot.options$lwd) ){ plot.options$lwd <- c(1,1) }
    
    ##colours. Order: col=c(Forecast,Actual)
    if( is.null(plot.options$col) ){ plot.options$col <- c("red","blue") }
    
    ##text for legend:
    if( is.null(plot.options$legend.text) ){
      if( spec == "variance" ){
        plot.options$legend.text <- c("Forecast", "Squared residuals")
      }else{
        plot.options$legend.text <- c("Forecast", "Actual")
      }
    }
    
    ##whether to include retained fitted or not:
    if( is.null(plot.options$fitted) ){ plot.options$fitted <- FALSE }

    ##should predictions start at origin?:
    if( is.null(plot.options$start.at.origin) ){
      plot.options$start.at.origin <- TRUE
    }

    ##add dot at forecast origin?:
    if( is.null(plot.options$dot.at.origin) ){    
      plot.options$dot.at.origin <- TRUE
    }
    
    ##add vertical line at forecast origin?:
    if( is.null(plot.options$line.at.origin) ){
      plot.options$line.at.origin <- FALSE
    }
      
    ##check if( !is.null(shades.of.grey) ):
    if( !is.null(plot.options[["shades.of.grey"]]) ){
      message(
        "argument 'shades.of.grey' has changed name to 'shades',\n",
        "'shades.of.grey' will be deprecated in future versions"
      )
      if( is.null(plot.options[["shades"]]) ){
        plot.options[["shades"]] <- plot.options[["shades.of.grey"]]
      }
    }


    ##start preparing:
    ##----------------
    
    ##select the shades of grey for the ci's:
    if( is.null(plot.options$shades) ){
      shadesOfGrey <- 40:90 #1 to 100 is possible
      shadesOfGrey <- quantile(shadesOfGrey, probs=ciLevelsArg) 
      shadesOfGrey <- shadesOfGrey[length(shadesOfGrey):1] #invert, i.e. last to first
      plot.options$shades <- round(as.numeric(shadesOfGrey))
    }
    greySelection <- paste0("grey", plot.options$shades)
    
    ##make dataForPlot:
    dataForPlot <- matrix(NA, n.ahead, 6)
    colnames(dataForPlot) <- c("MeanActual", "MeanFitted",
      "MeanPrediction", "ResidualsSquared", "VarianceFitted",
      "VariancePrediction")
    if(!is.null(plot.options$newmactual)){
      dataForPlot[1:length(plot.options$newmactual),"MeanActual"] <-
        plot.options$newmactual
    }
    if(!is.null(plot.options$newvactual)){
      dataForPlot[1:length(plot.options$newvactual),"ResidualsSquared"] <-
        plot.options$newvactual
    }
    if(!is.null(yhat)){
      dataForPlot[,"MeanPrediction"] <- coredata(yhat)
    }
    if(!is.null(sd2hat)){
      dataForPlot[,"VariancePrediction"] <- coredata(sd2hat)
    }
    retainedData <- matrix(NA, plot.options$keep, NCOL(dataForPlot))
    colnames(retainedData) <- colnames(dataForPlot)
    retainedData[,"MeanActual"] <-
      tail(coredata(yInSample), n=plot.options$keep)
    retainedData[,"ResidualsSquared"] <-
      tail(coredata(object$residuals)^2, n=plot.options$keep)
    retainedData[,"MeanFitted"] <-
      tail(coredata(object$mean.fit), n=plot.options$keep)
    retainedData[,"VarianceFitted"] <-
      tail(coredata(object$var.fit), n=plot.options$keep)
    if( plot.options$start.at.origin ){ ##let predictions start.at.origin:      
      retainedData[NROW(retainedData),"MeanPrediction"] <- 
        retainedData[NROW(retainedData),"MeanActual"]
      retainedData[NROW(retainedData),"VariancePrediction"] <- 
        retainedData[NROW(retainedData),"VarianceFitted"]
    }
    if( !plot.options$start.at.origin && plot.options$fitted ){
      retainedData[,"MeanPrediction"] <- retainedData[,"MeanFitted"]
    }
    dataForPlot <- rbind(retainedData, dataForPlot)
    tmpIndx <- c(tail(index(yInSample), n=plot.options$keep), newindex)
    dataForPlot <- zoo(dataForPlot, order.by=tmpIndx)
        
    ##if spec="mean" or "both":
    ##-------------------------

    if( spec %in% c("mean","both") ){

      ##create polygon index:
      i1 <- ifelse(plot.options$start.at.origin, 0, 1)   
      polygonIndx <- c(NROW(dataForPlot)-n.ahead+i1):NROW(dataForPlot)
      polygonIndx <- c(polygonIndx, polygonIndx[c(length(polygonIndx):1)])
      polygonIndx <- index(dataForPlot)[polygonIndx]

      ##matrices with the ci's:
      mCiLowerValsMean <- cbind(coredata(mMeanQs[,as.character(ciLower)]))
      colnames(mCiLowerValsMean) <- as.character(ciLower)
      mCiUpperValsMean <- cbind(coredata(mMeanQs[,as.character(ciUpper)]))
      colnames(mCiUpperValsMean) <- as.character(ciUpper)
      mCiUpperValsMean <-
        mCiUpperValsMean[NROW(mCiUpperValsMean):1,] #invert (first to last, last to first)
      if(n.ahead==1){ ##ensure they are still matrices:
        mCiLowerValsMean <- rbind(mCiLowerValsMean)
        mCiUpperValsMean <- rbind(mCiUpperValsMean)
      }
            
      ##add actual value at forecast origin to ci matrices?:
      if( plot.options$start.at.origin ){
        actualValue <- retainedData[NROW(retainedData),"MeanActual"]  
        mCiLowerValsMean <- rbind(actualValue,mCiLowerValsMean)
        mCiUpperValsMean <- rbind(mCiUpperValsMean,actualValue)
      }
      
      ##y-axis (limits):
      if( is.null(plot.options$ylim) ){

        ylimArg <- c(coredata(mCiLowerValsMean[,1]),
          coredata(mCiUpperValsMean[,1]))
        if( plot.options$keep > 0 ){
          ylimArg <- c(ylimArg,
            tail(coredata(yInSample), n=plot.options$keep) )
          if(plot.options$fitted){
            ylimArg <- c(ylimArg,
              tail(coredata(object$mean.fit), n=plot.options$keep))
          }
        }
        if(!is.null(plot.options$newmactual)){
          ylimArg <- c(ylimArg, coredata(plot.options$newmactual))
        }
        ylimArg <- range(ylimArg)
        eps <- abs(ylimArg[2]-ylimArg[1])
        ylimArg[2] <- ylimArg[2] + eps*0.15 #add more space at the top
        ylimArg[1] <- ylimArg[1] - eps*0.05 #add more space at the bottom

      }else{ ylimArg <- plot.options$ylim }
           
      ##get current par-values:
      def.par <- par(no.readonly=TRUE)
  
      ##margins:
      par(mar=parMarVals) 
  
      ##plot actual values in white (i.e. create plot):
      plot.zoo(dataForPlot[,"MeanActual"], xlab="", ylab="",
        main=plot.options$main, lty=plot.options$lty[2],
        col="white", lwd=plot.options$lwd[2], ylim=ylimArg)

      ##add start line?:
      if( plot.options$line.at.origin ){
        startlineIndx <- rep( index(dataForPlot)[plot.options$keep], 2)
        eps <- abs(ylimArg[2]-ylimArg[1])
        startlineVals <- c(ylimArg[1]-eps*0.05, ylimArg[2]/1.2)
        polygon(startlineIndx, startlineVals, col="grey",
          border="grey", lwd=plot.options$lwd)
      }
      
      ##add ci's:
      for(i in 1:length(ciLevelsArg)){
        polygon( polygonIndx,
          c(mCiLowerValsMean[,i],mCiUpperValsMean[,i]),
          col=greySelection[i], border=greySelection[i] )  
      }
                  
      ##add horisontal lines?:
      if(!is.null(plot.options$hlines)){
        abline(h=plot.options$hlines, col="grey", lty=3)
      }

      ##add prediction:
      lines(dataForPlot[,"MeanPrediction"], lty=plot.options$lty[1],
        col=plot.options$col[1], lwd=plot.options$lwd[1])
    
      ##add actual:
      lines(dataForPlot[,"MeanActual"], lty=plot.options$lty[2],
        col=plot.options$col[2], lwd=plot.options$lwd[2], type="l")
        
      ##add fitted (pre-prediction):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"MeanFitted"], lty=plot.options$lty[2],
          col=plot.options$col[1], lwd=plot.options$lwd[1],
          type="l")
      }

      ##add point at forecast origin?:
      if( plot.options$dot.at.origin ){
        points(index(dataForPlot)[NROW(retainedData)],
          retainedData[NROW(retainedData),"MeanActual"],
          pch=19, col=plot.options$col[2], lwd=plot.options$lwd[2])
      }
      
      ##add actual values out-of-sample:
      if( !is.null(plot.options$newmactual) ){
        lines(dataForPlot[,"MeanActual"], lty=plot.options$lty[2],
          col=plot.options$col[2], lwd=plot.options$lwd[2],
          type="l")
      }

      ##add text closer to plot than xlab or ylab would do
      mtextValue <- ifelse(is.null(plot.options$ylab),
        "Mean", plot.options$ylab)
      mtext(mtextValue, side=2, line=2)
  
      ##add plot-legend:
      legend("top", lty=plot.options$lty, col=plot.options$col,
        lwd=plot.options$lwd, legend=plot.options$legend.text,
        bg="white", bty="n")
  
      ##add ci-legend:
      legendArg <- ciLevelsArg[length(ciLevelsArg):1]*100
      legendArg <- paste0(legendArg, "%")      
      legend("topright", lty=c(1,1), lwd=13, bty="n",
        col=greySelection[c(1,length(ciLevelsArg))],
        legend=legendArg[c(1,length(ciLevelsArg))])
      
      ##return to old par-values:
      par(def.par)

    } #end if(spec %in% c("mean","both"))


    ##if spec="variance":
    ##-------------------

    if( spec=="variance" ){

      ##create polygon index:
      i1 <- ifelse(plot.options$start.at.origin, 0, 1)   
      polygonIndx <- c(NROW(dataForPlot)-n.ahead+i1):NROW(dataForPlot)
      polygonIndx <- c(polygonIndx, polygonIndx[c(length(polygonIndx):1)])
      polygonIndx <- index(dataForPlot)[polygonIndx]

      ##matrices with the ci's:
      mCiLowerValsVar <- cbind(coredata(mVarianceQs[,as.character(ciLower)]))
      colnames(mCiLowerValsVar) <- as.character(ciLower)
      mCiUpperValsVar <- cbind(coredata(mVarianceQs[,as.character(ciUpper)]))
      colnames(mCiUpperValsVar) <- as.character(ciUpper)
      mCiUpperValsVar <-
        mCiUpperValsVar[NROW(mCiUpperValsVar):1,] #invert (first to last, last to first)

      ##add fitted value to forecast origin?:
      if( plot.options$start.at.origin ){
        fittedValue <- retainedData[NROW(retainedData),"VarianceFitted"]  
        mCiLowerValsVar <- rbind(fittedValue,mCiLowerValsVar)
        mCiUpperValsVar <- rbind(mCiUpperValsVar,fittedValue)
      }

      ##y-axis (limits):
      if(is.null(plot.options$ylim)){

        ylimArg <- c(coredata(mCiLowerValsVar[,1]),
          coredata(mCiUpperValsVar[,1]))
        if( plot.options$keep > 0 ){
          ylimArg <- c(ylimArg,
            tail(coredata(object$residuals)^2, n=plot.options$keep) )
          if(plot.options$fitted){
            ylimArg <- c(ylimArg,
              tail(coredata(object$var.fit), n=plot.options$keep))
          }
        }
        if(!is.null(plot.options$newvactual)){
          ylimArg <- c(ylimArg, coredata(plot.options$newvactual))
        }
        ylimArg <- range(ylimArg)
        eps <- abs(ylimArg[2]-ylimArg[1])
        ylimArg[2] <- ylimArg[2] + eps*0.15 #add more space at the top
        ylimArg[1] <- ylimArg[1] - eps*0.05 #add more space at the bottom
     
      }else{ ylimArg <- plot.options$ylim }
      
      ##get current par-values:
      def.par <- par(no.readonly=TRUE)
  
      ##margins:
      par(mar=parMarVals) 
  
      ##plot the actual values:
      plot.zoo(dataForPlot[,"ResidualsSquared"], xlab="", ylab="",
        main=plot.options$main, lty=plot.options$lty[2],
        col=plot.options$col[2], lwd=plot.options$lwd,
        ylim=ylimArg)
    
      ##add start line?:
      if( plot.options$line.at.origin ){
        startlineIndx <- rep( index(dataForPlot)[plot.options$keep], 2)
        eps <- abs(ylimArg[2]-ylimArg[1])
        startlineVals <- c(ylimArg[1]-eps*0.05, ylimArg[2]/1.2)
        polygon(startlineIndx, startlineVals, col="grey",
          border="grey", lwd=plot.options$lwd)
      }
  
      ##add ci's:
      for(i in 1:length(ciLevelsArg)){
        polygon( polygonIndx,
          c(mCiLowerValsVar[,i],mCiUpperValsVar[,i]),
          col=greySelection[i], border=greySelection[i] )  
      }

      ##add horisontal lines?:
      if(!is.null(plot.options$hlines)){
        abline(h=plot.options$hlines, col="grey", lty=3)
      }
  
      ##add prediction:
      lines(dataForPlot[,"VariancePrediction"], lty=plot.options$lty[1],
        col=plot.options$col[1], lwd=plot.options$lwd,
        type="l")
    
      ##add fitted (in-sample):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"VarianceFitted"], lty=plot.options$lty[2],
          lwd=plot.options$lwd, col=plot.options$col[1],
          type="l")
      }
  
      ##add point at forecast origin?:
      if(plot.options$dot.at.origin){
        points(index(dataForPlot)[NROW(retainedData)],
          retainedData[NROW(retainedData),"VarianceFitted"], 
          #OLD: fittedValue,
          pch=19, col=plot.options$col[1], lwd=plot.options$lwd)
      }

      ##add actual values of residuals squared out-of-sample:
      if( !is.null(plot.options$newvactual) ){
        lines(dataForPlot[,"ResidualsSquared"], lty=plot.options$lty[2],
          col=plot.options$col[2], lwd=plot.options$lwd,
          type="l")
      }

      ##add text closer to plot than xlab or ylab would do
      mtextValue <- ifelse(is.null(plot.options$ylab),
        "Variance", plot.options$ylab)
      mtext(mtextValue, side=2, line=2)

      ##add plot-legend:
      legend("top", lty=plot.options$lty, col=plot.options$col,
        lwd=plot.options$lwd, legend=plot.options$legend.text,
        bty="n")
      
      ##add ci-legend:
      legendArg <- ciLevelsArg[length(ciLevelsArg):1]*100
      legendArg <- paste0(legendArg, "%")      
      legend("topright", lty=c(1,1), lwd=13, bty="n",
        col=greySelection[c(1,length(ciLevelsArg))],
        legend=legendArg[c(1,length(ciLevelsArg))])
      
      ##return to old par-values:
      par(def.par)

    } #end if(spec %in% c("mean","both"))

  } #end if(plotArg)

      
  ##-----------------------
  ## 7 if return=TRUE
  ##-----------------------

  if(return){

    ##change colnames on quantiles:
    if(!is.null(mMeanQs)){ colnames(mMeanQs) <- paste0("y",probsArg) }
    if(!is.null(mVarianceQs)){ colnames(mVarianceQs) <- paste0("sd2",probsArg) }
    
    ##return everything:
    if(verbose){
      result <- NULL
      if(!is.null(yhat)){ result <- cbind(yhat) }
      if(!is.null(mMeanQs)){
        if(is.null(result)){ result <- mMeanQs }else{ result <- cbind(result,mMeanQs) }
      }
      if(!is.null(mY)){
        if(is.null(result)){ result <- mY }else{ result <- cbind(result,mY) }
      }
      if(!is.null(sd2hat)){
        if(is.null(result)){ result <- sd2hat }else{ result <- cbind(result,sd2hat) }
      }
      if(!is.null(mVarianceQs)){
        if(is.null(result)){ result <- mVarianceQs }else{ result <- cbind(result,mVarianceQs) }
      }
      if(!is.null(mEpsilon)){
        if(is.null(result)){ result <- mEpsilon }else{ result <- cbind(result,mEpsilon) }
      }
      if(!is.null(mZhat)){
        if(is.null(result)){ result <- mZhat }else{ result <- cbind(result, mZhat) }
      }
    } #end if(verbose)
    
    ##do not return everything:
    if(!verbose){

      resultMean <- NULL
      resultVariance <- NULL
      
      ##mean specification:
      if( spec %in% c("mean","both" ) ){
        resultMean <- yhat
        if( !is.null(probs) || !is.null(ci.levels) ){
          resultMean <- cbind(yhat,mMeanQs)
        }
      }

      ##mean specification:
      if( spec %in% c("variance","both" ) ){
        resultVariance <- sd2hat
        if( !is.null(probs) || !is.null(ci.levels) ){
          resultVariance <- cbind(sd2hat,mVarianceQs)
        }
      }

      ##combine:
      if(is.null(resultMean)){ result <- resultVariance }
      if(is.null(resultVariance)){ result <- resultMean }
      if(!is.null(resultMean) && !is.null(resultVariance) ){
        result <- cbind(resultMean,resultVariance)
        colnames(result) <- c("yhat", "sd2hat")
      }
          
    } #end if(!verbose)

    ##return the result:
    return(result)

  } #end if(return)
  
} #close predict.arx  

##==================================================
## print estimation result
print.arx <- function(x, signif.stars=TRUE, ...)
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

  ##create goodness-of-fit matrix:
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
#NEW (suggested by Moritz)?:
#    printCoefmat(x$diagnostics, tst.ind=2,
#      signif.stars=signif.stars, has.Pvalue=TRUE)
    if( !is.null(x$gof) ){
      printCoefmat(x$gof, digits=6, signif.stars=FALSE)
    }
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
    tmpEst <- ols(vY, mX, tol=tol, LAPACK=LAPACK, method=olsMethod)
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
rsquared <- function(object, adjusted=FALSE, ...)
{
  classOK <- class(object) %in% c("arx", "gets", "isat")
  if(!classOK){ message("object not of class 'arx', 'gets' or 'isat'") }
  if( class(object) == "gets" ){
    specType <- switch(as.character(object$call)[1],
      getsm="mean", getsv="variance")
  }
  if( class(object) == "gets" && specType=="variance" ){
    result <- NA
#OLD:
#    Rsquared <- NA
  }else{
    TSS <- sum( (object$aux$y - mean(object$aux$y))^2 )
    residsTrimmed <- na.trim(object$residuals)
    RSS <- sum(residsTrimmed^2)
    Rsquared <- 1 - RSS/TSS
    if(adjusted){
      result <- 1 - (1-Rsquared)*(object$n-1)/(object$n-object$k)
    }else{
      result <- Rsquared
    }
  }
  return(result)
} #close rsquared function

##==================================================
## summarise output
summary.arx <- function(object, ...)
{
  summary.default(object)
} #end summary.arx

##==================================================
## LaTeX code (equation form)
toLatex.arx <- function(object, ...)
{
  printtex(object, ...)
} #end toLatex.arx

##==================================================
VaR <- function(object, level=0.99, type=7, ...)
{
  ##check whether class is valid:
  classType <- class(object)
  if( !classType %in% c("arx", "gets") ){
    stop("object not of class 'arx' or 'gets'")
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
  specOriginal <- spec
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

  ##create result:
  result <- NULL
  
  ##if mean:
  if(spec=="mean"){
    result <- object$vcov.mean
  }

  ##if variance:
  if(spec=="variance"){
    result <- object$vcov.var
  }

#  ##check and change if 0 x 0?:
#  if(all(dim(result)==0)){ result <- NULL }

  ##if user-specified estimator:
  if( !is.null(object$aux$user.estimator) && is.null(specOriginal) ){
    result <- object$vcov
    if( is.null(colnames(result)) ){
      colnames(result) <- names(coef.arx(object))
      rownames(result) <- colnames(result)
    }
  }
  
  return(result)
} #close vcov.arx


####################################################
##3 GETS FUNCTIONS
####################################################

##==================================================
## Multi-path GETS modelling of mean specification
getsm <- function(object, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
  gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
  include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
  turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE)
{
  ## contents:
  ## 1 arguments
  ## 2 gets modelling
  ## 3 estimate specific
  ## 4 output
  
  ##------------------
  ## 1 arguments
  ##------------------
  
  ##check if mean equation:
  if( is.null(object$aux$mX) ){ stop("Mean equation empty") }

  ##check max.paths:
  if( !is.null(max.paths) && max.paths < 1){
    stop("'max.paths' cannot be smaller than 1")
  }

  ##diagnostics: determine ar and arch lags:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])

  ##user-defined diagnostics?:
  if( is.null(user.diagnostics) ){
    user.diagnostics <- object$call$user.diagnostics  
  }
  
  ##if( user-defined estimator ):
  if( !is.null(object$aux$user.estimator) ){
    user.estimator <- object$call$user.estimator
    if( is.null(plot) || identical(plot,TRUE) ){
      plot <- FALSE
      message("  user-defined estimator: 'plot' set to FALSE")
    }
  } #close if user estimator

  ##if( default estimator ):
  if( is.null(object$call$user.estimator) ){

    ##determine ols method:
    if( is.null(vcov.type) ){ vcov.type <- object$aux$vcov.type }
    vcovTypes <- c("a", "b", "ordinary", "white", "newey-west")
    olsMethod <- charmatch(vcov.type, vcovTypes)
    if( (olsMethod%in%c(3,4,5))==FALSE ){ stop("'vcov.type' invalid") }
    
    ##ols arguments:
    user.estimator <- list()
    user.estimator$name <- "ols"
    user.estimator$tol <- object$aux$tol 
    user.estimator$LAPACK <- object$aux$LAPACK
    user.estimator$method <- olsMethod
    
    ##variance specification:
    if( is.null(object$variance.result) ){
      user.estimator$variance.spec <- NULL
    }else{
      user.estimator$variance.spec <- list(vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg)
    }
    
  } #close if( default estimator )

  ##gof arguments:
  if( is.null(gof.function) ){

    ##determine info method:
    infoTypes <- c("sc","aic","aicc","hq")
    whichMethod <- charmatch(info.method[1], infoTypes)
    info.method <- infoTypes[ whichMethod ]
    
    ##make gof arguments:
    gof.function <- list(name="infocrit", method=info.method)
    gof.method <- "min"
    
  }
  
  ##------------------
  ## 2 gets modelling
  ##------------------

  ##out list:
  out <- list()
  out$time.started <- date()
  out$time.finished <- NA ##added below, towards the end
  out$call <- sys.call() #used by coef.arx

  ##add gum results and diagnostics to out:
  tmp <- matrix(0, NROW(object$mean.results), 2)
  colnames(tmp) <- c("reg.no.", "keep")
  tmp[,1] <- 1:NROW(tmp) #fill reg.no. column
  tmp[keep,2] <- 1 #fill keep column
  out$gum.mean <- cbind(tmp, object$mean.results)
  out$gum.variance <- object$variance.results
  out$gum.diagnostics <- object$diagnostics

  ##do the gets:
  est <- getsFun(object$aux$y, object$aux$mX,
    user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
    wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
    arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
    user.diagnostics=user.diagnostics, gof.function=gof.function,
    gof.method=gof.method, keep=keep, include.gum=include.gum,
    include.1cut=include.1cut, include.empty=include.empty,
    max.paths=max.paths, turbo=turbo, tol=tol, max.regs=NULL,
    print.searchinfo=print.searchinfo, alarm=alarm)
  est$time.started <- NULL
  est$time.finished <- NULL
  out$time.finished <- date()
  est$call <- NULL
  out <- c(out, est)

  ##---------------------
  ## 3 estimate specific
  ##---------------------

  ## if no search has been undertaken:
  if( is.null(out$terminals.results) ){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }

  ##if search has been undertaken:
  if( !is.null(out$terminals.results) ){

    if( length(out$specific.spec)>0 ){
      out$specific.spec <- sort(out$specific.spec)
    }

    ##prepare estimation:
    yadj <- zoo(object$aux$y, order.by=object$aux$y.index)
    if( length(out$specific.spec)==0 ){
      mXadj <- NULL
    }else{
      mXadj <- cbind(object$aux$mX[, out$specific.spec ])
      colnames(mXadj) <- object$aux$mXnames[ out$specific.spec ]
      mXadj <- zoo(mXadj, order.by=object$aux$y.index)
    }
    if(is.null(ar.LjungB)){ ar.LjungB <- object$aux$qstat.options[1] }
    if(is.null(arch.LjungB)){ arch.LjungB <- object$aux$qstat.options[2] }
    if( is.null(normality.JarqueB) ){
      normality.JarqueB <- FALSE
    }else{
      normality.JarqueB <- TRUE
    }
        
    ##if( default estimator ):
    if( is.null(object$call$user.estimator) ){
      ##estimate specific model:
      est <- arx(yadj, mxreg=mXadj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
        zero.adj=object$aux$zero.adj,
        vc.adj=object$aux$vc.adj, vcov.type=vcov.type,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, plot=FALSE)
    } #end if( default estimator )

    ##if( user-defined estimator ):
    if( !is.null(object$call$user.estimator) ){
      ##estimate specific:
      est <- arx(yadj, mxreg=mXadj, user.estimator=user.estimator,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, plot=FALSE)
    } #end if( user-defined estimator )

    ##delete, rename, add:
    est$call <- est$date <- NULL
    where.diagnostics <- which(names(est)=="diagnostics")
    if(length(where.diagnostics)>0){
      names(est)[where.diagnostics] <- "specific.diagnostics"
    }
    est$aux$y.name <- object$aux$y.name
    est$aux$call.gum <- object$call #used by predict.gets()
    est <- unclass(est)
    out <- c(out,est)

  } #end if( !is.null(out$terminals.results) )

  ##------------------
  ## 4 output
  ##------------------

  ##finalise and return:
  out <- c(list(date=date(), gets.type="getsm"), out)
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
  info.method=c("sc", "aic", "aicc", "hq"),
  gof.function=NULL, gof.method=NULL, keep=c(1),
  include.gum=FALSE, include.1cut=TRUE, include.empty=FALSE,
  max.paths=NULL, tol=1e-07, turbo=FALSE, print.searchinfo=TRUE,
  plot=NULL, alarm=FALSE)
{
  ### ARGUMENTS ###########

  ##obligatory:
  vc=TRUE
  vcov.type <- "ordinary"

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
  #if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }

  ##gof arguments:
  if( is.null(gof.function) ){
    ##determine info method:
    infoTypes <- c("sc","aic","aicc","hq")
    whichMethod <- charmatch(info.method[1], infoTypes)
    info.method <- infoTypes[ whichMethod ]
    ##make gof arguments:
    gof.function <- list(name="infocrit", method=info.method)
    gof.method <- "min"
  }


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
    warning("Regressor 1 included into 'keep'")
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
  est <- getsFun(loge2, mX,
    user.estimator=list(name="ols", untransformed.residuals=eadj,
    tol=object$aux$tol, LAPACK=object$aux$LAPACK, method=6),
    gum.result=NULL, t.pval=t.pval, wald.pval=wald.pval, do.pet=do.pet,
    ar.LjungB=ar.LjungB, arch.LjungB=arch.LjungB,
    normality.JarqueB=normality.JarqueB, user.diagnostics=user.diagnostics,
    gof.function=gof.function, gof.method=gof.method, keep=keep,
    include.gum=include.gum, include.1cut=include.1cut,
    include.empty=include.empty, max.paths=max.paths, turbo=turbo,
    tol=tol, max.regs=NULL, print.searchinfo=print.searchinfo,
    alarm=alarm)
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
  if( is.null(ar.LjungB) ){ ar.LjungB <- object$aux$qstat.options[1] }
  if( is.null(arch.LjungB) ){ arch.LjungB <- object$aux$qstat.options[2] }
  if( is.null(normality.JarqueB) ){
    normality.JarqueB <- FALSE
  }else{
    normality.JarqueB <- TRUE
  }

  ## estimate model:
  est <- arx(e, vc=TRUE, vxreg=vXadj,
    zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
    qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
    normality.JarqueB=normality.JarqueB,
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
  est$vcov.mean <- NULL
  est$aux$vxreg <- est$aux$vxreg.index <- NULL
  est$aux$y.name <- "e"

  ## finalise:
  est <- unclass(est)
  out <- c(out,est)

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
  if( is.null(spec) ){
    spec <- switch(object$gets.type, getsm="mean", getsv="variance")
  }
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
    stop("object not of class 'gets' or 'isat'")
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
  n.sim=5000, innov=NULL, probs=NULL, ci.levels=NULL, 
  quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
  plot.options=list(), ...)  
{

  ##create new object to add stuff to in order to use predict.arx()
  objectNew <- object


  ##-----------------------------------
  ## arguments mean-equation:
  ##-----------------------------------

  ##coefficients of mean spec in final model:
  coefsMean <- coef.arx(objectNew, spec="mean")

  ##there is no mean equation:
  if( length(coefsMean)==0 ){

    objectNew$call$mc <- NULL
    objectNew$call$ar <- NULL
    objectNew$call$ewma <- NULL
    objectNew$call$mxreg <- NULL

  }

  ##there is a mean equation:
  if( length(coefsMean)>0 ){

    ##initiate index counter (used for mxreg):
    indxCounter <- 0

    ##mc argument:
    mconstRetained <- "mconst" %in% names(coefsMean)
    if( mconstRetained ){
      objectNew$call$mc <- TRUE
      indxCounter <- indxCounter + 1
    }else{
      objectNew$call$mc <- NULL
    }
    
    ##ar argument:
    gumTerms <- eval(object$aux$call.gum$ar)
    gumNamesAr <- paste0("ar", gumTerms)
    whichRetained <- which( gumNamesAr %in% names(coefsMean) )
    if( length(whichRetained)==0 ){
      objectNew$call$ar <- NULL
    }else{
      objectNew$call$ar <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
        
    ##ewma argument:
    gumTerms <- eval(object$aux$call.gum$ewma)
    gumNamesEwma <- paste0("EqWMA(", gumTerms$length, ")")
    whichRetained <- which( gumNamesEwma %in% names(coefsMean) )
    if( length(whichRetained)==0 ){
      objectNew$call$ewma <- NULL
    }else{
      objectNew$call$ewma <-
        list( length=gumTerms$length[ whichRetained ] )
      indxCounter <- indxCounter + length(whichRetained)
    }

    ##mxreg argument:
    if(indxCounter==0){ whichRetainedCoefs <- coefsMean }
    if(indxCounter>0){ whichRetainedCoefs <- coefsMean[ -c(1:indxCounter) ] }
    if( length(whichRetainedCoefs)==0 ){
      objectNew$call$mxreg <- NULL
    }else{
      whichRetainedNames <- names(whichRetainedCoefs)
      objectNew$call$mxreg <- whichRetainedNames
#more correct (but not needed, since mxreg only needs to be non-NULL)?:
#      whichRetained <- which( object$aux$mXnames %in% whichRetainedNames )
#      mxreg <- cbind(object$aux$mX[, whichRetained ])
#      colnames(mxreg) <- whichRetainedNames
#      objectNew$call$mxreg <- xreg
    }

  } #end if( length(coefsMean)>0 )
  

  ##-----------------------------------
  ## arguments variance-equation:
  ##-----------------------------------

  ##coefficients of variance spec in final model:
  coefsVar <- coef.arx(objectNew, spec="variance")
  if( length(coefsVar)>0 ){ #remove Elnz2 estimate:
    coefsVar <- coefsVar[ -length(coefsVar) ]  
  }

  ##there is no variance equation:
  if( length(coefsVar)==0 ){

    objectNew$call$vc <- NULL
    objectNew$call$arch <- NULL
    objectNew$call$asym <- NULL
    objectNew$call$log.ewma <- NULL
    objectNew$call$vxreg <- NULL

  }

  ##there is a variance equation:
  if( length(coefsVar)>0 ){

    ##vc argument (always present in variance equations):
    objectNew$call$vc <- TRUE
    indxCounter <- 1 #used for vxreg
    
    ##arch argument:
    gumTerms <- eval(object$aux$call.gum$arch)
    gumNamesArch <- paste0("arch", gumTerms)
    whichRetained <- which( gumNamesArch %in% names(coefsVar) )
    if( length(whichRetained)==0 ){
      objectNew$call$arch <- NULL
    }else{
      objectNew$call$arch <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
    
    ##asym argument:
    gumTerms <- eval(object$aux$call.gum$asym)
    gumNamesAsym <- paste0("asym", gumTerms)
    whichRetained <- which( gumNamesAsym %in% names(coefsVar) )
    if( length(whichRetained)==0 ){
      objectNew$call$asym <- NULL
    }else{
      objectNew$call$asym <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
    
    ##log.ewma argument:
    gumTerms <- eval(object$aux$call.gum$log.ewma)
    gumNamesLogEwma <- paste0("logEqWMA(", gumTerms$length, ")")
    whichRetained <- which( gumNamesLogEwma %in% names(coefsVar) )
    if( length(whichRetained)==0 ){
      objectNew$call$log.ewma <- NULL
    }else{
      objectNew$call$log.ewma <-
        list( length=gumTerms$length[ whichRetained ] )
      indxCounter <- indxCounter + length(whichRetained)
    }

    ##vxreg argument:
    whichRetainedCoefs <- coefsVar[ -c(1:indxCounter) ]
    if( length(whichRetainedCoefs)==0 ){
      objectNew$call$vxreg <- NULL
    }else{
      whichRetainedNames <- names(whichRetainedCoefs)
      objectNew$call$vxreg <- whichRetainedNames
#more correct (but not needed, since vxreg only needs to be non-NULL)?:
#      whichRetained <- which( object$aux$vXnames %in% whichRetainedNames )
#      vxreg <- cbind(object$aux$vX[, whichRetained ])
#      colnames(vxreg) <- whichRetainedNames
#      objectNew$call$vxreg <- vxreg
    }

  } #end if( length(coefsVar)>0 )


  ##----------------------------------
  ## pass arguments on to predict.arx:
  ##----------------------------------

  result <- predict.arx(objectNew, spec=spec, n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=newvxreg, newindex=newindex,
    n.sim=n.sim, innov=innov, probs=probs, ci.levels=ci.levels,
    quantile.type=quantile.type, return=return, verbose=verbose,
    plot=plot, plot.options=plot.options)

  ##-------------------
  ## return forecasts:
  ##-------------------

  if(return){ return(result) }

} #close predict.gets

##==================================================
## print gets results
print.gets <- function(x, signif.stars=TRUE, ...)
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
  estType <- ifelse(is.null(x$aux$user.estimator),
      "Ordinary Least Squares (OLS)", "User defined")
  cat("Method:", estType, "\n")

  ##header - if mean:
  if( specType=="mean" ){
    vcovType <- "Unknown"
    if( !is.null(x$aux$vcov.type) ){
      vcovType <- switch(x$aux$vcov.type,
        ordinary = "Ordinary", white = "White (1980)",
        "newey-west" = "Newey and West (1987)")
    }
    cat("Variance-Covariance:", vcovType, "\n")
    if(!is.null(x$aux$y.n)){
      cat("No. of observations (mean eq.):", x$aux$y.n, "\n") }
  }

  ##header - if variance:
  if( specType=="variance" ){
    if(!is.null(x$aux$loge2.n)){
      cat("No. of observations (variance eq.):",
        x$aux$loge2.n, "\n") }
  }

  ##header - sample info:
  if( !is.null(x$residuals) ){
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
    printCoefmat(x$gum.mean, tst.ind=c(1,2),
      signif.stars=signif.stars)
#OLD:
#    printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2),
#      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }
  if( !is.null(x$gum.variance) ){
    cat("\n")
    cat("GUM log-variance equation:\n")
    cat("\n")
    if(specType=="mean"){
      printCoefmat(x$gum.variance, signif.stars=FALSE)
    }
    if(specType=="variance"){
      printCoefmat(x$gum.variance, tst.ind=c(1,2),
        signif.stars=signif.stars)
#OLD:
#      printCoefmat(x$gum.variance, dig.tst=0, tst.ind=c(1,2),
#        signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
    }
  }
  if( !is.null(x$gum.diagnostics) ){
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
  
  ##specific mean model:
  if( specType=="mean" && !is.null(x$terminals.results) ){
#OLD: if( specType=="mean" && !is.null(x$specific.spec) ){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if( !is.null(x$mean.results) ){
      printCoefmat(x$mean.results, signif.stars=signif.stars)
#      printCoefmat(x$mean.results, signif.stars=FALSE)
    }
    if( length(x$specific.spec)==0 ){
#OLD: if( x$specific.spec[1]==0 ){
      cat("  the empty model\n")
    }
  }

  ##specific log-variance model:
  if( !is.null(x$variance.results) ){
    cat("\n")
    cat("SPECIFIC log-variance equation:\n")
    cat("\n")
    printCoefmat(x$variance.results, signif.stars=signif.stars)
#    printCoefmat(x$variance.results, signif.stars=FALSE)
#    printCoefmat(x$variance.results, dig.tst=0, tst.ind=c(1,2),
#      signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }

  ##diagnostics and fit:
  if( !is.null(x$specific.diagnostics) ){

    #fit-measures:
    mGOFnames <- "SE of regression"
    mGOF <- sigma.gets(x) 
    if( specType == "mean" ){
      mGOFnames <- c(mGOFnames, "R-squared")
      mGOF <- rbind(mGOF, rsquared(x))
    }
    logl <- logLik.arx(x)
    mGOFnames <- c(mGOFnames,
      paste0("Log-lik.(n=", attr(logl,"n"), ")") )
    mGOF <- rbind(mGOF, as.numeric(logl))
    rownames(mGOF) <- mGOFnames
    colnames(mGOF) <- ""

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
  summary.default(object)
} #end summary.gets

##==================================================
## extract terminal models
terminals <- function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$terminals)
  }else{
    stop("object not of class 'gets' or 'isat'")
  }
} #end terminals

##==================================================
## LaTeX code (equation form)
toLatex.gets <- function(object, ...)
{
  printtex(object, ...)
} #end toLatex.gets

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
    message("STATA code to estimate the model:\n")
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
      message("STATA code to estimate the model:\n")
      message(" ", out$regress, "\n")
    }
  } #end if(!is.null(file))

  ##out:
  if(return){ return(out) }

} #close stata()

##==================================================
##generate latex-code (equation form):
printtex <- function(x, fitted.name=NULL, xreg.names=NULL,
  digits=4, intercept=TRUE, gof=TRUE, diagnostics=TRUE, nonumber=FALSE,
  nobs="T")
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

  ##coef names:
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
  txtAddEq <- ifelse(gof+diagnostics>0, " \\\\[2mm]", "")
  eqtxt <- paste0("  ", yName, " &=& ", eqtxt, "", txtAddNonumber,
    txtAddEq, " \n")

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
        " \\qquad ", nobs, " = ", iT, " \\nonumber \\\\ \n", sep="")
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

##==================================================
##convert to model of class 'lm':
as.lm <- function(object)
{

  ##what kind of class?:
  objectClass <- class(object)
  classOK <-
    ifelse( objectClass %in% c("arx","gets","isat"), TRUE, FALSE)

  ##class not OK:
  if(!classOK){
    stop("'object' must be of class 'arx', 'gets' or 'isat'")
  }

  ##class OK:
  if(classOK){
    y <- object$aux$y
    x <- object$aux$mX
    colnames(x) <- object$aux$mXnames
    yx <- data.frame(y, x)
    result <- lm(formula = y ~ . - 1, data = yx)
  }
  
  ##return result:
  return(result)
  
}
