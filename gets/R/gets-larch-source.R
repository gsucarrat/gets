###########################################################
## This file contains the source of the larch() functions.
##
## First created 20 May 2024, Madrid.
##
## Contents:
##
## TO DO: larchSim()
## larchEstfun()
## larch()
## coef.larch()         #extraction functions
## fitted.larch()       #(all are S3 methods)
## gets.larch()
## logLik.larch()
## model.matrix.larch()
## nobs.larch()
## plot.larch()
## predict.larch()
## print.larch()
## residuals.larch()
## summary.larch()
## toLatex.larch()
## vcov.larch()
##
###########################################################

###########################################################
## The function larchSim() simulates from a heterogenous
## log-ARCH-X model
###########################################################

## TO DO: larchSim <- function(n, ...)


###########################################################
## The function larchEstfun() estimates a log-ARCH-X model
###########################################################

larchEstfun <- function(loge2, x, e, vcov.type=c("robust", "hac"),
  tol=1e-07)
{
  ##initialise:
  result <- list()
  result$n <- length(loge2)
  result$k <- ncol(x)
  result$df <- result$n - result$k
  
  ##step 1: ols estimation and fit
  out <- ols(loge2, x, tol=tol, method=2)

  ##step 2: log-arch estimation, fit, residuals and log-likelihood
  uhat <- out$residuals
  e2 <- e^2
  zstar2 <- e2/exp(out$fit)
  Ezstar2 <- mean( zstar2 )
  lnEzstar2 <- log(Ezstar2) #\widehat{\tau}
  result$coefficients <- out$coefficients
  result$coefficients[1] <- out$coefficients[1] + lnEzstar2
  result$fit <- exp(out$fit) * Ezstar2
  result$residuals <- e/sqrt(result$fit)
  logl <- -result$n*log(2*pi)/2 - sum(log(result$fit))/2 - sum(e2/result$fit)/2

  ##step 3: compute vcov
  if( !is.null(vcov.type) ){

    ##vcov.type argument:
    types <- c("robust", "hac")
    whichType <- charmatch(vcov.type[1], types)
    vcov.type <- types[ whichType ]

    ##build G matrix
    Extx <- crossprod(x)/result$n
    Extxinv <- result$n*out$xtxinv
    Ezstar2tx <- rbind( colMeans( zstar2 * x )/Ezstar2 )
    mG <- cbind( Extxinv, rep(0, result$k))
    mG <- rbind( mG, cbind(-Ezstar2tx %*% Extxinv, 1/Ezstar2) ) 

    ##build Sigma_w matrix
    w <- cbind(cbind(uhat*x), zstar2 - Ezstar2) #the w's
    Ewtw <- crossprod(w)/result$n
    if( vcov.type=="robust" ){
      mSigma_w <- Ewtw
    }
    if( vcov.type=="hac" ){
      iL <- round(result$n^(1/4), digits=0) ##lag truncation (bandwidth)
      vWeight <- 1 - 1:iL/(iL+1) ##kernel weights
      vWeightSqrt <- sqrt(vWeight)
      mS0 <- Ewtw    
      mSum <- 0
      for(l in 1:iL){
        mSadjw <- w*vWeightSqrt[l]
        mSadjwNo1 <- mSadjw[-c(1:l),]
        mSadjwNo2 <- mSadjw[-c(c(result$n-l+1):result$n),]
        mSum <- mSum +
          crossprod(mSadjwNo1, mSadjwNo2)/c(result$n-l) +
          crossprod(mSadjwNo2, mSadjwNo1)/c(result$n-l)
      }
      mSigma_w <- mS0 + mSum
    }  

    ##finalise vcov:
    mSigma_phi <- mG %*% mSigma_w %*% t(mG)
    mA <- matrix(0,result$k,result$k)
    diag(mA) <- 1
    mA <- cbind(mA, rep(0, NROW(mA)))
    mA[1,NCOL(mA)] <- 1
    vcovasym <- mA %*% mSigma_phi %*% t(mA) #asymptotic vcov
    result$vcov <- vcovasym/result$n #finite sample vcov
    result$vcov.type <- vcov.type

  } #close if( !is.null(vcov.type) )
    
  ##complete and return result:
  result$logl <- logl
  return(result)
  
} #close larchEstfun()


###########################################################
## The function larch() estimates a log-ARCH-X model
###########################################################

larch <- function(e, vc=TRUE, arch=NULL, harch=NULL, asym=NULL,
  asymind=NULL, log.ewma=NULL, vxreg=NULL, zero.adj=NULL,
  vcov.type=c("robust", "hac"), #bandwidth=NULL
  qstat.options=NULL, normality.JarqueB=FALSE, #user.estimator=NULL, user.diagnostics=NULL,
  tol=1e-07, singular.ok=TRUE, plot=NULL)
{
  ## contents:
  ##
  ## 1 initiate
  ## 2 create aux list
  ## 3 estimation
  ## 4 prepare result
  ## 5 finalise and return result


  ##-----------------------------------
  ## 1 initiate
  ##-----------------------------------

  ##record call:
  sysCall <- sys.call()

  ##check argument(s):
  if( !vc==TRUE ){ stop("argument 'vc' must be TRUE") }
  
  ##regressand, regressors:
  e.name <- deparse(substitute(e))
  tmp <- regressorsVariance(e, vc=TRUE, arch=arch, harch=harch, asym=asym, 
    asymind=asymind, log.ewma=log.ewma, vxreg=vxreg, prefix="v", 
    zero.adj=zero.adj, return.regressand=TRUE, return.as.zoo=TRUE,
    na.trim=TRUE, na.omit=FALSE)
  e <- as.zoo(e)
  e <- zoo( as.vector(coredata(e)), order.by=index(e)) #convert to vector (recall: regressorsVariance() does not accept NCOL(e) > 1)
  e <- window(e, start=index(tmp)[1], end=index(tmp)[NROW(tmp)]) 

  ##singularity ok?:
  if( singular.ok && NCOL(tmp)>2 ){
    tmpx <-
      colnames(dropvar(tmp[,-1], tol=tol, LAPACK=FALSE, silent=TRUE))
    droppedXs <- which( (colnames(tmp[,-1]) %in% tmpx)==FALSE )
    if( length(droppedXs)>0 ){
      droppedXsNames <- colnames(tmp)[c(1+droppedXs)]
      tmp <- tmp[,-c(1+droppedXs)]
      warning(
        "Regressor(s) removed due to singularity:\n",
        paste(" ", droppedXsNames)
      )            
    } #end if( length(droppedXs)>0 )
  } #end if( singular.ok )

  ##determine vcov type:
  types <- c("robust", "hac")
  whichType <- charmatch(vcov.type[1], types)
  vcov.type <- types[ whichType ]

  ##determine qstat.options:
  if( is.null(qstat.options) ){
    ar.lag <- 1
    if( is.null(arch) ){ arch.lag <- 1 }else{ arch.lag <- max(arch)+1 }
    qstat.options <- c(ar.lag, arch.lag)
  }

  
  ##-----------------------------------
  ## 2 create aux list
  ##-----------------------------------

  ##aux: auxiliary list, also to be used by getsv
  aux <- list()
#  aux$e.original <- e.original
  aux$e <- coredata(e)
  aux$e2 <- aux$e^2
  aux$e.index <- index(e)
  aux$e.name <- e.name #recorded above, in the beginning
  aux$loge2 <- coredata(tmp[,1])
  aux$vX <- cbind(coredata(tmp[,-1]))
  aux$vXnames <- colnames(tmp)[-1]
  colnames(aux$vX) <- NULL
#  aux$vXncol <- NCOL(aux$vX)
#  aux$vc <- TRUE #obligatory, but may change in the future
  aux$zero.adj <- zero.adj
#do we really need these?:
#  aux$arch <- arch
#  aux$harch <- harch
#  aux$asym <- asym
#  aux$asymind <- asymind
#  aux$log.ewma <- log.ewma
#  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$normality.JarqueB <- normality.JarqueB
  aux$tol <- tol


  ##-----------------------------------
  ## 3 estimation
  ##-----------------------------------

  out <- larchEstfun(aux$loge2, aux$vX, aux$e, vcov.type=vcov.type,
    tol=tol)
  names(out$coefficients) <- aux$vXnames
  colnames(out$vcov) <- aux$vXnames
  rownames(out$vcov) <- aux$vXnames
  where <- which( names(out) == "fit" )
  names(out)[ where ] <- "fitted" #rename
  aux <- c(aux,out)

  ##build 'results' data frame:
  s.e. <- sqrt(as.vector(diag(aux$vcov)))
  t.stat <- aux$coefficients/s.e.
  p.val <- pt(abs(t.stat), aux$k, lower.tail=FALSE)*2
  aux$results <-
    as.data.frame(cbind(aux$coefficients, s.e., t.stat, p.val))
  colnames(aux$results) <- c("coef", "std.error", "t-stat", "p-value")
  rownames(aux$results) <- aux$vXnames
  names(aux$coefficients) <- aux$vXnames


  ##-----------------------------------
  ## 4 prepare result
  ##-----------------------------------

  ##diagnostics:
  if( "residuals" %in% names(aux) ){
    ar.LjungBarg <- c(qstat.options[1],0)
    arch.LjungBarg <- c(qstat.options[2],0)
    if( normality.JarqueB==FALSE ){
      normality.JarqueBarg <- NULL
    }else{
      normality.JarqueBarg <- as.numeric(normality.JarqueB)
    }
  }else{
    ar.LjungBarg <- NULL
    arch.LjungBarg <- NULL
    normality.JarqueBarg <- NULL
  }
  aux$diagnostics <- diagnostics(aux,
    ar.LjungB=ar.LjungBarg, arch.LjungB=arch.LjungBarg,
    normality.JarqueB=normality.JarqueBarg,
    user.fun=NULL, verbose=TRUE)
    
  ##add zoo-indices:
  if( "fitted" %in% names(aux) ){
    aux$fitted <- zoo(aux$fitted, order.by=aux$e.index)
  }
  if( "residuals" %in% names(aux) ){
    aux$residuals <- zoo(aux$residuals, order.by=aux$e.index)
  }


  ##-----------------------------------
  ## 5 finalise and return result
  ##-----------------------------------

  versionTxt <- paste0("gets ", packageVersion("gets"), " under ",
    version$version.string)
  aux <-
    c(list(call=sysCall, date=date(), version=versionTxt), aux)
  class(aux) <- "larch"

  ##plot:
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if( plot ){ plot.larch(aux) }

  ##return result:
  return(aux)

} #close larch() function


##########################################################
## coef.larch()
###########################################################

coef.larch <- function(object, ...){ return(object$coefficients) }


###########################################################
## fitted.larch()
###########################################################

fitted.larch <- function(object, ...){ return(object$fitted) }


###########################################################
## gets.larch()
###########################################################

gets.larch <- function(x, t.pval=0.05, wald.pval=t.pval, do.pet=TRUE,
  ar.LjungB=NULL, arch.LjungB=NULL, normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc", "aic", "aicc", "hq"),
  gof.function=NULL, gof.method=NULL, keep=c(1), include.gum=FALSE,
  include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
  turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE, ...)
{
  ## contents:
  ## 1 arguments
  ## 2 gets modelling
  ## 3 estimate specific
  ## 4 result
  
  ##------------------
  ## 1 arguments
  ##------------------

  ##diagnostics: ar argument
  if( !is.null(ar.LjungB) && is.vector(ar.LjungB, mode="double") ){
      ar.LjungB <- list(lag=ar.LjungB[1], pval=ar.LjungB[2])
  }
  if( !is.null(ar.LjungB) && is.null(ar.LjungB$lag) ){
    ar.LjungB$lag <- x$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  ##(NULL if ar.LjungB is NULL)

  ##diagnostics: arch argument
  if( !is.null(arch.LjungB) && is.vector(arch.LjungB, mode="double") ){
      arch.LjungB <- list(lag=arch.LjungB[1], pval=arch.LjungB[2])
  }
  if( !is.null(arch.LjungB) && is.null(arch.LjungB$lag) ){
    arch.LjungB$lag <- x$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
  ##(NULL if arch.LjungB is NULL)

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
  out$time.finished <- NA
  out$gum.call <- x$call
  mX <- x$vX
  colnames(mX) <- x$vXnames
  if( !(1 %in% keep) ){
    keep <- union(1,keep)
    warning("Regressor 1 included into 'keep'")
  }

  ##add gum results and diagnostics to out:
  tmp <- matrix(0, NROW(x$results), 2)
  colnames(tmp) <- c("reg.no.", "keep")
  tmp[,1] <- 1:NROW(tmp) #fill reg.no. column
  tmp[keep,2] <- 1 #fill keep column
  out$gum.results <- cbind(tmp, x$results)
  out$gum.diagnostics <- x$diagnostics

  ##print start model (gum) info:
  if( print.searchinfo ){
    if( !is.null(out$gum.results) ){
      cat("\n")
      cat("GUM log-variance equation:\n")
      cat("\n")
      printCoefmat(out$gum.results, cs.ind=c(3,4), tst.ind=c(5),
        signif.stars=TRUE, P.values=TRUE)
    }
    if( !is.null(out$gum.diagnostics) ){
      cat("\n")
      cat("Diagnostics:\n")
      cat("\n")
      printCoefmat(out$gum.diagnostics, tst.ind=2, signif.stars=TRUE)
      cat("\n")
    }
  } #end if( print.searchinfo )

  ##do the gets:
  est <- getsFun(x$loge2, mX, user.estimator=list(name="larchEstfun",
    e=x$e, vcov.type=x$vcov.type, tol=x$tol),
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

  ##print paths, terminals and retained regressors:
  if( print.searchinfo && !is.null(est$terminals.results) ){

    ##paths:
    if( length(est$paths)>0 ){
      cat("\n")
      for(i in 1:length(est$paths)){
        txt <- paste0(est$paths[[i]], collapse=" ")
        txt <- paste0("  Path ", i, ": ", txt)    
        cat(txt, "\n")
      }
    }

    ##print terminals:
    cat("\n")
    cat("Terminal models:\n")
    cat("\n")
    print(est$terminals.results)    

    ##retained regressors:
    cat("\n")
    cat("Retained regressors (final model):\n")
    if( length(est$specific.spec)==0 ){
      cat("  none\n")
    }else{
      cat(paste0("  ", x$vXnames[as.numeric(est$specific.spec)]), "\n")
    }
    
  } #end   if( print.searchinfo && !is.null(est$terminals.results) )

  ##messages:
  if( print.searchinfo && !is.null(est$messages) ){
    message("\n", appendLF=FALSE)
    message("Messages:", appendLF=TRUE)
    message("\n", appendLF=FALSE)
    message(est$messages)
  }

  ##---------------------
  ## 3 estimate specific
  ##---------------------

  ##if there are no terminals, then use the gum as specific:
  if( is.null(est$terminals.results) ){
    out$messages <- paste0(out$messages,
      "- No terminal models, so the final model equals the GUM")
    out <- c(out,x)
  }  
  
  ##if there are terminals, then estimate specific:
  if( !is.null(est$terminals.results) ){

    ## prepare estimation:
    e <- zoo(cbind(x$e), order.by=x$e.index)
    colnames(e) <- "e"
    specificadj <- setdiff(out$specific.spec, 1)
    if(length(specificadj)==0){
      vXadj <- NULL
    }else{
      vXadj <- cbind(x$vX[,specificadj])
      colnames(vXadj) <- x$vXnames[specificadj]
      vXadj <- zoo(vXadj, order.by=x$e.index)
    }
    if( is.null(ar.LjungB) ){ ar.LjungB <- x$qstat.options[1] }
    if( is.null(arch.LjungB) ){ arch.LjungB <- x$qstat.options[2] }
    if( is.null(normality.JarqueB) ){
      normality.JarqueB <- FALSE
    }else{
      normality.JarqueB <- TRUE
    }
  
    ## estimate specific model:
    est <- larch(e, vc=TRUE, vxreg=vXadj, zero.adj=x$zero.adj,
      vcov.type=x$vcov.type, qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
      normality.JarqueB=normality.JarqueB, tol=x$tol, singular.ok=FALSE,
      plot=FALSE)
      
    ## delete, rename and change various stuff:
    est$call <- NULL
    est$date <- NULL
    est$e.name <- x$e.name
    est <- unclass(est)
  
    ##build new call (needed for predict.larch()):
    newCall <- x$call
    coefNames <- names(coef(est))
    indxCounter <- 1 #needed for vxreg
      
    ##arch argument for new call:
    gumTerms <- eval(x$call$arch)
    gumNames <- paste0("arch", gumTerms)
    whichRetained <- which( gumNames %in% coefNames )
    if( length(whichRetained)==0 ){
      newCall$arch <- NULL
    }else{
      newCall$arch <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
  
    ##harch argument for new call:
    gumTerms <- eval(x$call$harch)
    gumNames <- paste0("harch", gumTerms)
    whichRetained <- which( gumNames %in% coefNames )
    if( length(whichRetained)==0 ){
      newCall$harch <- NULL
    }else{
      newCall$harch <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
  
    ##asym argument for new call:
    gumTerms <- eval(x$call$asym)
    gumNames <- paste0("asym", gumTerms)
    whichRetained <- which( gumNames %in% coefNames )
    if( length(whichRetained)==0 ){
      newCall$asym <- NULL
    }else{
      newCall$asym <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
  
    ##asymind argument for new call:
    gumTerms <- eval(x$call$asymind)
    gumNames <- paste0("asymind", gumTerms)
    whichRetained <- which( gumNames %in% coefNames )
    if( length(whichRetained)==0 ){
      newCall$asymind <- NULL
    }else{
      newCall$asymind <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
  
    ##log.ewma argument for new call:
    gumTerms <- eval(x$call$log.ewma)
    gumNames <- paste0("logEqWMA(", gumTerms, ")")
    whichRetained <- which( gumNames %in% coefNames )
    if( length(whichRetained)==0 ){
      newCall$log.ewma <- NULL
    }else{
      newCall$log.ewma <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }
  
    ##vxreg argument for new call:
    whichRetainedCoefs <- coefNames[ -c(1:indxCounter) ]
    if( length(whichRetainedCoefs)==0 ){
      newCall$vxreg <- NULL
    }else{
      whichRetained <- which( est$vXnames %in% whichRetainedCoefs )
      newCall$vxreg <- est$vXnames[whichRetained] 
  #more correct, but predict.larch() only needs that vxreg is not NULL:
  #    vxreg <- cbind( est$vX[, whichRetained ] )
  #    colnames(vxreg) <- est$vXnames[whichRetained] 
  #    newCall$vxreg <- vxreg
    }              
  
    ##add new call to 'est', merge with 'out':
    est <- c(list(call=newCall), est)
    out <- c(out,est)

  } #close if( !is.null(est$terminals.results) )


  ##------------------
  ## 4 result
  ##------------------

  out$time.finished <- date()
  out <- c(list(date=out$time.finished), out)
  class(out) <- "larch"
  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.larch(out) }
  return(out)
  
} #close gets.larch() function


###########################################################
## logLik.larch()
###########################################################

logLik.larch <- function(object, ...){ return(object$logl) }


###########################################################
## model.matrix.larch()
###########################################################

model.matrix.larch <- function(object, response=FALSE,
  as.zoo=TRUE, ...)
{
  result <- NULL
  result <- object$vX
  colnames(result) <- object$vXnames
  if( response==TRUE ){
    loge2 <- cbind(object$loge2)
    colnames(loge2) <- "loge2"
    result <- cbind(loge2,result)
  }
  if( as.zoo==TRUE ){ result <- zoo(result, order.by=object$e.index) }
  return(result)
} #close model.matrix.larch()


###########################################################
## nobs.larch()
###########################################################

nobs.larch <- function(object, ...){ return(object$n) }


###########################################################
## print.larch() prints the estimation result
###########################################################

plot.larch <- function(x, col=c("red","blue"), lty=c("solid","solid"),
  lwd=c(1,1), ...)
{
  ##arguments:
  ##----------
  
  ##col:
  if( length(col)==1 ){
    col=rep(col,2)
  }else if (length(lty)>2){
    print("'col' needs two arguments only, but more provided. First two used.")
    col=col[1:2]
  }

  ##lty:
  if( length(lty)==1 ){
    lty=rep(lty,2)
  }else if (length(lty)>2){
    print("'lty' needs two arguments only, but more provided. First two used.")
    lty=lty[1:2]
  }

  ##lwd:
  if( length(lwd)==1 ){
    lwd=rep(lwd,2)
  }else if (length(lwd)>2){
    print("'lwd' needs two arguments only, but more provided. First two used.")
    lwd=lwd[1:2]
  }

  ##do the plotting:
  ##----------------

  vfitted <- fitted.larch(x)
  vactual <- zoo(x$e^2, order.by=x$e.index)
  actual.name <- x$e.name
  residsStd <- residuals.larch(x)

  ##get current par-values:
  def.par <- par(no.readonly=TRUE)

  ##set new par values for plot:
  par(mfrow=c(2,1))

  #set the plot margins:
  par(mar=c(2,2,0.5,0.5))

  ##plot actual values (e^2):
  minValue <- min(vactual, vfitted, na.rm=TRUE)
  maxValue <- 1.1*max(vactual, vfitted, na.rm=TRUE) #1.1: adds 10% on top
  if( is.regular(vactual) ) {
    plot(vactual, main = "", ylim=range(minValue, maxValue),
      type="l", ylab="", xlab="",col=col[2])
  }else{
    plot(as.Date(index(vactual)), coredata(vactual), main = "",
       ylim=range(minValue,maxValue), type="l", ylab="", xlab="", col=col[2])
  }

  ##plot fitted values:
  if( is.regular(vfitted) ) {
    lines(vfitted, col=col[1])
  } else {
    lines(as.Date(index(vfitted)), coredata(vfitted), col=col[1])
  }
  legend("topleft", lty=lty, lwd=lwd, ncol=2, col=col[c(2,1)],
    legend=c("actual squared","fitted variance"), bty="n")

  ##plot standardised residuals:
  minValue <- min(residsStd, na.rm=TRUE)
  maxValue <- 1.1*max(residsStd, na.rm=TRUE) #1.1: adds 10% on top
  if( is.regular(residsStd) ) {
    plot(residsStd, ylim=range(minValue,maxValue), type="h", col=col[1])
  }else{
    plot(as.Date(index(residsStd)), coredata(residsStd), 
      ylim=range(minValue,maxValue), type="h", col=col[1])
  }
  abline(0,0)
  legend("topleft", lty=1, col=col[1],
    legend=c("standardised residuals"), bty="n")

  #return to old par-values:
  par(def.par)

} #close plot.larhx()


###########################################################
## predict.larch() generate predictions
###########################################################

predict.larch <- function(object, n.ahead=12, newvxreg=NULL, newindex=NULL,
  n.sim=NULL, innov=NULL, probs=NULL, quantile.type=7, verbose=FALSE, ...)
{
  ## contents:
  ## 1 initialise
  ## 2 prepare terms
  ## 3 simulate innov
  ## 4 prepare prediction
  ## 5 generate predictions
  ## 6 prepare output
  ## 7 newindex
  ## 8 return result

  ##-----------------------
  ## 1 initialise
  ##-----------------------

  ##name of object:
  objectName <- deparse(substitute(object))

  ##check n.ahead:
  if(n.ahead < 1){ stop("n.ahead must be 1 or greater") }

  ##probs argument:
  if( !is.null(probs) ){
    if( any(probs <= 0) || any(probs >= 1) ){
      stop("the values of 'probs' must be between 0 and 1")
    }
    probs <- union(probs,probs) #ensure values are distinct/not repeated
    probs <- probs[order(probs, decreasing=FALSE)] #re-order to increasing
  }
  probsArg <- probs
  
  ##obtain zero.adj value:
  zerosWhere <- which( object$e==0 )
  if( length(zerosWhere)>0 ){
    zero.adj <- exp(object$loge2)[ zerosWhere[1] ]
  }else{
    zero.adj <- quantile(object$e2, 0.1, na.rm=TRUE) 
  }  


  ##-----------------------
  ## 2 prepare terms
  ##-----------------------

  ##record coef estimates:
  coefs <- as.numeric(coef.larch(object))
  
  ##vc:
  vconst <- as.numeric(coefs[1])
  
  ##arch:
  archMax <- 0
  archIndx <- 1
  if( "arch" %in% names(object$call) ){
    archEval <- eval(object$call$arch)
    archIndx <- 1:length(archEval) + 1
    archMax <- max(archEval)
    archCoefs <- rep(0,archMax)
    archCoefs[archEval] <- as.numeric(coefs[archIndx])
  }

  ##harch:
  harchMax <- 0
  harchIndx <- max(archIndx)
  if( "harch" %in% names(object$call) ){
    harchEval <- eval(object$call$harch)
    harchIndx <- 1:length(harchEval) + max(archIndx)
    harchMax <- max(harchEval)
    harchCoefs <- rep(0,harchMax)
    harchCoefs <- as.numeric(coefs[harchIndx])
  }

  ##asym:
  asymMax <- 0
  asymIndx <- max(harchIndx)
  if( "asym" %in% names(object$call) ){
    asymEval <- eval(object$call$asym)
    asymIndx <- 1:length(asymEval) + max(harchIndx)
    asymMax <- max(asymEval)
    asymCoefs <- rep(0,asymMax)
    asymCoefs[asymEval] <- as.numeric(coefs[asymIndx])
  }

  ##asymind:
  asymindMax <- 0  
  asymindIndx <- max(asymIndx)
  if( "asymind" %in% names(object$call) ){
    asymindEval <- eval(object$call$asymind)
    asymindIndx <- 1:length(asymindEval) + max(asymIndx)
    asymindMax <- max(asymindEval)
    asymindCoefs <- rep(0,asymindMax)
    asymindCoefs[asymindEval] <- as.numeric(coefs[asymindIndx])
  }
    
  ##log.ewma:
  logewmaMax <- 0
  logewmaIndx <- max(asymindIndx)
  if( "log.ewma" %in% names(object$call) ){
    logewmaEval <- eval(object$call$log.ewma)
    if(is.list(logewmaEval)){ logewmaEval <- logewmaEval$length }
    logewmaIndx <- 1:length(logewmaEval) + max(asymindIndx)
    logewmaMax <- max(logewmaEval)
    logewmaCoefs <- as.numeric(coefs[logewmaIndx])
  }
  
  ##backcast length:
  backcastMax <- max(archMax,harchMax,asymMax,asymindMax,logewmaMax)
  
  ##vxreg:
  vxreghat <- rep(0, n.ahead + backcastMax)
  if( !is.null(object$call$vxreg) ){

    ##check newvxreg:
    if( is.null(newvxreg) ){ stop("'newvxreg' is NULL") }
    if( NROW(newvxreg)!=n.ahead ){ stop("NROW(newvxreg) must equal n.ahead") }

    ##newmxreg:
    newvxreg <- coredata(cbind(as.zoo(newvxreg)))
    colnames(newvxreg) <- NULL

    ##vxreghat (predictions):
    vxregIndx <- c(max(logewmaIndx)+1):length(coefs)
    vxreghat <-  newvxreg %*% coefs[vxregIndx]
    vxreghat <- c(rep(0,backcastMax),vxreghat)

  } #end vxreg


  ##-----------------------
  ## 3 simulate innov
  ##-----------------------

  ##determine n.sim value:
  if( is.null(n.sim) && is.null(probs) ){
    if( backcastMax==0 ){
      n.sim <- 1
    }else{ 
      n.sim <- ifelse(n.ahead==1, 1, 5000)
    }
  }
  if( is.null(n.sim) && !is.null(probs) ){ n.sim <- 5000 }
     
  ##simulated innovations:
  mZhat <- NULL

  ##the classic bootstrap (innov not provided):
  if( is.null(innov) ){
    zhat <- coredata(residuals.larch(object))
    draws <- runif(n.ahead*n.sim, min=0.5+.Machine$double.eps,
                   max=length(zhat)+0.5+.Machine$double.eps)
    draws <- round(draws, digits=0)
    zhat <- zhat[draws]
  }
    
  ##if user-provided innov:
  if( !is.null(innov) ){
    if(length(innov)!=n.ahead*n.sim){ stop("length(innov) must equal n.ahead*n.sim") }
    zhat <- as.numeric(innov)
  }

  ##matrix of innovations:
  mZhat <- matrix(zhat, n.ahead, n.sim)
  mZhat <- rbind(matrix(NA, backcastMax, NCOL(mZhat)), mZhat) ##modify no. of rows:
  
  
  ##-----------------------
  ## 4 prepare prediction
  ##-----------------------
      
  ##prepare logsd2 predictions:
  logsd2hat <- rep(NA, n.ahead + backcastMax)
  logsd2hat.n <- length(logsd2hat)
  logsd2Fit <- log(coredata(fitted.larch(object)))
  if( backcastMax>0 ){
    logsd2hat[1:backcastMax] <- 
      logsd2Fit[c(length(logsd2Fit)-backcastMax+1):length(logsd2Fit)]
  }
  mLogsd2Hat <- matrix(NA, logsd2hat.n, n.sim) #matrix of logsd2 predictions
  mLogsd2Hat[,1:NCOL(mLogsd2Hat)] <- logsd2hat  #fill with backcast values
  
  ##prepare loge2 predictions:
  loge2hat <- rep(NA, n.ahead + backcastMax)
  loge2hat.n <- length(loge2hat)
  loge2Past <- object$loge2
  if( backcastMax>0 ){
    loge2hat[1:backcastMax] <- 
      loge2Past[c(length(loge2Past)-backcastMax+1):length(loge2Past)]
  }
  mLoge2Hat <- matrix(NA, loge2hat.n, n.sim) #matrix of loge2 predictions
  mLoge2Hat[,1:NCOL(mLoge2Hat)] <- loge2hat  #fill with backcast values
  
  ##prepare epsilon predictions:
  epshat <- rep(NA, n.ahead + backcastMax)
  epshat.n <- length(epshat)
  epsPast <- object$e
  if( backcastMax>0 ){
    epshat[1:backcastMax] <- 
      epsPast[c(length(epsPast)-backcastMax+1):length(epsPast)]
  }
  mEpsilon <- matrix(NA, epshat.n, n.sim) #matrix of epsilon predictions
  mEpsilon[,1:NCOL(mEpsilon)] <- epshat  #fill with backcast values

  ##prepare harch:
  if( harchMax>0 ){
    mLogHarchHat <- matrix(NA, n.ahead+backcastMax, length(harchCoefs))
    colnames(mLogHarchHat) <- object$vXnames[harchIndx]
    mLogHarchHat[1:backcastMax,] <- object$vX[c(NROW(object$vX)-backcastMax+1):NROW(object$vX),harchIndx]
    mLogHarchHat <- as.matrix(mLogHarchHat)
  }
  
  ##prepare asym and/or asymind:
  if( asymMax>0 || asymindMax>0 ){
    zhatIneg <- rep(NA, n.ahead + backcastMax)
    zhatIneg.n <- length(zhatIneg)
    zhatFit <- coredata(residuals(object))
    zhatIneg[1:backcastMax] <- zhatFit[c(length(zhatFit)-backcastMax+1):length(zhatFit)]
    zhatIneg <- as.numeric(zhatIneg<0)
    mZhatIneg <- matrix(NA, zhatIneg.n, n.sim)
    mZhatIneg[,1:NCOL(mZhatIneg)] <- zhatIneg
    mZhatIneg[c(backcastMax+1):NROW(mZhatIneg),] <- matrix(as.numeric(zhat<0),NROW(zhat),NCOL(zhat))
  }
    
  ##prepare log.ewma:
  if( logewmaMax>0 ){
    mLogEwmaHat <- matrix(NA, n.ahead+backcastMax, length(logewmaCoefs))
    colnames(mLogEwmaHat) <- object$vXnames[logewmaIndx]
    mLogEwmaHat[1:backcastMax,] <- object$vX[c(NROW(object$vX)-backcastMax+1):NROW(object$vX),logewmaIndx]
    mLogEwmaHat <- as.matrix(mLogEwmaHat)
  }
  
  ##-------------------------
  ## 5 generate predictions
  ##-------------------------
  
  archTerm <- 0
  harchTerm <- 0
  asymTerm <- 0
  asymindTerm <- 0
  logewmaTerm <- 0
  
  ##loop over n.sim:
  for(j in 1:NCOL(mLogsd2Hat)){
  
    ##loop over backcast+n.ahead:
    for(i in c(backcastMax+1):NROW(mLogsd2Hat) ){
    
      ##compute terms:
      if( archMax>0 ){
        archTerm <- sum( archCoefs*mLoge2Hat[c(i-1):c(i-archMax),j] )
      }
      if( harchMax>0 ){
        for(k in 1:NCOL(mLogHarchHat)){
          sumTerm <- sum( mEpsilon[c(i-harchEval[k]):c(i-1),j] )^2
          sumTerm <- ifelse(sumTerm==0, zero.adj, sumTerm)
          mLogHarchHat[i,k] <- log(sumTerm)
        }
        harchTerm <- sum( coefs[harchIndx] * mLogHarchHat[i,] )
      }
      if( asymMax>0 ){
        asymTerm <- sum( asymCoefs*mLoge2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
      }
      if( asymindMax>0 ){
        asymindTerm <- sum( asymindCoefs*mZhatIneg[c(i-1):c(i-asymindMax),j] )      
      }
      if( logewmaMax>0 ){
        for(k in 1:NCOL(mLogEwmaHat)){
          meanTerm <- mean( mEpsilon[c(i-logewmaEval[k]):c(i-1),j]^2 )
          meanTerm <- ifelse(meanTerm==0, zero.adj, meanTerm)
          mLogEwmaHat[i,k] <- log(meanTerm)
        }
        logewmaTerm <- sum( coefs[logewmaIndx] * mLogEwmaHat[i,] )
      }

      ##compute predictions:
      mLogsd2Hat[i,j] <-
        vconst + archTerm + harchTerm + asymTerm + asymindTerm + 
        logewmaTerm + vxreghat[i]
      mEpsilon[i,j] <- exp(mLogsd2Hat[i,j]/2)*mZhat[i,j]
      if( exp(mLogsd2Hat[i,j]/2)==Inf ){
        message("message: one or more predictions are Inf")
      }
      mLoge2Hat[i,j] <- ifelse( mEpsilon[i,j]==0, log(zero.adj), log(mEpsilon[i,j]^2) )

    } ##end for(i)

  } ##end for(j)

  
  ##-----------------------
  ## 6 prepare output
  ##-----------------------

  ##variance predictions:
  mSd2Hat <- exp( mLogsd2Hat[c(logsd2hat.n-n.ahead+1):logsd2hat.n,] )
  if( n.ahead==1 ){ #rbind() needed when n.ahead=1
    mSd2Hat <- rbind(mSd2Hat)
  }else{ #cbind() needed in case n.sim=1 (when n.ahead>1)
    mSd2Hat <- cbind(mSd2Hat)
  } 
  sd2hat <- as.vector(rowMeans(mSd2Hat))
  if( verbose ){ colnames(mSd2Hat) <- paste0("mSd2Hat.", seq(1,NCOL(mSd2Hat))) }
  
  ##innovations:
  mZhat <- rbind(mZhat[c(logsd2hat.n-n.ahead+1):logsd2hat.n,])
  if( verbose ){ colnames(mZhat) <- paste0("mZhat.", seq(1,n.sim)) }

  ##epsilon:
  mEpsilon <- rbind(mEpsilon[c(logsd2hat.n-n.ahead+1):logsd2hat.n,])
  if( verbose ){ colnames(mEpsilon) <- paste0("mEpsilon.", seq(1,n.sim)) }
  
  ##quantiles of variance predictions:
  mVarianceQs <- NULL
  if( !is.null(probsArg) ){
    mVarianceQs <- matrix(NA, n.ahead, length(probsArg))
    for(i in 1:NROW(mSd2Hat)){
      mVarianceQs[i,] <- quantile(mSd2Hat[i,], probs=probsArg, type=quantile.type)
    }
    colnames(mVarianceQs) <- paste0("q", probsArg)
  }

#FOR THE FUTURE?:
#  ##quantiles of epsilon:
#  mEpsilonQs <- NULL
#  if( !is.null(probsArg) ){
#    mEpsilonQs <- matrix(NA, n.ahead, length(probsArg))
#    for(i in 1:NROW(mEpsilon)){
#      mEpsilonQs[i,] <- quantile(mEpsilon[i,], probs=probsArg, type=quantile.type)
#    }
#    colnames(mEpsilonQs) <- paste0("e", probsArg)
#  }

#  ##quantiles of epsilon^2:
#  TBA


  ##-----------------------
  ## 7 newindex
  ##-----------------------
  
  ##in-sample:
  eInSample <- zoo(object$e, order.by=object$e.index)

  #newindex user-provided:
  if( !is.null(newindex) ){
    if( n.ahead!=length(newindex) ){
      stop("length(newindex) must equal 'n.ahead'")
    }
    newindexInSample <- any( newindex %in% object$e.index )
  }else{ newindexInSample <- FALSE }

  #in-sample index regular:
  if( is.null(newindex) && is.regular(eInSample, strict=TRUE) ){
    endCycle <- cycle(eInSample)
    endCycle <- as.numeric(endCycle[length(endCycle)])
    endYear <- floor(as.numeric(object$e.index[length(object$e)]))
    eFreq <- frequency(eInSample)
    ehataux <- rep(NA, n.ahead+1)
    eDeltat <- deltat(eInSample)
    if( eDeltat==1 && eFreq==1 ){
      ehataux <- zoo(ehataux, order.by=seq(endYear, endYear+n.ahead, by=1))
      eAsRegular <- FALSE
    }else{
      ehataux <- zooreg(ehataux, start=c(endYear, endCycle), frequency=eFreq)
      eAsRegular <- TRUE
    }
    ehataux <- ehataux[-1]
    newindex <- index(ehataux)
  }

  ##neither user-provided nor regular:
  if( is.null(newindex) ){ newindex <- 1:n.ahead }

  ##add index to results:
  if( !is.null(sd2hat) ){ sd2hat <- zoo(sd2hat, order.by=newindex) }
  if( verbose ){
    if( !is.null(mSd2Hat) ){ mSd2Hat <- zoo(mSd2Hat,  order.by=newindex) }
    if( !is.null(mEpsilon) ){ mEpsilon <- zoo(mEpsilon, order.by=newindex) }
    if( !is.null(mZhat) ){ mZhat <- zoo(mZhat, order.by=newindex) }
    if( !is.null(mVarianceQs) ){ mVarianceQs <- zoo(mVarianceQs, order.by=newindex) }
##FOR THE FUTURE?:
#    if( !is.null(mEpsilonQs) ){ mEpsilonQs <- zoo(mEpsilonQs, order.by=newindex) }
#    if( !is.null(mEpsilon2Qs) ){ mEpsilon2Qs <- zoo(mEpsilon2Qs, order.by=newindex) }
  }

      
  ##-----------------------
  ## 8 return result
  ##-----------------------

  ##return variance predictions only:
  if( !verbose ){
    if( is.null(probs) ){
      result <- sd2hat
    }else{    
      result <- cbind(sd2hat,mVarianceQs)
      colnames(result)[1] <- "sd2"
    }
  } #close if( !verbose )

  ##return everything:
  if( verbose ){
    result <- cbind(sd2hat)
    colnames(result) <- "sd2"
    if( !is.null(mVarianceQs) ){ result <- cbind(result, mVarianceQs) }
#FOR THE FUTURE?:
#    if( !is.null(mEpsilonQs) ){ result <- cbind(result,mEpsilonQs) }
    if( !is.null(mSd2Hat) ){ result <- cbind(result, mSd2Hat) }
    if( !is.null(mEpsilon) ){ result <- cbind(result, mEpsilon) }
    if( !is.null(mZhat) ){ result <- cbind(result, mZhat) }
  } #close if( verbose )
    
  ##return the result:
  return(result)
  
} #close predict.larch()  


###########################################################
## print.larch() prints the estimation result
###########################################################

print.larch <- function(x, signif.stars=TRUE, verbose=FALSE, ...)
{
  ##obtain entry names:
  xNames <- names(x)

  ##-------------
  ## gets results
  ##-------------
  
  ##check if there is gets modelling:
  thereIsGetsModelling <- ifelse( "gum.call" %in% xNames, TRUE, FALSE)

  ##if gets modelling:
  if( verbose && thereIsGetsModelling ){

    ##print paths, terminals and retained regressors:
    if( !is.null(x$terminals.results) ){

      ##print gum results:
      if( !is.null(x$gum.results) ){
        cat("\n")
        cat("GUM log-variance equation:\n")
        cat("\n")
        printCoefmat(x$gum.results, cs.ind=c(1,2), tst.ind=c(3),
          signif.stars=TRUE, P.values=TRUE)
      }

      ##print gum diagnostics:
      if( !is.null(x$gum.diagnostics) ){
        cat("\n")
        cat("Diagnostics:\n")
        cat("\n")
        printCoefmat(x$gum.diagnostics, tst.ind=2, signif.stars=TRUE)
        cat("\n")
      }

      ##paths:
      if( length(x$paths)>0 ){
        cat("\n")
        cat("Paths searched:\n")
        cat("\n")
        for(i in 1:length(x$paths)){
          txt <- paste0(x$paths[[i]], collapse=" ")
          txt <- paste0("  Path ", i, ": ", txt)    
          cat(txt, "\n")
        }
      }
  
      ##print terminals:
      cat("\n")
      cat("Terminal models:\n")
      cat("\n")
      print(x$terminals.results)    
  
      ##retained regressors:
      cat("\n")
      cat("Retained regressors (final model):\n")
      if( length(x$specific.spec)==0 ){
        cat("  none\n")
      }else{
        cat(paste0(c("", names(x$specific.spec), "\n"), collapse="  "))
      }
      
    } #end if( print.searchinfo )
  
    ##messages:
    if(!is.null(x$messages)){
      message("\n", appendLF=FALSE)
      message("Messages:", appendLF=TRUE)
      message("\n", appendLF=FALSE)
      message(x$messages)
    }
        
  } #close if( verbose && thereIsGetsModelling )


  ##-------------------------
  ## larch estimation results
  ##-------------------------
  
  ##check if there are larch estimation results:
  thereAreResults <- ifelse("results" %in% xNames, TRUE, FALSE)

  ##header - first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Dependent var.:", x$e.name, "\n")
  cat("Variance-Covariance:", switch(x$vcov.type, robust = "Robust (default)",
    "hac" = "HAC (Newey and West, 1987)"), "\n")
  cat("No. of observations:", x$n, "\n")

  ##header - sample info:
  if( "residuals" %in% xNames ){
    indexResiduals <- index(x$residuals)
    isRegular <- is.regular(x$residuals, strict=TRUE)
    isCyclical <- frequency(x$residuals) > 1
    if(isRegular && isCyclical){
      cycleResiduals <- cycle(x$residuals)
      startYear <- floor(as.numeric(indexResiduals[1]))
      startAsChar <- paste(startYear,
        "(", cycleResiduals[1], ")", sep="")
      endYear <- floor(as.numeric(indexResiduals[length(indexResiduals)]))
      endAsChar <- paste(endYear,
        "(", cycleResiduals[length(indexResiduals)], ")", sep="")
    }else{
      startAsChar <- as.character(indexResiduals[1])
      endAsChar <- as.character(indexResiduals[length(indexResiduals)])
    }
    cat("Sample:", startAsChar, "to", endAsChar, "\n")
  } #end if( "residuals" %in% xNames )

  ##print results:
  if( thereAreResults ){
    cat("\n")
    cat("Log-variance equation:\n")
    cat("\n")
    printCoefmat(x$results, signif.stars=signif.stars)
  }else{  
    cat("\n")
    cat("   No estimation results\n")
  }

  ##create goodness-of-fit matrix:
  gof <- matrix(NA, 1, 1)
  rownames(gof) <- paste0("Log-lik.(n=", x$n, "):")
  colnames(gof) <- ""
  gof[1,1] <- as.numeric(x$logl)

  ##print diagnostics and fit:
  if( !is.null(x$diagnostics) ) {
    cat("\n")
    cat("Diagnostics and fit:\n")
    cat("\n")
    printCoefmat(x$diagnostics, tst.ind=2,
      signif.stars=signif.stars, has.Pvalue=TRUE)
    printCoefmat(gof, digits=6, signif.stars=signif.stars)
  }

} #end print.larch()


###########################################################
## residuals.larch()
###########################################################

residuals.larch <- function(object, ...){ return(object$residuals) }


###########################################################
## summary.larch()
###########################################################

summary.larch <- function(object, ...){ summary.default(object) }


###########################################################
## toLatex.larch()
###########################################################

toLatex.larch <- function(object, ...)
{
  printtex(object, fitted.name="\\ln\\sigma_t^2", ...)
}


###########################################################
## vcov.larch()
###########################################################

vcov.larch <- function(object, ...){ return(object$vcov) }
