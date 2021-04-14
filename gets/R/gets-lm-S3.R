
isat.lm <- function(lmobject, ar=NULL, ewma=NULL, iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=NULL,
                    ratio.threshold=0.8, max.block.size=30, t.pval=0.001,
                    wald.pval=t.pval, vcov.type=c("ordinary", "white", "newey-west"),
                    do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL,
                    normality.JarqueB=NULL, info.method=c("sc", "aic", "hq"), 
                    user.diagnostics=NULL, user.estimator=NULL, gof.function=NULL, 
                    gof.method=c("min","max"), include.gum=NULL,
                    include.1cut=FALSE, include.empty=FALSE, max.paths=NULL,
                    parallel.options=NULL, turbo=FALSE, tol=1e-07, LAPACK=FALSE,
                    max.regs=NULL, print.searchinfo=TRUE, plot=NULL, alarm=FALSE){
  
  # Checks
  if(!is.null(lmobject$weights)){stop("Usage of weights is not yet implemented in isat. Please estimate the lm object without weights.")}
  
  data <- stats::model.frame(lmobject)
  y <- data.frame(y = data[,1])
  names(y) <- names(data)[1] 
  
  
  x_vars <- data[,2:ncol(data)]
  mxreg <- stats::model.matrix(lmobject$terms, data)
  # Deal with intercept
  mc <- ifelse(attr(lmobject$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept
  
  out <- gets::isat(y = y, mxreg = mxreg, mc = mc, 
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
}

getsm.lm <- function(lmobject, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
                     do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
                     arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
                     user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
                     gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
                     include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
                     turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE){
  arx_object <- arx.lm(lmobject)
  getsm(arx_object, t.pval, wald.pval, vcov.type,do.pet, ar.LjungB,arch.LjungB, 
        normality.JarqueB,user.diagnostics, info.method,gof.function, gof.method, keep, include.gum,
        include.1cut, include.empty, max.paths, tol,turbo, print.searchinfo, plot, alarm)
}


arx.lm <- function(lmobject, ar=NULL, ewma=NULL, 
                   vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
                   zero.adj=0.1, vc.adj=TRUE,
                   vcov.type=c("ordinary", "white", "newey-west"),
                   qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
                   user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, plot=NULL){
  # Checks
  if(!is.null(lmobject$weights)){stop("Usage of weights is not yet implemented in arx. Please estimate the lm object without weights.")}
  
  data <- stats::model.frame(lmobject)
  y <- data.frame(y = data[,1])
  names(y) <- names(data)[1] 
  
  
  x_vars <- data[,2:ncol(data)]
  mxreg <- stats::model.matrix(lmobject$terms, data)
  
  # Deal with intercept
  mc <- ifelse(attr(lmobject$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept
  out <- arx(y = y, mxreg = mxreg, mc = mc, 
             ar, ewma, vc, arch, asym, log.ewma, vxreg,
             zero.adj, vc.adj,
             vcov.type,
             qstat.options, normality.JarqueB, user.estimator,
             user.diagnostics, tol, LAPACK, plot)
  return(out)
}


isat <- function(x) {
  UseMethod("isat")
}
getsm <- function(x) {
  UseMethod("getsm")
}
arx <- function(x) {
  UseMethod("arx")
}
