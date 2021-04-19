isat.lm <- function(y, ar=NULL, ewma=NULL, iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=NULL,
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
  
  
  x_vars <- data[,2:ncol(data)]
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
}


isat.arx <- function(
  y,
  ar = NULL,
  ewma = NULL,
  iis = FALSE,
  sis = TRUE,
  tis = FALSE,
  uis = FALSE,
  blocks = NULL,
  ratio.threshold = 0.8,
  max.block.size = 30,
  t.pval = 0.001,
  wald.pval = t.pval,
  vcov.type = c("ordinary", "white", "newey-west"),
  do.pet = FALSE,
  ar.LjungB = NULL,
  arch.LjungB = NULL,
  normality.JarqueB = NULL,
  info.method = c("sc", "aic", "hq"),
  user.diagnostics = NULL,
  user.estimator = NULL,
  gof.function = NULL,
  gof.method = c("min", "max"),
  include.gum = NULL,
  include.1cut = FALSE,
  include.empty = FALSE,
  max.paths = NULL,
  parallel.options = NULL,
  turbo = FALSE,
  tol = 1e-07,
  LAPACK = FALSE,
  max.regs = NULL,
  print.searchinfo = TRUE,
  plot = NULL,
  alarm = FALSE,
  ...
){
  
  # warnings and checks (mainly for variance specification)
  if(!is.null(y$variance.results)){
    warning("Input object contains variance specification. Note that 'isat' is not configured for variance specifications.\nVariance specification in 'isat' are dropped.")
  }
  
  # Check if one of these arguments is explicitly supplied to the function
  # if not, then check if the original item has this arguemnt supplied
  # if it does, take the setting of the original object
  # if it does not, then take the default
  if(missing(ar)){ar <- if(is.null(y$aux$arguments[["ar"]])) {NULL} else{y$aux$arguments[["ar"]]}}
  if(missing(vcov.type)){vcov.type <- y$aux[["vcov.type"]]}
  if(missing(normality.JarqueB)){normality.JarqueB <- if(is.null(y$aux$arguments[["normality.JarqueB"]])){FALSE}else{y$aux$arguments[["normality.JarqueB"]]}}
  if(missing(user.estimator)){user.estimator <- if(is.null(y$aux$arguments[["user.estimator"]])){NULL}else{y$aux$arguments[["user.estimator"]]}}
  if(missing(user.diagnostics)){user.diagnostics <- if(is.null(y$aux$arguments[["user.diagnostics"]])){NULL}else{y$aux$arguments[["user.diagnostics"]]}}
  if(missing(LAPACK)){LAPACK <- if(is.null(y$aux$arguments[["LAPACK"]])){FALSE}else{y$aux$arguments[["LAPACK"]]}}
  if(missing(plot)){plot <- if(is.null(y$aux$arguments[["plot"]])){NULL}else{y$aux$arguments[["plot"]]}}
  if(missing(tol)){tol <- if(is.null(y$aux$arguments[["tol"]])){1e-07}else{y$aux$arguments[["tol"]]}}
  
  mxreg <- y$aux$mX
  colnames(mxreg) <- y$aux$mXnames
  
  out <- isat.default(
    y$aux$y,
    FALSE,
    ar,
    ewma,
    mxreg,
    iis,
    sis,
    tis,
    uis,
    blocks,
    ratio.threshold,
    max.block.size,
    t.pval,
    wald.pval,
    vcov.type,
    do.pet,
    ar.LjungB,
    arch.LjungB,
    normality.JarqueB,
    info.method,
    user.diagnostics,
    user.estimator,
    gof.function,
    gof.method,
    include.gum,
    include.1cut,
    include.empty,
    max.paths,
    parallel.options,
    turbo,
    tol,
    LAPACK,
    max.regs,
    print.searchinfo,
    plot,
    alarm
  )
  # if(plot){
  #   plot(out)
  # }
  
  return(out)
}



arx.isat <- function(y, ar=NULL, ewma=NULL,
                     vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
                     zero.adj=0.1, vc.adj=TRUE,
                     vcov.type=c("ordinary", "white", "newey-west"),
                     qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
                     user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, plot=NULL, ...){
  
  # Check if one of these arguments is explicitly supplied to the function
  # if not, then check if the original item has this arguemnt supplied
  # if it does, take the setting of the original object
  # if it does not, then take the default
  if(missing(ar)){ar <- if(is.null(y$aux$arguments[["ar"]])) {NULL} else{y$aux$arguments[["ar"]]}}
  if(missing(vcov.type)){vcov.type <- y$aux[["vcov.type"]]}
  if(missing(normality.JarqueB)){normality.JarqueB <- if(is.null(y$aux$arguments[["normality.JarqueB"]])){FALSE}else{y$aux$arguments[["normality.JarqueB"]]}}
  if(missing(user.estimator)){user.estimator <- if(is.null(y$aux$arguments[["user.estimator"]])){NULL}else{y$aux$arguments[["user.estimator"]]}}
  if(missing(user.diagnostics)){user.diagnostics <- if(is.null(y$aux$arguments[["user.diagnostics"]])){NULL}else{y$aux$arguments[["user.diagnostics"]]}}
  if(missing(LAPACK)){LAPACK <- if(is.null(y$aux$arguments[["LAPACK"]])){FALSE}else{y$aux$arguments[["LAPACK"]]}}
  if(missing(plot)){plot <- if(is.null(y$aux$arguments[["plot"]])){NULL}else{y$aux$arguments[["plot"]]}}
  if(missing(tol)){tol <- if(is.null(y$aux$arguments[["tol"]])){1e-07}else{y$aux$arguments[["tol"]]}}
  
  browser()
  dep_var <- matrix(y$aux$y)
  colnames(dep_var) <- y$aux$y.name 
  
  mxreg <- y$aux$mX
  
  out <- arx(
    dep_var,
    FALSE,
    ar,
    ewma,
    mxreg = mxreg,
    vc,
    arch,
    asym,
    log.ewma,
    vxreg,
    zero.adj,
    vc.adj,
    vcov.type,
    qstat.options,
    normality.JarqueB,
    user.estimator = user.estimator,
    user.diagnostics,
    tol,
    LAPACK,
    FALSE
  )
  
  if(!is.null(plot) && plot){
    plot(out)
  }
  
  return(out)
}



gets.isat <- function(x, t.pval=0.05, wald.pval=t.pval, vcov.type = NULL,
                      do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
                      arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
                      user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
                      gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
                      include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
                      turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE, ...)
{
  
  # Check if one of these arguments is explicitly supplied to the function
  # if not, then check if the original item has this argument supplied
  # if it does, take the setting of the original object
  # if it does not, then take the default
  if(missing(t.pval)){t.pval <- if(is.null(x$aux$arguments[["t.pval"]])) {0.05} else{x$aux$arguments[["t.pval"]]}}
  if(missing(vcov.type)){vcov.type <- x$aux[["vcov.type"]]}
  if(missing(normality.JarqueB)){normality.JarqueB <- if(is.null(x$aux$arguments[["normality.JarqueB"]])){FALSE}else{x$aux$arguments[["normality.JarqueB"]]}}
  if(missing(user.diagnostics)){user.diagnostics <- if(is.null(x$aux$arguments[["user.diagnostics"]])){NULL}else{x$aux$arguments[["user.diagnostics"]]}}
  if(missing(plot)){plot <- if(is.null(x$aux$arguments[["plot"]])){NULL}else{x$aux$arguments[["plot"]]}}
  if(missing(tol)){tol <- if(is.null(x$aux$arguments[["tol"]])){1e-07}else{x$aux$arguments[["tol"]]}}
  
  
  object <- arx(x, plot = FALSE)

  out <- getsm(
    object,
    t.pval,
    wald.pval,
    vcov.type,
    do.pet,
    ar.LjungB,
    arch.LjungB,
    normality.JarqueB,
    user.diagnostics,
    info.method,
    gof.function,
    gof.method,
    keep,
    include.gum,
    include.1cut,
    include.empty,
    max.paths,
    tol,
    turbo,
    print.searchinfo,
    plot,
    alarm
  )
  return(out)
}


arx.lm <- function(y, ar=NULL, ewma=NULL,
                   vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
                   zero.adj=0.1, vc.adj=TRUE,
                   vcov.type=c("ordinary", "white", "newey-west"),
                   qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
                   user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, plot=NULL, ...){
  # Checks
  if(!is.null(y$weights)){stop("Usage of weights is not yet implemented in arx. Please estimate the lm y without weights.")}
  
  data <- stats::model.frame(y)
  dep_var <- matrix(data[,1])
  colnames(dep_var) <- names(data)[1] 
  
  
  x_vars <- data[,2:ncol(data)]
  
  mxreg <- stats::model.matrix(y$terms, data)
  invalid_variables_in_lm <- names(y$coef)[is.na(y$coef)]
  mxreg <- mxreg[,!colnames(mxreg) %in% invalid_variables_in_lm]
  
  # Deal with intercept
  mc <- ifelse(attr(y$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept
  
  out <- arx.default(
    dep_var,
    mc,
    ar,
    ewma,
    mxreg,
    vc,
    arch,
    asym,
    log.ewma,
    vxreg,
    zero.adj,
    vc.adj,
    vcov.type,
    qstat.options,
    normality.JarqueB,
    user.estimator,
    user.diagnostics,
    tol,
    LAPACK,
    plot
  )
  return(out)
}




gets.lm <- function(x, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
                    do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
                    arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
                    user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
                    gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
                    include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
                    turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE, ...)
{
  arx_object <- arx.lm(x)
  
  out <- getsm(
    arx_object,
    t.pval,
    wald.pval,
    vcov.type,
    do.pet,
    ar.LjungB,
    arch.LjungB,
    normality.JarqueB,
    user.diagnostics,
    info.method,
    gof.function,
    gof.method,
    keep,
    include.gum,
    include.1cut,
    include.empty,
    max.paths,
    tol,
    turbo,
    print.searchinfo,
    plot,
    alarm
  )
  return(out)
}



