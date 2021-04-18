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
  if(missing(ar)){ar <- if(is.null(y$call[["ar"]])) {NULL} else{y$call[["ar"]]}}
  if(missing(vcov.type)){vcov.type <- y$aux[["vcov.type"]]}
  if(missing(normality.JarqueB)){normality.JarqueB <- if(is.null(y$call[["normality.JarqueB"]])){FALSE}else{y$call[["normality.JarqueB"]]}}
  if(missing(user.estimator)){user.estimator <- if(is.null(y$call[["user.estimator"]])){NULL}else{y$call[["user.estimator"]]}}
  if(missing(user.diagnostics)){user.diagnostics <- if(is.null(y$call[["user.diagnostics"]])){NULL}else{y$call[["user.diagnostics"]]}}
  if(missing(LAPACK)){LAPACK <- if(is.null(y$call[["LAPACK"]])){FALSE}else{y$call[["LAPACK"]]}}
  if(missing(plot)){plot <- if(is.null(y$call[["plot"]])){NULL}else{y$call[["plot"]]}}
  if(missing(tol)){tol <- if(is.null(y$call[["tol"]])){1e-07}else{y$call[["tol"]]}}
  
  out <- isat.default(
    y = y$aux$y,
    mxreg = y$aux$mX,
    mc = FALSE,
    ar = ar,
    ewma = ewma,
    iis = iis,
    sis = sis,
    tis = tis,
    uis = uis,
    blocks = blocks,
    ratio.threshold = ratio.threshold,
    max.block.size = max.block.size,
    t.pval = t.pval,
    wald.pval = wald.pval,
    vcov.type = vcov.type,
    do.pet = do.pet,
    ar.LjungB = ar.LjungB,
    arch.LjungB = arch.LjungB,
    normality.JarqueB = normality.JarqueB,
    info.method = info.method,
    user.diagnostics = user.diagnostics,
    user.estimator = user.estimator,
    gof.function = gof.function,
    gof.method = gof.method,
    include.gum = include.gum,
    include.1cut = include.1cut,
    include.empty = include.empty,
    max.paths = max.paths,
    parallel.options = parallel.options,
    turbo = turbo,
    tol = tol,
    LAPACK = LAPACK,
    max.regs = max.regs,
    print.searchinfo = print.searchinfo,
    plot = FALSE,
    alarm = alarm
  )
  if(plot){
    plot(out)
  }
  
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
  if(missing(ar)){ar <- if(is.null(y$call[["ar"]])) {NULL} else{y$call[["ar"]]}}
  if(missing(vcov.type)){vcov.type <- y$aux[["vcov.type"]]}
  if(missing(normality.JarqueB)){normality.JarqueB <- if(is.null(y$call[["normality.JarqueB"]])){FALSE}else{y$call[["normality.JarqueB"]]}}
  if(missing(user.estimator)){user.estimator <- if(is.null(y$call[["user.estimator"]])){NULL}else{y$call[["user.estimator"]]}}
  if(missing(user.diagnostics)){user.diagnostics <- if(is.null(y$call[["user.diagnostics"]])){NULL}else{y$call[["user.diagnostics"]]}}
  if(missing(LAPACK)){LAPACK <- if(is.null(y$call[["LAPACK"]])){FALSE}else{y$call[["LAPACK"]]}}
  if(missing(plot)){plot <- if(is.null(y$call[["plot"]])){NULL}else{y$call[["plot"]]}}
  if(missing(tol)){tol <- if(is.null(y$call[["tol"]])){1e-07}else{y$call[["tol"]]}}
  
  
  dep_var <- data.frame(y$aux$y)
  names(dep_var) <- y$aux$y.name 
  
  mxreg <- y$aux$mX
  
  out <- arx(y = dep_var, mxreg = mxreg, mc = FALSE, 
             ar = ar, ewma = ewma,
             vc = vc, arch = arch, asym = asym, log.ewma = log.ewma, vxreg = vxreg,
             zero.adj = zero.adj, vc.adj = vc.adj,
             vcov.type = vcov.type,
             qstat.options = qstat.options, normality.JarqueB = normality.JarqueB, user.estimator = user.estimator,
           user.diagnostics = user.diagnostics, tol = tol, LAPACK = LAPACK, plot = FALSE)
  if(plot){
    plot(arx)
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
  
  if(missing(vcov.type)){vcov.type <- x$aux$vcov.type}
  if(missing(info.method)){info.method <- "sc"}
  
  object <- arx(x, vcov.type = vcov.type)
 
  out <- getsm(object,
               t.pval = t.pval, wald.pval = wald.pval, vcov.type = vcov.type,do.pet = do.pet, 
               ar.LjungB = ar.LjungB,arch.LjungB = arch.LjungB,normality.JarqueB = normality.JarqueB,
               user.diagnostics = user.diagnostics, info.method = info.method,
               gof.function = gof.function, gof.method = gof.method, keep = keep, 
               include.1cut= include.1cut, include.empty = include.empty, 
               max.paths = max.paths, tol = tol,turbo = turbo, print.searchinfo = print.searchinfo, 
               plot = plot, alarm = alarm)
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
  dep_var <- data.frame(y = data[,1])
  names(dep_var) <- names(data)[1] 
  
  
  x_vars <- data[,2:ncol(data)]
  mxreg <- stats::model.matrix(y$terms, data)
  
  # Deal with intercept
  mc <- ifelse(attr(y$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept
  
  out <- arx.default(
    y = dep_var,
    mxreg = mxreg,
    mc = mc,
    ar = ar,
    ewma = ewma,
    vc = vc,
    arch = arch,
    asym = asym,
    log.ewma = log.ewma,
    vxreg = vxreg,
    zero.adj = zero.adj,
    vc.adj = vc.adj,
    vcov.type = vcov.type,
    qstat.options = qstat.options,
    normality.JarqueB = normality.JarqueB,
    user.estimator =  user.estimator,
    user.diagnostics = user.diagnostics,
    tol = tol,
    LAPACK = LAPACK,
    plot =  plot
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
  arx_object <- arx(x)
  out <- getsm(object = arx_object, t.pval = t.pval, wald.pval = wald.pval,vcov.type =  vcov.type,do.pet = do.pet, ar.LjungB = ar.LjungB,
               arch.LjungB = arch.LjungB,normality.JarqueB = normality.JarqueB,user.diagnostics = user.diagnostics, info.method = info.method,
               gof.function = gof.function, gof.method = gof.method, keep = keep, include.gum = include.gum,
               include.1cut = include.1cut, include.empty = include.empty, max.paths = max.paths, tol = tol,
               turbo = turbo, print.searchinfo = print.searchinfo, plot = plot, alarm = alarm)
  return(out)
}



