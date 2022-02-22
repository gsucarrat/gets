distorttest <- function(x, coef="all"){
  
  # Checks
  
  if (is.null(x)){stop("Object is NULL - please make sure to pass an isat object.")}
  if (class(x)!="isat"){ stop("x must be an isat object")}
  if (any(x$call$sis, x$call$tis, x$call$uis)==TRUE){stop("Test only valid for iis - not valid for sis, uis, or tis.")} 
  if (all(coef != "all") && any(!coef %in% names(coef(x)))){stop("The 'coef' variable or vector contains one or more regressors not in the isat object.")}
  if (x$call$iis != TRUE){stop("The isat object has not selected iis = TRUE. This is necessary for this test. Re-estimate isat with iis = TRUE.")}
  # WRONG CHECK - Therefore commented out: if (is.null(x$ISnames)){stop("IIS did not identify any Indicators (Outliers). Therefore OLS = IIS and no distortion is detectable.")}
  
  ISnames <- c(x$ISnames[grep("iis", x$ISnames)])
  noutl = length(ISnames)
  t.pval <- x$aux$t.pval
  # T <- x$n # MORITZ: Want to change this. Using T is not good practice - update: this is not used anywhere else so commented out
  
  alpha <-   t.pval
  c <- abs(qnorm(alpha/2))
  fc <- dnorm(c)
  psi <- pnorm(c) - pnorm(-c)
  tau <- psi - 2 * c * fc
  
  varsigmac = tau / psi
  rhoc1 = psi^2 / (4*(c^2)*(fc^2) + 4*tau*c*fc + tau)
  rhocinfinity = (psi - 2*c*fc)^2 / tau
  avar1 = ((2*c*fc - psi)^2 + 2*tau*(2*c*fc - psi) + tau) /psi^2
  avarinfinity = ((2*c*fc - psi)^2 + 2*tau*(2*c*fc - psi) + tau) / (2*c*fc - psi)^2
  
  
  x.date <- isatdates(x)$iis$date
  
  # changed by M-orca because there is a problem when the y.index is not starting at 1 e.g. if its years
  # September 2021
  # n.null <-  x$aux$y.index[!(x$aux$y.index %in% x.date)]
  n.null <-  !(x$aux$y.index %in% x.date)
  
  keep <- which(!(names(coef(x)) %in%  x$ISnames))
  
  y.null <-  x$aux$y[n.null]
  
  
  if (!is.null(x$call$ar)){
    ar.call <- eval(x$call$ar)
  }
  
  if (!is.null(x$call$ar)){
    mx.null <- x$aux$mX[(n.null-max(eval(x$call$ar))),keep] # M-orca note: this won't work with panels!!
  } else{
    mx.null <- x$aux$mX[n.null,keep]
  }
  
  
  y <-  x$aux$y
  mx <- x$aux$mX[,keep]
  
  nOLS <- length(x$aux$y)
  nOLS_rob <- length(x$aux$y) - length(x$ISnames)
  
  # This is the GUM
  # Save original arx mc warning setting and disable it here
  tmpmc <- getOption("mc.warning")
  options(mc.warning = FALSE)
  
  ols.y <- arx(y, mxreg=mx, mc=FALSE, plot=FALSE, ar=FALSE)
  
  # Set the old arx mc warning again
  options(mc.warning = tmpmc)
  
  if (NROW(coef)==1){
    if (coef=="all"){ #if testing on all coefficients
      
      betaOLS1 <- coef(x)[keep]
      betaOLS <- coef(ols.y)
      
      V <- (nOLS_rob/nOLS) * (rhoc1)^(-1)*(varsigmac)^(-1) * x$vcov.mean[keep, keep]
      
    } else { #testing on subset of coefficients
      
      betaOLS1 <- coef(x)[coef]
      betaOLS <- coef(ols.y)[coef]
      
      V <- (nOLS_rob/nOLS) * (rhoc1)^(-1)*(varsigmac)^(-1) * x$vcov.mean[coef, coef]
      
    }
    
  } else {
    
    betaOLS1 <- coef(x)[coef]
    betaOLS <- coef(ols.y)[coef]
    
    V <- (nOLS_rob/nOLS) * (rhoc1)^(-1)*(varsigmac)^(-1) * x$vcov.mean[coef, coef]
    
  }
  
  if (NROW(coef)==1){
    if ( coef=="all"){
      rel.df <- length(coef(x)[keep])  #degrees of freedom
    } else {
      rel.df <- length(coef)
    }
    
  } else {
    rel.df <- length(coef)
  }
  
  
  cf_diff <- (betaOLS1 - betaOLS) #difference between coefficients
  
  esigma2MinverseOLS1 = nOLS * rhoc1 * V
  eavarOLS10 = avar1 * esigma2MinverseOLS1
  
  
  if (any(!cf_diff==0)){ #if the difference between OLS and IIS coefficients is not zero
    
    HtestOLS10 = as.numeric(nOLS * t(cf_diff) %*% solve(eavarOLS10) %*% as.vector(cf_diff))
    
  } else {
    HtestOLS10 <- 0
  }
  
  p.test <- pchisq(HtestOLS10, df = rel.df, lower.tail = FALSE)
  #attr(cf_diff,"names") <- "Difference in Estimates"
  attr(HtestOLS10,"names") <- "X-squared"
  rval_chi <- list(
    #statistic = HtestOLS10, p.value = p.test, estimate=NULL, null.value = NULL, # Felix
    statistic = HtestOLS10, p.value = p.test, estimate=cf_diff, null.value = NULL, df = rel.df, # Moritz
    # alternative = NULL, # Felix
    alternative = "True difference between IIS and OLS estimates is not equal to 0.", # M-Orca attempt
    
    method="Jiao-Pretis-Schwarz Outlier Distortion Test", # Moritz
    #data.name="Difference between IIS and OLS Estimates", # Felix
    
    data.name=deparse(substitute(x)), # M-orca changed
    coef.diff = cf_diff, 
    var.diff = eavarOLS10, 
    iis=x, 
    ols=ols.y)
  attr(rval_chi, "class") <- c("htest", "distorttest")
  
  out <- return(rval_chi)
  
}


distorttestboot <- function(
  x,
  nboot,
  clean.sample = TRUE,
  parametric = FALSE,
  scale.t.pval = 1,
  parallel.options = NULL,
  quantiles = c(0.9, 0.95, 0.99),
  #parallel = FALSE,
  #ncore = detectCores()[1] - 1,
  ...
){
  
  
  # x <- is1
  # nboot=199
  #
  ####compute distortion on full model
  #dist.full <- distorttest(x)
  
  if(!identical(class(x),c("htest","distorttest"))){
    if (class(x)=="isat") {
      data.name <- deparse(substitute(x))
      x <- distorttest(x)
      x$data.name <- data.name
    } else {stop("Unrecognised class of x. The function needs either an isat object or the output of the distorttest() function.")}
  }
  
  
  dist.full <- x
  #coef_name <- attr(x$coef.diff, "names") # FUTURE DEVELOPMENT. Current error stopping below: 
  if(!identical(attr(x$coef.diff,"names"),x$iis$aux$mXnames[!x$iis$aux$mXnames %in% x$iis$ISnames])){stop("Bootstrap on distorttest currently must be coef = 'all'.")}
  x <- x$iis
  #names(dist.full$coef.diff)
  
  out.full <- outliertest(x)
  
  is0.date <- isatdates(x)$iis$index
  N <- x$aux$y.n
  # define the set over which to sample the bootstraps from as all observations apart from those where IIS identified outliers
  n.null <- seq(1:N)[!(seq(1:N) %in% is0.date)]
  
  boot.samples <- list()
  
  if (parametric){
    clean.sample <- TRUE
  }
  
  for (i in 1:nboot){
    if (clean.sample == TRUE){
      boot.samp <- sample(n.null, N, replace=TRUE)
    } else {
      boot.samp <- sample(seq(1:N), N, replace=TRUE)
    }
    boot.samples[[i]] <- boot.samp
    
  } #i closed
  
  
  boot.tpval <- x$aux$t.pval*scale.t.pval #bootstrap level of significance of selection
  
  
  
  ### function used in parallel loop
  dist.boot.temp <- function(y.boot, x.boot, boot.tpval, ...){
    is.boot <- isat(y.boot, mxreg=x.boot, mc=FALSE, t.pval=boot.tpval, iis=TRUE, sis=FALSE,  print.searchinfo=FALSE, ...)
    dist.boot <- distorttest(is.boot)
    out.boot <- outliertest(is.boot)
    dist.res <- c(
      dist.boot$coef.diff,
      dist.boot$statistic,
      
      out.boot$proportion$estimate,
      out.boot$proportion$statistic,
      out.boot$count$estimate,
      out.boot$count$statistic
    )
    names(dist.res) <- c(names(dist.boot$coef.diff), "dist", "prop", "prop.test", "count", "count.test")
    return(dist.res)
  }
  
  overall_fun <- function(i){
    boot.samp <-  boot.samples[[i]]
    
    if (parametric==TRUE){ #parametric bootstrap (residual resampling)
      x.boot <- x$aux$mX[, !(colnames(x$aux$mX) %in% x$ISnames)]
      res.boot <- as.vector(x$residuals)[boot.samp]
      y.boot <-   x.boot %*% coef(x)[!(colnames(x$aux$mX) %in% x$ISnames)] + res.boot
      
    } else { #nonparametric
      y.boot <- x$aux$y[boot.samp]
      x.boot <- x$aux$mX[boot.samp, !(colnames(x$aux$mX) %in% x$ISnames)]
    }
    
    tempMatrix =   dist.boot.temp(y.boot, x.boot, boot.tpval)
  }
  
  
  
  lapply(1:nboot, overall_fun)
  
  
  if (!is.null(parallel.options)){
    if(!is.numeric(parallel.options)){stop("parallel.options must either be NULL or a numeric value.")}
    
    if(is.numeric(parallel.options)){
      clusterSpec <- parallel.options
      OScores <- detectCores()
      if(parallel.options > OScores){
        stop("parallel.options > number of cores/threads")
      }
    }
    
    #ncore=7
    #  boot.tpval <- 0.05
    
    cl <- makeCluster(parallel.options) #not to overload your computer
    
    coefdist.sample.list <- parLapply(cl = cl, 
                                      X = 1:nboot,
                                      fun = overall_fun)
    
    coefdist.sample <- do.call(rbind, coefdist.sample.list)
    
  } else { #if not using parallel
    
    coefdist.sample.list <- lapply(
      X = 1:nboot,
      FUN = overall_fun)
    
    coefdist.sample <- do.call(rbind, coefdist.sample.list)
    
  } #parallel if closed
  
  # end.time <- Sys.time()
  # time.diff <-  end.time - start.time
  # print(paste("Boot Complete in", sep=""))
  # print(time.diff)
  
  # start.time <- Sys.time()
  
  coefdist.sample <- as.data.frame(coefdist.sample)
  coefdist.res <- data.frame(matrix(NA, nrow=nboot, ncol=1))
  names(coefdist.res) <- "boot"
  coefdist.res$boot <- seq(1:nboot)
  if (length(names(dist.full$coef.diff)) > 1){
    coefdist.res$L2 <- apply(coefdist.sample[,names(dist.full$coef.diff)], 1, function(x)  sqrt(sum(x^2)))
    coefdist.res$L1 <- apply(coefdist.sample[,names(dist.full$coef.diff)], 1, function(x)  sum(abs(x)))
  } else {
    coefdist.res$L2 <- sqrt(coefdist.sample[,1]^2)
    coefdist.res$L1 <- abs(coefdist.sample[,1])
  }
  
  coefdist.res$dist <- coefdist.sample$dist
  
  coefdist.res$prop <- coefdist.sample$prop
  coefdist.res$count <- coefdist.sample$count
  coefdist.res$prop.test <- coefdist.sample$prop.test
  coefdist.res$count.test <- coefdist.sample$count.test
  
  ############
  # Distortion
  ############
  L2.full <- sqrt(sum(dist.full$coef.diff^2))
  L1.full <- sum(abs(dist.full$coef.diff))
  dist.full.stat <- dist.full$statistic
  
  boot.q.L2 <- quantile(coefdist.res$L2, probs = quantiles) # sample estimate L2 for distortion
  boot.q.L1 <- quantile(coefdist.res$L1, probs = quantiles) # sample estimate L1 for distortion
  boot.q.dist <- quantile(coefdist.res$dist, probs = quantiles)
  
  boot.p.L2 <- sum(coefdist.res$L2 > L2.full)/nboot
  boot.p.L1 <- sum(coefdist.res$L1 > L1.full)/nboot
  boot.p.dist <- sum(coefdist.res$dist > dist.full.stat)/nboot
  
  ############
  # Proportion
  ############
  prop.full <- out.full$proportion$estimate
  prop.full.stat <- out.full$proportion$statistic
  
  count.full <- out.full$count$estimate
  count.full.stat <- out.full$count$statistic
  
  boot.q.prop <- quantile(coefdist.res$prop, probs = quantiles) # sample estimate for proportion
  boot.q.count <- quantile(coefdist.res$count, probs = quantiles) # sample estimate for count
  boot.q.prop.test <- quantile(coefdist.res$prop.test, probs = quantiles) # proportion stat
  boot.q.count.test <- quantile(coefdist.res$count.test, probs = quantiles) # count stat
  
  boot.p.prop <- 2*min(sum(coefdist.res$prop > prop.full)/nboot,  sum(coefdist.res$prop <= prop.full)/nboot) #check the p-values here.
  boot.p.count <- 2*min(sum(coefdist.res$count > count.full)/nboot, sum(coefdist.res$count <= count.full)/nboot)
  
  boot.p.prop.test <- sum(abs(coefdist.res$prop.test) > abs(prop.full.stat))/(nboot)
  boot.p.count.test <- sum(abs(coefdist.res$count.test) > abs(count.full.stat))/(nboot)
  
  ####
  out <- list()
  
  out$dist.full <- dist.full
  out$dist.full$method <- "Jiao-Pretis-Schwarz Outlier Distortion Test (Full Sample)"
  
  out$coefdist.res <- coefdist.res
  
  ############
  # Distortion
  ############
  out$L2.full <- L2.full
  out$boot.q.L2 <- boot.q.L2
  out$boot.p.L2 <- boot.p.L2
  out$L1.full <- L1.full
  out$boot.q.L1 <- boot.q.L1
  out$boot.p.L1 <- boot.p.L1
  out$boot.q.dist <- boot.q.dist
  out$boot.p.dist <- boot.p.dist
  
  
  ############
  # Proportion
  ############
  out$prop.full <- prop.full
  out$boot.q.prop <- boot.q.prop
  out$boot.p.prop <- boot.p.prop
  
  out$prop.full.stat <- prop.full.stat
  out$boot.q.prop.stat <- boot.q.prop.test
  out$boot.p.prop.stat <- boot.p.prop.test
  
  out$count.full <- count.full
  out$boot.q.count <- boot.q.count
  out$boot.p.count <- boot.p.count
  
  out$count.full.stat <- count.full.stat
  out$boot.q.count.stat <- boot.q.count.test
  out$boot.p.count.stat <- boot.p.count.test
  
  
  out$args <- list(parametric = parametric, 
                   nboot = nboot, 
                   scale.t.pval = scale.t.pval,
                   clean.sample = clean.sample)
  
  class(out) <- "distorttestboot"
  
  return(out)
  
}




print.distorttestboot <- function(x, print.proportion = FALSE, ...){
  
  print(x$dist.full)
  cat("\n")
  cat("------ Bootstrap Results ------")
  cat("\n")
  cat(paste0("Number of Bootstrap Replications: ",x$args$nboot))
  cat("\n")
  # cat("\n")
  # mOutl_p <- matrix(NA, 2, 2)
  # colnames(mOutl_p) <- c("Stat.", "p-value")
  # rownames(mOutl_p) <- c("Jiao-Pretis Bootstrapped Outlier Proportion", "Jiao-Pretis Bootstrapped Outlier Count")
  # mOutl_p[1,] <- c(x$prop.full.stat, x$boot.p.prop.test)
  # mOutl_p[2,] <- c(x$count.full.stat, x$boot.p.count.test)
  # cat("\n")
  # printCoefmat(mOutl_p, digits=6, signif.stars = TRUE,P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE) 
  
  cat("\n")
  cat("Distortion Bootstrap Results (Jiao-Pretis-Schwarz)")
  cat("\n")
  printCoefmat(matrix(c(
    c(x$L1.full,x$boot.q.L1, x$boot.p.L1),
    c(x$L2.full,x$boot.q.L2), x$boot.p.L2), nrow = 2, byrow = TRUE,
    dimnames = list(c("L1","L2"),c("Full Sample", paste0("Boot ",names(x$boot.q.L1)),"p-value"))), has.Pvalue = TRUE,P.values = TRUE)
  cat("\n")
  
  if(print.proportion){
    cat("Proportion Bootstrap Results (Jiao-Pretis)")
    cat("\n")
    print(matrix(c(
      c(x$prop.full,x$boot.q.prop),
      c(x$prop.full.stat,x$boot.q.prop.stat),
      c(x$count.full,x$boot.q.count),
      c(x$count.full.stat,x$boot.q.count.stat)), nrow = 4, byrow = TRUE,
      dimnames = list(c("Proportion","Proportion Estimate","Count","Count Estimate"),
                      c("Full Sample", paste0("Boot ",names(x$boot.q.prop))))))
    cat("\n")
  }
  

  
  cat(paste0("Cleaned Sample (always TRUE with parametric Bootstraps): ",x$args$clean.sample))
  cat("\n")
  cat(paste0("Parametric Bootstrap (residual resampling): ",x$args$parametric))
  cat("\n")
  cat(paste0("Scaled target p-val. rel. to initial p-val.: ",x$args$scale.t.pval))
}

