####################################################
## This file contains the isat-source of the gets
## package.
##
## Current version: 0.20
##
## CONTENTS:
##
## isat
## coef.isat        #extraction functions
## fitted.isat      #(some are S3 methods)
## logLik.isat
## plot.isat
## predict.isat
## print.isat
## residuals.isat
## summary.isat
## vcov.isat
##
## biascorr         #auxiliary functions:
## isattest
## isatvar
## isvarcor
## isvareffcor
## iim #make matrix of impulse indicators
## sim #make matrix of step indicators
## tim #make matrix of trend indicators
##
####################################################


####################################################
## ISAT FUNCTIONS
####################################################

##==================================================
## indicator saturation
isat <- function(y, mc=TRUE, ar=NULL, ewma=NULL, mxreg=NULL,
  iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=NULL,
  ratio.threshold=0.8, max.block.size=30, t.pval=0.001,
  wald.pval=t.pval, vcov.type=c("ordinary", "white", "newey-west"),
  do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL,
  normality.JarqueB=NULL, info.method=c("sc", "aic", "hq"), 
  user.diagnostics=NULL, user.estimator=NULL, gof.function=NULL, 
  gof.method=c("min","max"), include.gum=NULL,
  include.1cut=FALSE, include.empty=FALSE, max.paths=NULL,
  parallel.options=NULL, turbo=FALSE, tol=1e-07, LAPACK=FALSE,
  max.regs=NULL, print.searchinfo=TRUE, plot=NULL, alarm=FALSE)
{

  ##arguments:
  isat.call <- sys.call()
  vcov.type <- match.arg(vcov.type)
  info.method <- match.arg(info.method)
  gof.method <- match.arg(gof.method)
  
  ##name of regressand:
  y.name <- deparse(substitute(y))
  if( y.name[1] == "" ){ y.name <- "y" }
 
  ##determine qstat.options:
  if(is.null(ar)){
    qstat.options <- c(1,1)
  }else{
    qstat.options <- c(max(ar),1)
  }

  ##check include.gum argument:
  if(!is.null(include.gum)){
    warning("The 'include.gum' argument is ignored (temporarily deprecated in isat)")
  }
  include.gum <- TRUE

  ##make userEstArg:
  if(is.null(user.estimator)){ #default (ols):
    olsMethod <- switch(vcov.type,
      "ordinary"=3, "white"=4, "newey-west"=5)
    userEstArg <- list(name="ols", tol=tol, LAPACK=LAPACK,
      method=olsMethod)
    userEstArgArx <- NULL 
  }else{ #user-defined:
    userEstArg <- user.estimator
    userEstArgArx <- user.estimator
  }

  ##make gof.function argument:
  if(is.null(gof.function)){
    gofFunArg <- list(name="infocrit", method=info.method)
  }else{
    gofFunArg <- gof.function
  }
  
  ##max paths argument:
  if( !is.null(max.paths) && max.paths < 1 ){
    stop("'max.paths' cannot be smaller than 1")
  }

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

#OLD:
#  ##parallel.options argument:
#  if(!is.null(parallel.options)){
#    if(is.numeric(parallel.options)){
#      clusterSpec <- parallel.options
#    }
#    OScores <- detectCores()
#    if(parallel.options > OScores){
#      stop("parallel.options > number of cores/threads")
#    }
#    #to do: enable exportCluster argument?
#    #add: memory.limit()/memory.size() = max cores check?
#  }

  ##create regressors (no indicators), record info:
  mX <- regressorsMean(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg,
    return.regressand = TRUE, return.as.zoo = TRUE, na.trim = TRUE)
  y.n <- NROW(mX)
  y.index <- index(mX)
  y.index.as.char <- as.character(y.index)
  y <- coredata(mX[,1])
  #recall: y.name already defined above (in the beginning)
  if(NCOL(mX)==1){
    mX <-NULL
    mXnames <- NULL
    mXncol <- 0
    mxkeep <- NULL
  }else{
    mXnames <- colnames(mX)[-1]
    mX <- as.matrix(coredata(mX[,-1]))
    colnames(mX) <- mXnames
    mXncol <- NCOL(mX)
    mxkeep <- 1:mXncol
  }

  ##ar.LjungB argument:
  arLjungB <- NULL
  if(!is.null(ar.LjungB)){
    arLjungB <- c(NA, ar.LjungB$pval)
    if(is.null(ar.LjungB$lag)){
      arLjungB[1] <- qstat.options[1]
    }else{
      arLjungB[1] <- ar.LjungB$lag
    }
  }

  ##arch.LjungB argument:
  archLjungB <- NULL
  if(!is.null(arch.LjungB)){
    archLjungB <- c(NA, arch.LjungB$pval)
    if(is.null(arch.LjungB$lag)){
      archLjungB[1] <- qstat.options[2]
    }else{
      archLjungB[1] <- arch.LjungB$lag
    }
  }

  ##indicator saturation matrices:
  ISmatrices <- list()

  if(iis){ #impulse indicators
    mIIS <- matrix(0,y.n,y.n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste0("iis", y.index.as.char)
    ISmatrices <- c(ISmatrices,list(IIS=mIIS))
  }

  if(sis){ #step-shift indicators
    mSIS <-matrix(0,y.n,y.n) #replace by , y.n, y.n-1 ?
    loop.indx <- 1:y.n #replace by 2:y.n ?
    tmp <- function(i){ mSIS[i,1:i] <<- 1 }
    tmp <- sapply(loop.indx,tmp)
    colnames(mSIS) <- paste0("sis", y.index.as.char)
    mSIS <- mSIS[,-1]
    ISmatrices <- c(ISmatrices,list(SIS=mSIS))
  }

  if(tis){ #trend indicators
    mTIS <- matrix(0,y.n,y.n)
    v1n <- seq(1,y.n)
    loop.indx <- 1:y.n
    tmp <- function(i){
      mTIS[c(i:y.n),i] <<- v1n[1:c(y.n-i+1)]
    }
    tmp <- sapply(loop.indx,tmp)
    colnames(mTIS) <- paste0("tis", y.index.as.char)
    mTIS <- mTIS[,-1]
    ISmatrices <- c(ISmatrices,list(TIS=mTIS))
  }

  ##user-defined indicators/variables:
  ##----------------------------------

  #if uis is a matrix:
  if(!is.list(uis) && !identical(as.numeric(uis),0)){

    ##handle colnames:
    uis <- as.zoo(cbind(uis))
    uis.names <- colnames(uis)
    if(is.null(uis.names)){
      uis.names <- paste0("uisxreg", 1:NCOL(uis))
    }
    if(any(uis.names == "")){
      missing.colnames <- which(uis.names == "")
      for(i in 1:length(missing.colnames)){
       uis.names[i] <- paste0("uisxreg", missing.colnames[i])
      }
    }

    ##select sample:
    uis <- na.trim(uis, sides="both", is.na="any")
    uis.index.as.char <- as.character(index(uis))
    t1 <- which(uis.index.as.char==y.index.as.char[1])
    t2 <- which(uis.index.as.char
      == y.index.as.char[length(y.index.as.char)])
    uis <- coredata(uis)
    uis <- window(uis, start=t1, end=t2)
    uis <- cbind(coredata(as.zoo(uis)))
    colnames(uis) <- uis.names

    #check nrow(uis):
    if(nrow(uis) != y.n) stop("nrow(uis) is unequal to no. of observations")
    ISmatrices <- c(ISmatrices,list(UIS=uis))

  } #end if uis is a matrix

  ##if uis is a list of matrices:
  if(is.list(uis)){

    #check nrow(uis[[i]]):
    for(i in 1:length(uis)){
      uis[[i]] <- as.matrix(coredata(as.zoo(uis[[i]])))
      if(nrow(uis[[i]]) != y.n){
        stop(paste("nrow(uis[[",i,"]]) is unequal to no. of observations",
          sep=""))
      }
    } #end check nrow
    uis.names <- paste0("UIS", 1:length(uis))
    if(is.null(names(uis))){
      names(uis) <- uis.names
    }else{
      for(i in 1:length(uis)){
        if(names(uis)[i]==""){
          names(uis)[i] <- uis.names[i]
        }else{
          names(uis)[i] <- paste0(uis.names[i], ".", names(uis)[i])
        } #close if..else
      } #close for..loop
    }
    ISmatrices <- c(ISmatrices,uis)

    ##to do: check indices of matrix against index(y)?

  } #end if uis is a list of matrices

  ##check blocks:
  if(is.list(blocks)){
    if(length(ISmatrices)!=length(blocks)){
      stop("No. of IS matrices is unequal to length(blocks)")
    }
    blocks.is.list <- TRUE
    ISblocks <- blocks
  }else{
    blocks.is.list <- FALSE
    ISblocks <- list()
  }

  ##loop on ISmatrices:
  ##-------------------
  
  ISfinalmodels <- list()
  for(i in 1:length(ISmatrices)){

    ##blocks:
    if(!blocks.is.list){

      ncol.adj <- NCOL(ISmatrices[[i]])

      if(is.null(blocks)){
        blockratio.value <- ncol.adj/(ratio.threshold*ncol.adj - mXncol)
        blocksize.value <- ncol.adj/min(y.n*ratio.threshold, max.block.size)
        no.of.blocks <- max(2,blockratio.value,blocksize.value)
        no.of.blocks <- ceiling(no.of.blocks)
        no.of.blocks <- min(ncol.adj, no.of.blocks) #ensure blocks < NCOL
      }else{
        no.of.blocks <- blocks
      }

      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for(j in 1:no.of.blocks){
        if( blocksize*j <= ncol.adj ){
          partitions.t2[j] <- blocksize*j
        }
      }
      #check if last block contains last indicator:
      if(partitions.t2[length(partitions.t2)] < ncol.adj){
        partitions.t2 <- c(partitions.t2, ncol.adj)
      }
      blocksadj <- length(partitions.t2)
      partitions.t1 <- partitions.t2 + 1
      partitions.t1 <- c(1,partitions.t1[-blocksadj])

      tmp <- list()
      for(j in 1:blocksadj){
        tmp[[j]] <- partitions.t1[j]:partitions.t2[j]
      }
      ISblocks[[i]] <- tmp

    } #end if(!blocks.is.list)

    ##make blocks function for lapply/parLapply:
    ISblocksFun <- function(j, i, ISmatrices, ISblocks, mX,
      parallel.options, y, userEstArg, t.pval, wald.pval, do.pet,
      arLjungB, archLjungB, normality.JarqueB, user.diagnostics,
      gofFunArg, gof.method, mxkeep, include.gum, include.1cut,
      include.empty, max.paths, turbo, tol, LAPACK, max.regs,
      print.searchinfo){

      ##check if block contains 1 regressor:
      if( length(ISblocks[[i]][[j]])==1 ){
        tmp <- colnames(ISmatrices[[i]])[ ISblocks[[i]][[j]] ]
        mXis <- cbind(ISmatrices[[i]][, ISblocks[[i]][[j]] ])
        colnames(mXis) <- tmp
        mXis <- cbind(mX, mXis)
      }else{
        mXis <- cbind(mX,ISmatrices[[i]][, ISblocks[[i]][[j]] ])
      }

      ##apply dropvar:
      mXis <- dropvar(mXis, tol=tol, LAPACK=LAPACK,
        silent=print.searchinfo)

      ##print info:
      if(is.null(parallel.options)){
        if(print.searchinfo){
          message("\n", appendLF=FALSE)
          message(names(ISmatrices)[i],
            " block ", j, " of ", length(ISblocks[[i]]), ":",
            appendLF=TRUE)
          #message("\n", appendLF=FALSE)
        }
      }

      ##gum:
      getsis <- getsFun(y, mXis, untransformed.residuals=NULL,
        user.estimator=userEstArg, gum.result=NULL, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=arLjungB,
        arch.LjungB=archLjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, gof.function=gofFunArg,
        gof.method=gof.method, keep=mxkeep, include.gum=include.gum,
        include.1cut=include.1cut, include.empty=include.empty,
        max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
        max.regs=max.regs, print.searchinfo=print.searchinfo,
        alarm=FALSE)

      if(is.null(getsis$specific.spec)){
        ISspecific.models <- NULL
      }else{
        ISspecific.models <- names(getsis$specific.spec)
#For the future?:
#        ISgums[[j]] <- getsis$gum.mean
#        ISpaths[[j]] <- getsis$paths
#        ISterminals.results[[j]] <- getsis$terminals.results
      }

      ##return
      return(ISspecific.models)

    } #close ISblocksFun

    ##do gets on each block: no parallel computing
    if(is.null(parallel.options)){
      ISspecific.models <- lapply(1:length(ISblocks[[i]]),
        ISblocksFun, i, ISmatrices, ISblocks, mX, parallel.options,
        y, userEstArg, t.pval, wald.pval, do.pet, arLjungB,
        archLjungB, normality.JarqueB, user.diagnostics, gofFunArg,
        gof.method, mxkeep, include.gum, include.1cut,
        include.empty, max.paths, turbo, tol, LAPACK, max.regs,
        print.searchinfo)
    }

    ##do gets on each block: with parallel computing
    if(!is.null(parallel.options)){

      ##print info:
      if(print.searchinfo){
        message("\n", appendLF=FALSE)
        message("Preparing parallel computing...",
          appendLF=TRUE)
        message(names(ISmatrices)[i],
          " blocks to search in parallel: ", length(ISblocks[[i]]),
          appendLF=TRUE)
        message("Searching...", appendLF=TRUE)
        #message("\n", appendLF=FALSE)
      }

      blocksClust <- makeCluster(clusterSpec, outfile="") #make cluster
      clusterExport(blocksClust, clusterVarlist,
        envir=.GlobalEnv) #idea for the future?: envir=clusterEnvir
#OLD:
#      clusterExport(blocksClust,
#        c("dropvar", "getsFun", "ols", "infocrit", "diagnostics"),
#        envir=.GlobalEnv)
      ISspecific.models <- parLapply(blocksClust,
        1:length(ISblocks[[i]]), ISblocksFun, i, ISmatrices,
        ISblocks, mX, parallel.options, y, userEstArg, t.pval,
        wald.pval, do.pet, arLjungB, archLjungB, normality.JarqueB,
        user.diagnostics, gofFunArg, gof.method, mxkeep,
        include.gum, include.1cut, include.empty, max.paths, turbo,
        tol, LAPACK, max.regs, print.searchinfo)
      stopCluster(blocksClust)

    } #end if..

    ##print info:
    if(print.searchinfo){
      message("\n", appendLF=FALSE)
      message("GETS of union of retained ",
        names(ISmatrices)[i], " variables... ",
        appendLF=TRUE)
    }

    ##if no indicators retained from the blocks:
    if(length(ISspecific.models) == 0){
      isNames <- NULL
      ISfinalmodels[[i]] <- NULL
    }

    ##when indicators/variables(uis) retained from the blocks:
    if(length(ISspecific.models) > 0){

      isNames <- NULL

      #which indicators/variables(uis) retained?:
      for(j in 1:length(ISspecific.models)){
        #check if mean is non-empty:
        if(!is.null(ISspecific.models[[j]])){
          isNames <- union(isNames, ISspecific.models[[j]])
        }
      } #end for(j) loop
      isNames <- setdiff(isNames, mXnames)

      #redo gets with union of retained indicators:
      if(length(isNames) == 0){
        ISfinalmodels[[i]] <- mXnames
      }else{
        mXisNames <- c(mXnames, isNames)
        mXis <- cbind(mX,ISmatrices[[i]][,isNames])
        colnames(mXis) <- mXisNames
        mXis <- dropvar(mXis, tol=tol, LAPACK=LAPACK,
          silent=print.searchinfo)

        getsis <- getsFun(y, mXis, untransformed.residuals=NULL,
          user.estimator=userEstArg, gum.result=NULL, t.pval=t.pval,
          wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=arLjungB,
          arch.LjungB=archLjungB, normality.JarqueB=normality.JarqueB,
          user.diagnostics=user.diagnostics, gof.function=gofFunArg,
          gof.method=gof.method, keep=mxkeep, include.gum=include.gum,
          include.1cut=include.1cut, include.empty=include.empty,
          max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
          max.regs=max.regs, print.searchinfo=print.searchinfo,
          alarm=FALSE)

        ISfinalmodels[[i]] <- names(getsis$specific.spec)
      }

    } #end if(length(ISspecific.models > 0)

  } #end for(i) loop (on ISmatrices)

  ##add names to ISblocks:
  names(ISblocks) <- names(ISmatrices)


  ##gets of union of all variables:
  ##-------------------------------
  
  ##some info:
  if(print.searchinfo){
    message("\n", appendLF=FALSE)
    message("GETS of union of ALL retained variables...",
      appendLF=TRUE)
    #message("\n", appendLF=FALSE)
  }

  ##if final models estimated:
  if(length(ISfinalmodels)>0){

    mIS <- NULL #becomes a matrix

    #which indicators were retained?
    for(i in 1:length(ISfinalmodels)){
      isNames <- NULL
      #check if non-empty:
      if(!is.null(ISfinalmodels[[i]])){
        isNames <- setdiff(ISfinalmodels[[i]], mXnames)
      }
      if(length(isNames)>0){
        tmp <- cbind(ISmatrices[[i]][, isNames ])
        colnames(tmp) <- isNames
        mIS <- cbind(mIS, tmp)
      }
    } #end for loop

    mXis <- dropvar(cbind(mX,mIS), tol=tol, LAPACK=LAPACK,
      silent=print.searchinfo)

  } #end if(length(ISfinalmodels)>0)


  ##if no final models estimated:
  if(length(ISfinalmodels)==0){
    ISfinalmodels <- NULL
    if(is.null(mX)){ mXis <- NULL }else{
      mXis <- cbind(mX)
      colnames(mXis) <- mXnames
    }
  }


  ##make return object:
  ##-------------------

  ##do final gets:
  getsis <- getsFun(y, mXis, untransformed.residuals=NULL,
    user.estimator=userEstArg, gum.result=NULL, t.pval=t.pval,
    wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=arLjungB,
    arch.LjungB=archLjungB, normality.JarqueB=normality.JarqueB,
    user.diagnostics=user.diagnostics, gof.function=gofFunArg,
    gof.method=gof.method, keep=mxkeep, include.gum=include.gum,
    include.1cut=include.1cut, include.empty=include.empty,
    max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
    max.regs=max.regs, print.searchinfo=print.searchinfo,
    alarm=FALSE)
  ##messages from final gets:
  if( print.searchinfo && !is.null(getsis$messages)){ message(getsis$messages) }
       
  ##estimate final model:
  y <- zoo(y, order.by=y.index)
  if(is.null(getsis$specific.spec)){
    mXis <- NULL
  }else{
    mXisNames <- colnames(mXis)[getsis$specific.spec]
    mXis <- cbind(mXis[,getsis$specific.spec])
    colnames(mXis) <- mXisNames
    mXis <- zoo(mXis, order.by=y.index)
  }
  if(is.null(normality.JarqueB)){
    normalityArg <- FALSE
  }else{
    normalityArg <- as.numeric(normality.JarqueB)  
  }
  mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
    qstat.options=qstat.options, normality.JarqueB=normalityArg,
    user.estimator=userEstArgArx, user.diagnostics=user.diagnostics,
    tol=tol, LAPACK=LAPACK, plot=FALSE)
  mod$call <- NULL
   
  ##complete the returned object (result):
  ISnames <- setdiff(getsis$aux$mXnames, mXnames) #names of retained impulses
  if(length(ISnames)==0){ ISnames <- NULL }
  colnames(mod$aux$mX) <- mod$aux$mXnames #needed for predict.isat?
  getsis$gets.type <- "isat"
  getsis$call <- isat.call
  getsis <- c(list(ISfinalmodels=ISfinalmodels,
    ISnames=ISnames), getsis, mod)
  getsis$aux$t.pval <- t.pval #needed for biascorr
  class(getsis) <- "isat"
  if(alarm){ alarm() }
  if( is.null(plot) ){ #determine whether to plot or not
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.isat(getsis, coef.path=TRUE) }
  return(getsis)

} #close isat function

##==================================================
coef.isat <- function(object, ...)
{
  result <- object$coefficients
  if(!is.null(result)){ names(result) <- names(object$specific.spec) }
  return(result)
} #close coef.isat


##==================================================
## fitted values
fitted.isat <- function(object, ...)
{
  result <- object$mean.fit
  if(is.null(result)){ result <- object$fit }
  return(result)
} #end fitted.isat

##==================================================
logLik.isat <- function(object, ...)
{
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
} #end logLik.isat

##==================================================
## plot isat object
plot.isat <- function(x, col=c("red","blue"),
  lty=c("solid","solid"), lwd=c(1,1), coef.path=TRUE, ...)
{

#  ##check if mean quation:
#  if( is.null(x$mean.results) ){
#    cat("No mean equation to plot\n")
#  }

  ##if fitted mean:
  if(!is.null(x$mean.fit)){

    ##check line width:
    if(length(lwd)==1){
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd=rep(lwd,2)
    }else if (length(lwd)>2){
      print("lwd needs two arguments, but more provided. First two used.")
      lwd=lwd[1:2]
    }

    ##check line type:
    if(length(lty)==1){
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty=rep(lty,2)
    }else if (length(lwd)>2){
      print("lty needs two arguments, but more provided. First two used.")
      lty=lty[1:2]
    }

    ##check colour:
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
    }

    ##get fitted and actual values, and dependent variable name
    fitted <- x$mean.fit
    actual <- zoo(x$aux$y, order.by=x$aux$y.index)
    residuals <- x$std.residuals
    actual.name <- x$aux$y.name

    ##get current par-values, set new ones:
    def.par <- par(no.readonly=TRUE)
    par(mar=c(2,2,0.5,0.5))

    ##if isat or coef-path:
    if( (x$gets.type=="isat" || coef.path==TRUE) && length(x$ISnames)!=0 ){
      par(mfrow=c(3,1))
      is.x <- cbind(x$aux$mX[,x$aux$mXnames %in% x$ISnames])
      is.coef.ests <- coef.isat(x)[x$ISnames]
      coef.path.0 <- zoo(is.x %*% is.coef.ests, order.by=x$aux$y.index)
    }else{
      par(mfrow=c(2,1))
    }

    ##comment?
    if(is.regular(actual)) {
      plot(actual, main = "",ylim=range(min(actual,fitted,na.rm=TRUE),max(actual,fitted,na.rm=TRUE)),
           type="l",ylab="",xlab="",col=col[2])
    }else{
      plot(as.Date(index(actual)),coredata(actual), main = "",ylim=range(min(actual,fitted,na.rm=TRUE),max(actual,fitted,na.rm=TRUE)),
           type="l",ylab="",xlab="",col=col[2])
    }
    if(is.regular(fitted)) {
      lines(fitted,col=col[1])
    }else{
      lines(as.Date(index(fitted)),coredata(fitted),col=col[1])
    }

    legend("topleft",lty=lty,lwd=lwd,ncol=2,col=col[c(2,1)],legend=c(actual.name,"fitted"),bty="n")
    if(is.regular(residuals)) {
      plot(residuals,type="h",col=col[1])
    }
    else{
      plot(as.Date(index(residuals)),coredata(residuals),type="h",col=col[1])
    }

    abline(0,0)
    legend("topleft",lty=1,col=col[1],legend="standardised residuals",bty="n")
#    legend("topleft",lty=1,col=col[1],legend=c(paste(actual.name,"standardised residuals",sep=": ")),bty="n")

    ##coefficient path
    if( (x$gets.type=="isat" | coef.path==TRUE) & length(x$ISnames)!=0 ) {
      ## we only get standard error bars if TIS *not* run


      ###if tis is there and it is not null, then don't plot, else


      if(!is.null(as.list(x$call)$tis) && as.list(x$call)$tis==TRUE){
        message("\n", appendLF=FALSE)
        message("NB: Because TIS selected, coefficient standard errors invalid hence not plotted",
          appendLF=TRUE)
#OLD:
#        cat("\nNB: Because TIS selected, coefficient standard errors invalid hence not plotted\n", sep="")
        ylim.values <- range(coef.path.0)
        if(is.regular(coef.path.0)) {
          ylim.values <- range(coef.path.0)
          plot(coef.path.0,type="l",col=col[1],ylim=ylim.values)
        }else{
          plot(as.Date(index(coef.path.0)),coredata(coef.path.0),type="l",col=col[1],ylim=ylim.values)
        }
      }   else {

        coef.path.v <- isatvar(x)
        if(is.regular(coef.path.0)) {
          ylim.values <- range(min(coef.path.0-qt(0.975, NROW(coef.path.0))*coef.path.v$const.se),
                               max(coef.path.0+qt(0.975, NROW(coef.path.0))*coef.path.v$const.se))
          plot(coef.path.0,type="l",col=col[1],
               ylim=ylim.values)
          lines(coef.path.0+qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)
          lines(coef.path.0-qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)
        }else{
          plot(as.Date(index(coef.path.0)),coredata(coef.path.0),type="l",col=col[1],ylim=ylim.values)
          lines(as.Date(index(coef.path.0)),coredata(coef.path.0)+qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)
          lines(as.Date(index(coef.path.0)),coredata(coef.path.0)-qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)

        }

      }

      abline(0,0,lty=3)
      legend("topleft",lty=1,col=col[1],legend=c(paste(actual.name,"Coefficient Path",sep=": ")),bty="n")
    }

    #return to old par-values:
    par(def.par)
  }
} #plot.isat closed

##==================================================
## forecast up to n.ahead
predict.isat <- function(object, n.ahead=12, newmxreg=NULL,
  newindex=NULL, n.sim=2000, probs=NULL, ci.levels=NULL,
  quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
  plot.options=list(), ...)

{
  ## 1 arguments of mean-equation
  ## 2 make plot.options argument
  ## 3 pass arguments on to predict.arx
  ## 4 return forecasts

  ##create new object to add stuff to in order to use predict.arx()
  objectNew <- object

  ##-----------------------------------
  ## 1 arguments of mean-equation:
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
    gumTerms <- eval(object$call$ar)
#OLD:
#    gumTerms <- eval(object$aux$call.gum$ar)
    gumNamesAr <- paste0("ar", gumTerms)
    whichRetained <- which( gumNamesAr %in% names(coefsMean) )
    if( length(whichRetained)==0 ){
      objectNew$call$ar <- NULL
    }else{
      objectNew$call$ar <- gumTerms[ whichRetained ]
      indxCounter <- indxCounter + length(whichRetained)
    }

    ##ewma argument:
    gumTerms <- eval(object$call$ewma)
#OLD:
#    gumTerms <- eval(object$aux$call.gum$ewma)
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
#      objectNew$call$mxreg <- mxreg
    }

  } #end if( length(coefsMean)>0 )

  ##here: introduce the automated detection of indicators and how they
  ##should modify newmxreg?

  ##----------------------------------
  ## 2 make plot.options argument:
  ##----------------------------------
  
  if(is.null(plot.options$start.at.origin)){
    plot.options$start.at.origin <- FALSE
  }
  if(is.null(plot.options$line.at.origin)){
    plot.options$line.at.origin <- TRUE
  }
  if(is.null(plot.options$fitted)){
    plot.options$fitted <- TRUE
  }
  
  
  ##-------------------------------------
  ## 3 pass arguments on to predict.arx:
  ##-------------------------------------

  innov <- rnorm(n.ahead*n.sim) #force normal errors
  result <- predict.arx(objectNew, spec="mean", n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=NULL, newindex=newindex,
    n.sim=n.sim, innov=innov, probs=probs, ci.levels=ci.levels,
    quantile.type=quantile.type, return=return, verbose=verbose,
    plot=plot, plot.options=plot.options)

  ##---------------------
  ## 4 return forecasts:
  ##---------------------

  if(return){ return(result) }

} #close predict.isat

##==================================================
# ## print isat results
print.isat <- function(x, ...)
{

  ##messages from final gets:
  #if(!is.null(x$messages)){ message(x$messages) }
  
  ##header:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Dependent var.:", x$aux$y.name, "\n")
  cat("Method: Ordinary Least Squares (OLS)\n")
  cat("Variance-Covariance:", switch(x$aux$vcov.type,
                                     ordinary = "Ordinary", white = "White (1980)",
                                     "newey-west" = "Newey and West (1987)"), "\n")
  
  ##header - sample info:
  cat("No. of observations (mean eq.):", x$aux$y.n, "\n")
  tmp <- zoo(x$aux$y, order.by=x$aux$y.index)
  
  indexTrimmed <- index(na.trim(tmp))
  isRegular <- is.regular(tmp, strict=TRUE)
  isCyclical <- frequency(tmp) > 1
  if(isRegular && isCyclical){
    cycleTrimmed <- cycle(na.trim(tmp))
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

####### START the part commented out 17 July 2019 by G-man:  
#
#  ##gum:
#  if(specType=="mean"){
#    cat("\n")
#    cat("GUM mean equation:\n")
#    cat("\n")
#    printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2),
#                 signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
#  }
#  if(!is.null(x$gum.variance)){
#    cat("\n")
#    cat("GUM log-variance equation:\n")
#    cat("\n")
#    printCoefmat(x$gum.variance, signif.stars=FALSE)
#  }
#  cat("\n")
#  cat("Diagnostics and fit:\n")
#  cat("\n")
#  printCoefmat(x$gum.diagnostics, dig.tst=0, tst.ind=2,
#               signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
#  
#
#  ##paths:
#  cat("\n")
#  cat("Paths searched: \n")
#  cat("\n")
#  if(is.null(x$paths)){
#    print(NULL)
#  }else{
#    for(i in 1:length(x$paths)){
#      cat("path",i,":",x$paths[[i]],"\n")
#    }
#  } #end if(is.null(x$paths))
#  
#  ##terminal models and results:
#  if(!is.null(x$terminals)){
#    cat("\n")
#    cat("Terminal models: \n")
#    cat("\n")
#    for(i in 1:length(x$terminals)){
#      cat("spec",i,":",x$terminals[[i]],"\n")
#    }
#  }
#  if(!is.null(x$terminals.results)){
#    cat("\n")
#    printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
#                 signif.stars=FALSE)
#  }
#
####### END the part commented out 17 July 2019 by G-man  
  
  ##specific model:
  if(!is.null(x$specific.spec)){
    cat("\n")
    cat("SPECIFIC mean equation:\n")
    cat("\n")
    if(!is.null(x$mean.results)){
      print(x$mean.results)
      #OLD: DOES NOT WORK IN A PREDICTABLE WAY!
      #      printCoefmat(x$mean.results, signif.stars=FALSE,
      #        P.values=FALSE, has.Pvalues=FALSE)
    }
    if(x$specific.spec[1]==0){
      cat("empty\n")
    }
    ##in the future: use estimate.specific=FALSE more directly?
    if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
      cat("Not estimated\n")
    }
  }
  
  ##diagnostics and fit:
  if(!is.null(x$diagnostics)){
    
    #fit-measures:
    mGOF <- matrix(NA, 3, 1)
    rownames(mGOF) <- c("SE of regression", "R-squared",
                        paste0("Log-lik.(n=", x$n, ")"))
    colnames(mGOF) <- ""
    mGOF[1,1] <- sigma.isat(x) #OLD: sqrt(x$sigma2)
    mGOF[2,1] <- rsquared(x) #OLD: x$specific.diagnostics[4,1]
    mGOF[3,1] <- as.numeric(logLik.isat(x)) #OLD: x$logl
    #mGOF[4,1] <- outliertest(x)$#x$logl #OLD: as.numeric(logLik.arx(x))
    
    cat("\n")
    cat("Diagnostics and fit:\n")
    cat("\n")
    printCoefmat(x$diagnostics, dig.tst=0, tst.ind=2,
                 signif.stars=FALSE)
    if(!is.null(x$call$iis)){
      if (x$call$iis==TRUE){
        outltest <- outliertest(x)
        mOutl <- matrix(NA, 2, 2)
        colnames(mOutl) <- c("Stat.", "p-value")
        rownames(mOutl) <- c("Jiao-Pretis Prop.", "Jiao-Pretis Count")
        mOutl[1,] <- c(outltest$prop$statistic, outltest$prop$p.value)
        mOutl[2,] <- c(outltest$count$statistic, outltest$count$p.value)
        cat("\n")
        printCoefmat(mOutl, digits=6, signif.stars = FALSE) 
        #cat("\n")
      }
    }
    printCoefmat(mGOF, digits=6, signif.stars=FALSE)
    
  }
    
} #end print.isat

#OLD:
##==================================================
# ## print isat results
# print.isat <- function(x, ...)
# {
#   ##specification type:
#   specType <- "mean"
# 
#   ##header:
#   cat("\n")
#   cat("Date:", x$date, "\n")
#   cat("Dependent var.:", x$aux$y.name, "\n")
#   cat("Method: Ordinary Least Squares (OLS)\n")
#   cat("Variance-Covariance:", switch(x$aux$vcov.type,
#     ordinary = "Ordinary", white = "White (1980)",
#     "newey-west" = "Newey and West (1987)"), "\n")
# 
#   ##header - sample info:
#   cat("No. of observations (mean eq.):", x$aux$y.n, "\n")
#   tmp <- zoo(x$aux$y, order.by=x$aux$y.index)
# 
#   indexTrimmed <- index(na.trim(tmp))
#   isRegular <- is.regular(tmp, strict=TRUE)
#   isCyclical <- frequency(tmp) > 1
#   if(isRegular && isCyclical){
#     cycleTrimmed <- cycle(na.trim(tmp))
#     startYear <- floor(as.numeric(indexTrimmed[1]))
#     startAsChar <- paste(startYear,
#       "(", cycleTrimmed[1], ")", sep="")
#     endYear <- floor(as.numeric(indexTrimmed[length(indexTrimmed)]))
#     endAsChar <- paste(endYear,
#       "(", cycleTrimmed[length(indexTrimmed)], ")", sep="")
#   }else{
#     startAsChar <- as.character(indexTrimmed[1])
#     endAsChar <- as.character(indexTrimmed[length(indexTrimmed)])
#   }
#   cat("Sample:", startAsChar, "to", endAsChar, "\n")
# 
#   ##gum:
#   if(specType=="mean"){
#     cat("\n")
#     cat("GUM mean equation:\n")
#     cat("\n")
#     printCoefmat(x$gum.mean, dig.tst=0, tst.ind=c(1,2),
#       signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
#   }
#   if(!is.null(x$gum.variance)){
#     cat("\n")
#     cat("GUM log-variance equation:\n")
#     cat("\n")
#     printCoefmat(x$gum.variance, signif.stars=FALSE)
#   }
#   cat("\n")
#   cat("Diagnostics and fit:\n")
#   cat("\n")
#   printCoefmat(x$gum.diagnostics, dig.tst=0, tst.ind=2,
#     signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
# 
#   ##paths:
#   cat("\n")
#   cat("Paths searched: \n")
#   cat("\n")
#   if(is.null(x$paths)){
#     print(NULL)
#   }else{
#     for(i in 1:length(x$paths)){
#       cat("path",i,":",x$paths[[i]],"\n")
#     }
#   } #end if(is.null(x$paths))
# 
#   ##terminal models and results:
#   if(!is.null(x$terminals)){
#     cat("\n")
#     cat("Terminal models: \n")
#     cat("\n")
#     for(i in 1:length(x$terminals)){
#       cat("spec",i,":",x$terminals[[i]],"\n")
#     }
#   }
#   if(!is.null(x$terminals.results)){
#     cat("\n")
#     printCoefmat(x$terminals.results, dig.tst=0, tst.ind=c(3,4),
#       signif.stars=FALSE)
#   }
# 
#   ##specific model:
#   if(specType=="mean" && !is.null(x$specific.spec)){
#     cat("\n")
#     cat("SPECIFIC mean equation:\n")
#     cat("\n")
#     if(!is.null(x$mean.results)){
#       print(x$mean.results)
# #OLD: DOES NOT WORK IN A PREDICTABLE WAY!
# #      printCoefmat(x$mean.results, signif.stars=FALSE,
# #        P.values=FALSE, has.Pvalues=FALSE)
#     }
#     if(x$specific.spec[1]==0){
#       cat("empty\n")
#     }
# ##in the future: use estimate.specific=FALSE more directly?
#     if(x$specific.spec[1]!=0 && is.null(x$mean.results)){
#       cat("Not estimated\n")
#     }
#   }
# 
#   ##diagnostics and fit:
#   if(!is.null(x$specific.diagnostics)){
# 
#     #fit-measures:
#     mGOF <- matrix(NA, 3, 1)
#     rownames(mGOF) <- c("SE of regression", "R-squared",
#       paste("Log-lik.(n=", length(na.trim(x$std.residuals)), ")", sep=""))
#     colnames(mGOF) <- ""
#     mGOF[1,1] <- sigma.isat(x) #OLD: sqrt( RSS/(nobs-DFs) )
#     mGOF[2,1] <- rsquared(x) #OLD: x$specific.diagnostics[4,1]
#     mGOF[3,1] <- as.numeric(logLik.arx(x))
# 
#     cat("\n")
#     cat("Diagnostics and fit:\n")
#     cat("\n")
#     printCoefmat(x$specific.diagnostics, dig.tst=0, tst.ind=2,
#       signif.stars=FALSE)
#     printCoefmat(mGOF, digits=6, signif.stars=FALSE)
# 
#   }
# 
#   ##messages:
#   if(!is.null(x$messages)){
#     message("\n", appendLF=FALSE)
#     message(x$messages)
#   }
# 
# } #end print.isat

##==================================================
residuals.isat <- function(object, std=FALSE, ...)
{
  if(is.null(object)){
    result <- NULL
  }else{
    if(std){
      result <- object$residuals/sigma.isat(object)
    }else{
      result <- object$residuals
    }
  }
  return(result)
} #end residuals.isat

##==================================================
## SE of regression
sigma.isat <- function(object, ...)
{
  if(is.null(object$residuals)){
    result <- NULL
  }else{
    RSS <- sum(object$residuals^2)
    result <- sqrt(RSS/(object$n - object$k))
  }
  return(result)
} #close sigma.isat

###==================================================
### summarise output
summary.isat <- function(object, ...)
{
  summary.default(object)
} #end summary.isat

##==================================================
vcov.isat <- function(object, ...)
{
  result <- object$vcov #also works if object$vcov.mean exists
  if(!is.null(result)){
    if(is.null(colnames(result))){
      colnames(result) <- names(object$specific.spec)
    }
    if(is.null(rownames(result))){
      rownames(result) <- names(object$specific.spec)
    }
  }
  return(result)  
} #end vcov.isat

##==================================================
## needs orthogonal variables, can be applied directly to the coefficient path though in an isat only model
### purposefully not entirely usefriendly - needs many arguments, but it should be an expert user function
biascorr <- function(b, b.se, p.alpha, T){

  c_alpha <- abs(qt((1-(1-p.alpha))/2, T))
  bt <- b/b.se

  dr <- (dnorm(c_alpha-bt)-dnorm(-c_alpha-bt))/ (1-pnorm(c_alpha-bt) + pnorm(-c_alpha-bt))
  dtbar <- bt - dr
  drbar <- (dnorm(c_alpha-dtbar)-dnorm(-c_alpha-dtbar))/ (1-pnorm(c_alpha-dtbar) + pnorm(-c_alpha-dtbar))

  #only correct if significant

  b_1step <- b
  b_1step[abs(bt)>c_alpha] <- b*(1-(dr/bt))[abs(bt)>c_alpha]

  b_2step <- b
  b_2step[abs(bt)>c_alpha] <- b*(1-(drbar/bt))[abs(bt)>c_alpha]

  b_corr <- cbind(b, b_1step, b_2step)
  colnames(b_corr) <- c("beta", "beta.1step", "beta.2step")

  return(b_corr)
}

##==================================================
####### isat.test: forecast bias test as based on the working paper with a plot
### loads an isat object, and then conducts the bias analysis on it


### loads an isat object, and then conducts the bias analysis on it
#### Now with added efficiency and consistency correction
#### also on coefficient path of other break variables


isattest <- function(x, hnull=0, lr=FALSE, ci.pval = 0.99,
  plot=NULL, plot.turn = FALSE, conscorr=FALSE, effcorr=FALSE,
  mcor = 1, biascorr=FALSE, mxfull = NULL, mxbreak=NULL)
{
  trend.incl <- FALSE
  if(!is.null(as.list(x$call)$tis)){
    if (as.list(x$call)$tis) {
      stop("isat.test currently not implemented for trend-indicator saturation")
      trend.incl <- TRUE
    }
  }

  arcall <- as.list(x$call)$ar
  x.var <- isatvar(x,lr=lr, conscorr=conscorr, effcorr=effcorr, mcor = mcor, mxfull = mxfull, mxbreak=mxbreak)

  #misy1.var


  if (biascorr==TRUE){

    if (!is.null(as.list(x$call)$mxreg) | !is.null(arcall) | trend.incl){

      biascorr <- FALSE
      message("Bias Correction not applicable in isat regression with additional non-step covariates. Has been set to FALSE.")
    }
  }


  T <- dim(x$aux$mX)[1]
  N <- dim(x$aux$mX)[2]

  crit <- abs(qt((1-ci.pval)/2, T-N))
  bias.low <- matrix(0, T, 1)
  bias.high <- matrix(0, T, 1)

  ci.low <- matrix(0, T, 1)
  ci.high <- matrix(0, T, 1)
  x.mean <- matrix(0, T, 1)

  if (lr == TRUE & !is.null(arcall))
  {
    x.is.lr <- x.var$lr.path
    x.is.const <- x.var$const.path

  } else {

    x.is.lr <- NA

    if (biascorr){

      xbias <-biascorr(b=x.var$const.path, b.se=x.var$const.se, p.alpha = x$aux$t.pval, T=length(x.var$const.path))
      x.is.const <- xbias$beta.2step

    } else {
      x.is.const <- x.var$const.path

    }


  }


  if (lr == TRUE & !is.null(arcall))
  {

    ci.low <- x.var$lr.path-crit*x.var$lr.se
    ci.high <- x.var$lr.path+crit*x.var$lr.se

    bias.low[which((ci.low) > hnull)] <- 1
    bias.high[which((ci.high) < hnull)] <- 1

    bias.low <- bias.low*(x.var$lr.path-hnull)
    bias.high <- bias.high*(x.var$lr.path-hnull)

    x.mean <- x.var$lr.path

  } else {

    ci.low <- x.is.const-crit*x.var$const.se
    ci.high <- x.is.const+crit*x.var$const.se

    bias.low[which((ci.low) > hnull)] <- 1
    bias.high[which((ci.high) < hnull)] <- 1


    bias.low <- bias.low*(x.is.const-hnull)
    bias.high <- bias.high*(x.is.const-hnull)

    x.mean <- x.is.const

  }

  ###determining the turning points
  time <- x$aux$y.index
  bias.sum.ar <- bias.low+bias.high
  lr.path.d <- diff(bias.sum.ar)


  if(all(lr.path.d==0)){
    plot.turn <- FALSE
    turn.ar <- NULL
  } else {
    turn.ar <- time[which(lr.path.d != 0)]+1
  }

  turn.ar.y <- bias.sum.ar[which(lr.path.d != 0)]

  turn.x.lab <- turn.ar
  turn.x <- turn.ar

  fitted <- x$mean.fit
  actual <- zoo(x$aux$y, order.by=x$aux$y.index)

  ylabel_a <- "Coefficient"
  ylabel_b <- "Bias"



  Ylim_main <- c(min(actual, na.rm=TRUE)*1.2,max(actual, na.rm=TRUE)*1.2)
  Ylim_bias <- c(min(bias.high, na.rm=TRUE)*1.2,max(bias.low, na.rm=TRUE)*1.2)

  ##plot argument:
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }

  if (plot){
    par(mfrow=c(2,1), mar = c(2, 4,1,3))
    plot(time, x.mean, ylim=Ylim_main, col="red", main=NULL, xlab=NA, ylab=ylabel_a, sub=NA, type="l")

    lines(ci.low, col="red", lty=2)
    lines(ci.high, col="red", lty=2)

    if (is.null(mxbreak))
    {
      lines(actual, col="blue")
    }
    abline(a =hnull, b=0, col="black", lty=3, lwd=2)
    plot(time, bias.low, type="h", col="red", ylim=Ylim_bias, main=NULL, xlab=NA, ylab=ylabel_b, sub=NA)
    lines(bias.high, type="h", col="red")

    if ( plot.turn ){
      text(turn.x.lab, y=turn.ar.y, x=turn.x, pos=4, offset=-0.5, cex=0.8)
    }
  }


  if (lr==TRUE & !is.null(arcall))
  {
    mean.var <- cbind(ci.low, ci.high, bias.low,  bias.high)
    colnames(mean.var) <- c("ci.low", "ci.high", "bias.high",  "bias.low")


  } else {
    mean.var <- cbind( ci.low, ci.high, bias.low,  bias.high)
    colnames(mean.var) <- c("ci.low", "ci.high", "bias.high",  "bias.low")

  }


  return(mean.var)

} #isattest function closed



##==================================================
####### isatvar: function to extract the variance of the coefficient path
#takes an isat gets object, returns:
# - the coefficient path (both relative to the full sample coefficient and the constant path itself)
# - the variance and standard errors of the coefficient path
# - if lr is specified and the isat object has AR terms, then also returns the LR coefficient path and its variance and standard errors
#input: "x" (an isat results object)

#new arguments:
#-conscorr applies the consistency correction from Johansen and Nielsen 2016
#-effcorr applies the efficiency correction from Johansen and Nielsen 2016
#-mxfull allows a full-sample break variable to be specified which is selected over using uis, so it can construct the MIS variance

isatvar <- function(x, lr=FALSE, conscorr=FALSE, effcorr=FALSE, mcor = 1,  mxfull = NULL, mxbreak=NULL)
{



  if (lr == TRUE && !is.null(mxfull)){
    message("Warning: LR currently not defined with user-specified break variables. LR set to FALSE")
    lr <- FALSE
  }



  if(!is.null(as.list(x$call)$uis) && is.null(mxfull)){

    if (!is.null(as.list(x$call)$iis) || !is.null(as.list(x$call)$sis) ){
      message("Warning: uis specified but no mxfull variable given. Using mconst instead.")
      mxfull <- "mconst"
    } else {
      stop("uis specified but no mxfull variable given")
    }


  }


  if (is.null(mxfull)){   ##if no full variable is specified, then the constant is used
    mxfull <- "mconst"
  }


  if (!is.null(x$mean.fit)){

    if (!is.null(mxbreak))
    {

      ISnames <- x$ISnames[grep(mxbreak, x$ISnames)]
      ####################################
    } else {

      mxbreak_iis <- "iis"
      mxbreak_sis <- "sis"

      ISnames <- c(x$ISnames[grep(mxbreak_iis, x$ISnames)], x$ISnames[grep(mxbreak_sis, x$ISnames)])
    }


    if (!is.null(ISnames)){


      var.rel <- c( which(substr(x$aux$mXnames,1,6) %in% mxfull), which(x$aux$mXnames %in% ISnames))


      if (conscorr==TRUE){
        x$vcov.mean <- x$vcov.mean*as.numeric(isvarcor(x$aux$t.pval,1)[2]^2)
      }

      if (effcorr==TRUE){

        if (!is.null(x$keep)) {
          x$vcov.mean[x$keep, x$keep] <- x$vcov.mean[x$keep, x$keep] * as.numeric(isvareffcor(x$aux$t.pval, 1, mcor)[2]^2)
        }
      }

      #coefficient path function
      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% ISnames])
      is.coef.ests <- coef.isat(x)[ISnames]

      if(is.null(mxfull)){
        coef.path <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
      } else {
        is.x[is.x!=0] <- 1 #the break variable is not necessarily an indicator variable
        coef.path <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
      }




      colnames(coef.path) <- "coef.path"

      const <- coef.isat(x)[mxfull]
      const.path <- coef.path + const
      colnames(const.path) <- "const.path"

    } else {
      var.rel <- which(substr(x$aux$mXnames,1,6) %in% mxfull)

      #coefficient path function
      const <- coef.isat(x)[mxfull]

      coef.path <- zoo(0, order.by = x$aux$y.index)

      const.path <- coef.path + const

    }

    vcov.rel <- x$vcov.mean[var.rel,var.rel]


    dim.var <- NCOL(vcov.rel)
    const.var <- matrix(NA, dim.var, 1 )

    #construct a matrix to multiply by the variances
    dim.in <- NCOL(x$aux$mX[,var.rel])

    const.mat <- matrix(NA, NROW(x$aux$mX[,var.rel]), dim.in)

    if(is.null(mxfull)){
      indic.mat <- x$aux$mX[,var.rel]

    } else {
      indic.mat <- x$aux$mX[,var.rel]
      indic.mat[indic.mat!=0] <- 1

    }



    if (dim.var > 1) #if there are indicators retained
    {

      #order the breaks in correct order so the covariance matrix and relevance order is correct

      if(is.null(mxfull)){
        order.mat <- apply(indic.mat[,], 2, function(x) min(which(x==1)))  #find where each indicator is first one, for sorting

      } else { #if using a custom break variable
        order.mat <- apply(indic.mat[,], 2, function(x) min(which(x!=0)))
      }




      indic.mat <- indic.mat[,order(order.mat)]  #order the indicators
      vcov.rel <- vcov.rel[order(order.mat),order(order.mat)] #order the covariance matrix

      for (i in 1:dim.var)
      {
        const.var[i] <- sum(vcov.rel[1:i,1:i])   #sum over the expanding variance covariance matrix to get the variance of the sums of coefficients
      }

      const.mat <- indic.mat

      for (j in 2:(dim.in))
      {

        const.mat[which(rowSums(as.matrix(const.mat[,(j):(dim.in)]))>0),j-1] <- 0     #puts zeros in the appropriate places to make sure the correct s.e. is applied for each point in time

      }

      ind.var.mat <- const.mat %*% const.var   #the variance as it applies to each subsection

    } else { #just the constant remains

      const.var <- vcov.rel
      const.mat <- indic.mat
      ind.var.mat <- const.mat * const.var

    }

    ind.se.mat <- sqrt(ind.var.mat)   #the standard errors of the coefficient as it changes with the SIS breaks

    ####the coefficient path of the LR mean

    if(lr){
      if (!is.null(as.list(x$call)$ar)){

        if (!is.null(x$mean.fit)){

          vcov.rel.tot <- x$vcov.mean



          coef.rel <- coef.isat(x)

          ###coefficient path of LR mean

          arcall <- as.list(x$call)$ar
          arnames <- paste("ar",eval(as.expression(arcall)), sep="")

          ar.coefs <- coef.isat(x)[arnames]
          ar.sum <- sum(ar.coefs)
          lr.path <- const.path/(1-ar.sum)

          ###variance of the coefficient path

          ar.var <- vcov.rel.tot[arnames,arnames] #variance of ar() terms
          armu.cov <- vcov.rel.tot[c(mxfull, ISnames), arnames] #covariance of ar(1) and const + sis


          if (!is.null(ISnames)) {
            var.rel <- c(which(substr(x$aux$mXnames, 1, 6) %in%
                                 mxfull), which(x$aux$mXnames %in% ISnames))
          } else {
            var.rel <- which(substr(x$aux$mXnames, 1, 6) %in%
                               mxfull)
          }

          vcov.rel <- x$vcov.mean[var.rel, var.rel]
          dim.var <- NCOL(vcov.rel)
          const.var <- matrix(NA, dim.var, 1)
          dim.in <- NCOL(x$aux$mX[, var.rel])
          indic.mat <- x$aux$mX[, var.rel]
          const.mat <- matrix(NA, NROW(x$aux$mX[, var.rel]), dim.in)
          armu.cov.sum <- matrix(NA, dim.var, 1)

          dim.ar <- NCOL(ar.var)

          if (dim.var > 1) {
            order.mat <- apply(indic.mat[, ], 2, function(x) min(which(x ==
                                                                         1)))
            indic.mat <- indic.mat[, order(order.mat)]
            vcov.rel <- vcov.rel[order(order.mat), order(order.mat)]   #ordering the variables for cumulative summing for covariances of sis

            #do the same for the autoregressive terms

            if (dim.ar > 1) #if more than one ar term
            {
              armu.cov <- armu.cov[order(order.mat), ]

            } else {
              armu.cov <- armu.cov[order(order.mat)]
            }

            for (i in 1:dim.var) {
              const.var[i] <- sum(vcov.rel[1:i, 1:i])
              if (dim.ar > 1) #if more than one ar term
              {
                armu.cov.sum[i] <- sum(armu.cov[1:i,])
              } else {
                armu.cov.sum[i] <- sum(armu.cov[1:i])
              }
            }

            const.mat <- indic.mat
            for (j in 2:(dim.in)) {
              const.mat[which(rowSums(as.matrix(const.mat[,
                                                          (j):(dim.in)])) > 0), j - 1] <- 0
            }

            #variance of ar
            ar.var.sum <- sum(ar.var)
            ar.var.part <- ((const.path^2)/((1-ar.sum)^4))*ar.var.sum

            #covariance and variance part of sis
            ind.var.mat <- const.mat %*% const.var # this one should be ok
            covar.is.part <- (1/(1-ar.sum)^2)*ind.var.mat

            #covariance part of ar and sis
            ind.ar.covar.mat <- const.mat %*% armu.cov.sum
            covar.aris.part <- 2*((const.path)/(1-ar.sum)^2)*(1/(1-ar.sum))*ind.ar.covar.mat

            lr.mean.var <- ar.var.part + covar.is.part + covar.aris.part
            lr.mean.se <- sqrt(lr.mean.var)

          } else { #if no steps retained

            ar.var.sum <- sum(ar.var)
            ar.var.part <- ((const.path^2)/((1-ar.sum)^4))*ar.var.sum

            ind.var.mat <-  vcov.rel
            covar.is.part <- (1/(1-ar.sum)^2)*ind.var.mat

            #covariance part of ar and sis

            armu.cov.sum <- sum(armu.cov)
            ind.ar.covar.mat <- armu.cov.sum
            covar.aris.part <- 2*((const.path)/(1-ar.sum)^2)*(1/(1-ar.sum))*ind.ar.covar.mat

            lr.mean.var <- ar.var.part + covar.is.part + covar.aris.part
            lr.mean.se <- sqrt(lr.mean.var)

          }

        } #if isnullxmean closed

        if(!is.null(as.list(x$call)$tis)){
          if (as.list(x$call)$tis) { #if TIS

            ind.var.mat <- NA
            ind.se.mat <- NA
            lr.path <- NA
            lr.mean.var <- NA
            lr.mean.se <- NA
          }
        }
        const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat, lr.path, lr.mean.var, lr.mean.se)
        colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se", "lr.path", "lr.var", "lr.se")

      }  else { #if there are no ar

        if(!is.null(as.list(x$call)$tis)){
          if (as.list(x$call)$tis) { #if TIS
            ind.var.mat <- NA
            ind.se.mat <- NA
          }
        }
        const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat)
        colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se")

      }  #if no ar closed
    } else {   #if no lr

      if(!is.null(as.list(x$call)$tis)){
        if (as.list(x$call)$tis) {
          ind.var.mat <- NA
          ind.se.mat <- NA
        }
      }
      const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat)
      colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se")



    }#if lr closed

    const.varse <- zoo(const.varse , order.by=x$aux$y.index)
    return(const.varse)

  } ##if (is null) closed

} #end isatvar


##==================================================
#### Impulse Indicator Consistency correction function for the sigma estimate
#Equations (9) and (15) in Johansen and Nielsen (2016)

isvarcor <- function(t.pval, sigma) {

  alpha <- t.pval
  sigmals <- sigma

  c <- abs(qnorm(alpha/2))

  psi <- pnorm(c) - pnorm(-c)
  tau <- psi - 2*c*dnorm(c)

  xi_sq <- tau/psi

  xi <- sqrt(xi_sq)
  corrxi <- 1/xi

  sigmacorr <- sigmals * corrxi

  object <-  data.frame(cbind(sigmacorr, corrxi))
  names(object) <- c("sigma.cor", "corxi")
  return(object)
}

##==================================================

###### IS Efficiency Correction for the fixed regressors
### Equation (37) in Johansen and Nielsen (2016)

isvareffcor <- function(t.pval, se, m=1) {

  alpha <- t.pval

  c <- abs(qnorm(alpha/2))

  psi <- pnorm(c) - pnorm(-c)
  tau <- psi - 2*c*dnorm(c)

  rhobeta <- 2*c*dnorm(c)/psi
  etam <- (((1-rhobeta^m)/((1-rhobeta)*psi))^2 + 2 * (1-rhobeta^m)/((1-rhobeta)*psi)*rhobeta^m)*tau+rhobeta^(2*m)

  se_cor <- se*sqrt(etam)

  output <-  data.frame(cbind(se_cor, sqrt(etam) ))
  names(output) <- c("se.cor", "eta.m")
  return(output)
}

#OLD VERSIONS OF biascorr...etc.:
###==================================================
### needs orthogonal variables, can be applied directly to the coefficient path though in an isat only model
#### purposefully not entirely usefriendly - needs many arguments, but it should be an expert user function
#biascorr <- function(b, b.se, p.alpha, T){
#
#  c_alpha <- abs(qt((1-(1-p.alpha))/2, T))
#  bt <- b/b.se
#
#  dr <- (dnorm(c_alpha-bt)-dnorm(-c_alpha-bt))/ (1-pnorm(c_alpha-bt) + pnorm(-c_alpha-bt))
#  dtbar <- bt - dr
#  drbar <- (dnorm(c_alpha-dtbar)-dnorm(-c_alpha-dtbar))/ (1-pnorm(c_alpha-dtbar) + pnorm(-c_alpha-dtbar))
#
#  #only correct if significant
#
#  b_1step <- b
#  b_1step[abs(bt)>c_alpha] <- b*(1-(dr/bt))[abs(bt)>c_alpha]
#
#  b_2step <- b
#  b_2step[abs(bt)>c_alpha] <- b*(1-(drbar/bt))[abs(bt)>c_alpha]
#
#  b_corr <- cbind(b, b_1step, b_2step)
#  colnames(b_corr) <- c("beta", "beta.1step", "beta.2step")
#
#  return(b_corr)
#}
#
###==================================================
######## isat.test: forecast bias test as based on the working paper with a plot
#### loads an isat object, and then conducts the bias analysis on it
#
#isattest <- function(x, hnull=0, lr=FALSE, ci.pval = 0.99, plot=TRUE, plot.turn = FALSE, biascorr=FALSE){
#
#  trend.incl <- FALSE
#  if(!is.null(as.list(x$call)$tis)){
#    if (as.list(x$call)$tis) {
#      stop("isat.test currently not implemented for trend-indicator saturation")
#      trend.incl <- TRUE
#    }
#  }
#
#  arcall <- as.list(x$call)$ar
#  x.var <- isatvar(x,lr=lr)
#
#  if (biascorr==TRUE){
#
#    if (!is.null(as.list(x$call)$mxreg) | !is.null(arcall) | trend.incl){
#
#      biascorr <- FALSE
#      print("Bias Correction not applicable in isat regression with additional non-step covariates. Has been set to FALSE.")
#    }
#  }
#
#
#  T <- dim(x$aux$mX)[1]
#  N <- dim(x$aux$mX)[2]
#
#  crit <- abs(qt((1-ci.pval)/2, T-N))
#  bias.low <- matrix(0, T, 1)
#  bias.high <- matrix(0, T, 1)
#
#  ci.low <- matrix(0, T, 1)
#  ci.high <- matrix(0, T, 1)
#  x.mean <- matrix(0, T, 1)
#
#  if (lr == TRUE & !is.null(arcall))
#  {
#    x.is.lr <- x.var$lr.path
#    x.is.const <- x.var$const.path
#
#  } else {
#
#    x.is.lr <- NA
#
#    if (biascorr){
#
#      #is2$aux$t.pval
#
#      xbias <-biascorr(b=x.var$const.path, b.se=x.var$const.se, p.alpha = x$aux$t.pval, T=length(x.var$const.path))
#
#      #xbias <-biascorr(b=x.var$const.path, b.se=x.var$const.se, p.alpha = as.list(x$call)$t.pval, T=length(x.var$const.path))
#      x.is.const <- xbias$beta.2step
#
#    } else {
#      x.is.const <- x.var$const.path
#
#    }
#
#
#  }
#
#
#  if (lr == TRUE & !is.null(arcall))
#  {
#
#    ci.low <- x.var$lr.path-crit*x.var$lr.se
#    ci.high <- x.var$lr.path+crit*x.var$lr.se
#
#    bias.low[which((ci.low) > hnull)] <- 1
#    bias.high[which((ci.high) < hnull)] <- 1
#
#    bias.low <- bias.low*(x.var$lr.path-hnull)
#    bias.high <- bias.high*(x.var$lr.path-hnull)
#
#    x.mean <- x.var$lr.path
#
#  } else {
#
#    ci.low <- x.is.const-crit*x.var$const.se
#    ci.high <- x.is.const+crit*x.var$const.se
#
#    bias.low[which((ci.low) > hnull)] <- 1
#    bias.high[which((ci.high) < hnull)] <- 1
#
#
#    bias.low <- bias.low*(x.is.const-hnull)
#    bias.high <- bias.high*(x.is.const-hnull)
#
#    x.mean <- x.is.const
#
#  }
#
#  ###determining the turning points
#  time <- x$aux$y.index
#  bias.sum.ar <- bias.low+bias.high
#  lr.path.d <- diff(bias.sum.ar)
#
#
#  if(all(lr.path.d==0)){
#    plot.turn <- FALSE
#    turn.ar <- NULL
#  } else {
#    turn.ar <- time[which(lr.path.d != 0)]+1
#  }
#
#  turn.ar.y <- bias.sum.ar[which(lr.path.d != 0)]
#
#  turn.x.lab <- turn.ar
#  turn.x <- turn.ar
#
#  fitted <- x$mean.fit
#  actual <- zoo(x$aux$y, order.by=x$aux$y.index)
#
#  ylabel_a <- "Series"
#  ylabel_b <- "Bias"
#
#
#  par(mfrow=c(2,1), mar = c(2, 4,1,3))
#  Ylim_main <- c(min(actual, na.rm=TRUE)*1.2,max(actual, na.rm=TRUE)*1.2)
#  Ylim_bias <- c(min(bias.high, na.rm=TRUE)*1.2,max(bias.low, na.rm=TRUE)*1.2)
#
#  if (plot){
#
#    plot(time, x.mean, ylim=Ylim_main, col="blue", title(main=NULL, xlab=NULL), xlab=NA, ylab=ylabel_a, sub=NA, type="l")
#    lines(ci.low, col="blue", lty=2)
#    lines(ci.high, col="blue", lty=2)
#    lines(actual)
#    abline(a =hnull, b=0, col="black", lty=3, lwd=2)
#
#    plot(time, bias.low, type="h", col="red", ylim=Ylim_bias, title(main=NULL, xlab=NULL), xlab=NA, ylab=ylabel_b, sub=NA)
#    lines(bias.high, type="h", col="blue")
#
#    if ( plot.turn ){
#      text(turn.x.lab, y=turn.ar.y, x=turn.x, pos=4, offset=-0.5, cex=0.8)
#    }
#
#  }
#
#
#  if (lr==TRUE & !is.null(arcall))
#  {
#    mean.var <- cbind(ci.low, ci.high, bias.low,  bias.high)
#    colnames(mean.var) <- c("ci.low", "ci.high", "bias.high",  "bias.low")
#
#
#  } else {
#    mean.var <- cbind( ci.low, ci.high, bias.low,  bias.high)
#    colnames(mean.var) <- c("ci.low", "ci.high", "bias.high",  "bias.low")
#
#  }
#
#
#  return(mean.var)
#
#} #isat.test function closedakes an isat gets object, returns:
## - the coefficient path (both relative to the full sample coefficient and the constant path itself)
## - the variance and standard errors of the coefficient path
## - if lr is specified and the isat object has AR terms, then also returns the LR coefficient path and its variance and standard errors
##input: "x" (an isat results object)
#
#
#isatvar <- function(x, lr=FALSE)
#{
#
#  if (!is.null(x$mean.fit)){
#    if (!is.null(x$ISnames)){
#      var.rel <- c( which(substr(x$aux$mXnames,1,6) %in% "mconst"), which(x$aux$mXnames %in% x$ISnames))  #vector of where the constant and is terms are
#
#      #coefficient path function
#      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% x$ISnames])
#      is.coef.ests <- coef.isat(x)[x$ISnames]
#      coef.path <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
#      colnames(coef.path) <- "coef.path"
#
#      const <- coef.isat(x)["mconst"]
#      const.path <- coef.path + const
#      colnames(const.path) <- "const.path"
#
#    } else {
#      var.rel <- which(substr(x$aux$mXnames,1,6) %in% "mconst")
#
#      #coefficient path function
#      const <- coef.isat(x)["mconst"]
#
#      coef.path <- zoo(0, order.by = x$aux$y.index)
#
#      const.path <- coef.path + const
#
#    }
#
#    vcov.rel <- x$vcov.mean[var.rel,var.rel]
#    dim.var <- NCOL(vcov.rel)
#    const.var <- matrix(NA, dim.var, 1 )
#
#    #construct a matrix to multiply by the variances
#    dim.in <- NCOL(x$aux$mX[,var.rel])
#    indic.mat <- x$aux$mX[,var.rel]
#    const.mat <- matrix(NA, NROW(x$aux$mX[,var.rel]), dim.in)
#
#
#    if (dim.var > 1) #if there are indicators retained
#    {
#
#      #order the breaks in correct order so the covariance matrix and relevance order is correct
#      order.mat <- apply(indic.mat[,], 2, function(x) min(which(x==1)))  #find where each indicator is first one, for sorting
#
#      indic.mat <- indic.mat[,order(order.mat)]  #order the indicators
#      vcov.rel <- vcov.rel[order(order.mat),order(order.mat)] #order the covariance matrix
#
#      for (i in 1:dim.var)
#      {
#        const.var[i] <- sum(vcov.rel[1:i,1:i])   #sum over the expanding variance covariance matrix to get the variance of the sums of coefficients
#      }
#
#      const.mat <- indic.mat
#
#      for (j in 2:(dim.in))
#      {
#
#        const.mat[which(rowSums(as.matrix(const.mat[,(j):(dim.in)]))>0),j-1] <- 0     #puts zeros in the appropriate places to make sure the correct s.e. is applied for each point in time
#
#      }
#
#      ind.var.mat <- const.mat %*% const.var   #the variance as it applies to each subsection
#
#    } else { #just the constant remains
#
#      const.var <- vcov.rel
#      const.mat <- indic.mat
#      ind.var.mat <- const.mat * const.var
#
#    }
#
#    ind.se.mat <- sqrt(ind.var.mat)   #the standard errors of the coefficient as it changes with the SIS breaks
#
#    ####the coefficient path of the LR mean
#
#    if(lr){
#      if (!is.null(as.list(x$call)$ar)){
#
#        if (!is.null(x$mean.fit)){
#
#          vcov.rel.tot <- x$vcov.mean
#          coef.rel <- coef.isat(x)
#
#          ###coefficient path of LR mean
#
#          arcall <- as.list(x$call)$ar
#          arnames <- paste("ar",eval(as.expression(arcall)), sep="")
#
#          ar.coefs <- coef.isat(x)[arnames]
#          ar.sum <- sum(ar.coefs)
#          lr.path <- const.path/(1-ar.sum)
#
#          ###variance of the coefficient path
#
#          ar.var <- vcov.rel.tot[arnames,arnames] #variance of ar() terms
#          armu.cov <- vcov.rel.tot[c("mconst", x$ISnames), arnames] #covariance of ar(1) and const + sis
#
#
#          if (!is.null(x$ISnames)) {
#            var.rel <- c(which(substr(x$aux$mXnames, 1, 6) %in%
#                                 "mconst"), which(x$aux$mXnames %in% x$ISnames))
#          } else {
#            var.rel <- which(substr(x$aux$mXnames, 1, 6) %in%
#                               "mconst")
#          }
#
#          vcov.rel <- x$vcov.mean[var.rel, var.rel]
#          dim.var <- NCOL(vcov.rel)
#          const.var <- matrix(NA, dim.var, 1)
#          dim.in <- NCOL(x$aux$mX[, var.rel])
#          indic.mat <- x$aux$mX[, var.rel]
#          const.mat <- matrix(NA, NROW(x$aux$mX[, var.rel]), dim.in)
#          armu.cov.sum <- matrix(NA, dim.var, 1)
#
#          dim.ar <- NCOL(ar.var)
#
#          if (dim.var > 1) {
#            order.mat <- apply(indic.mat[, ], 2, function(x) min(which(x ==
#                                                                         1)))
#            indic.mat <- indic.mat[, order(order.mat)]
#            vcov.rel <- vcov.rel[order(order.mat), order(order.mat)]   #ordering the variables for cumulative summing for covariances of sis
#
#            #do the same for the autoregressive terms
#
#            if (dim.ar > 1) #if more than one ar term
#            {
#              armu.cov <- armu.cov[order(order.mat), ]
#
#            } else {
#              armu.cov <- armu.cov[order(order.mat)]
#            }
#
#            for (i in 1:dim.var) {
#              const.var[i] <- sum(vcov.rel[1:i, 1:i])
#              if (dim.ar > 1) #if more than one ar term
#              {
#                armu.cov.sum[i] <- sum(armu.cov[1:i,])
#              } else {
#                armu.cov.sum[i] <- sum(armu.cov[1:i])
#              }
#            }
#
#            const.mat <- indic.mat
#            for (j in 2:(dim.in)) {
#              const.mat[which(rowSums(as.matrix(const.mat[,
#                                                          (j):(dim.in)])) > 0), j - 1] <- 0
#            }
#
#            #variance of ar
#            ar.var.sum <- sum(ar.var)
#            ar.var.part <- ((const.path^2)/((1-ar.sum)^4))*ar.var.sum
#
#            #covariance and variance part of sis
#            ind.var.mat <- const.mat %*% const.var # this one should be ok
#            covar.is.part <- (1/(1-ar.sum)^2)*ind.var.mat
#
#            #covariance part of ar and sis
#            ind.ar.covar.mat <- const.mat %*% armu.cov.sum
#            covar.aris.part <- 2*((const.path)/(1-ar.sum)^2)*(1/(1-ar.sum))*ind.ar.covar.mat
#
#            lr.mean.var <- ar.var.part + covar.is.part + covar.aris.part
#            lr.mean.se <- sqrt(lr.mean.var)
#
#          } else { #if no steps retained
#
#            ar.var.sum <- sum(ar.var)
#            ar.var.part <- ((const.path^2)/((1-ar.sum)^4))*ar.var.sum
#
#            ind.var.mat <-  vcov.rel
#            covar.is.part <- (1/(1-ar.sum)^2)*ind.var.mat
#
#            #covariance part of ar and sis
#
#            armu.cov.sum <- sum(armu.cov)
#            ind.ar.covar.mat <- armu.cov.sum
#            covar.aris.part <- 2*((const.path)/(1-ar.sum)^2)*(1/(1-ar.sum))*ind.ar.covar.mat
#
#            lr.mean.var <- ar.var.part + covar.is.part + covar.aris.part
#            lr.mean.se <- sqrt(lr.mean.var)
#
#          }
#
#        } #if isnullxmean closed
#
#        if(!is.null(as.list(x$call)$tis)){
#          if (as.list(x$call)$tis) { #if TIS
#
#            ind.var.mat <- NA
#            ind.se.mat <- NA
#            lr.path <- NA
#            lr.mean.var <- NA
#            lr.mean.se <- NA
#          }
#        }
#        const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat, lr.path, lr.mean.var, lr.mean.se)
#        colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se", "lr.path", "lr.var", "lr.se")
#
#      }  else { #if there are no ar
#
#        if(!is.null(as.list(x$call)$tis)){
#          if (as.list(x$call)$tis) { #if TIS
#            ind.var.mat <- NA
#            ind.se.mat <- NA
#          }
#        }
#        const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat)
#        colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se")
#
#      }  #if no ar closed
#    } else {   #if no lr
#
#      if(!is.null(as.list(x$call)$tis)){
#        if (as.list(x$call)$tis) {
#          ind.var.mat <- NA
#          ind.se.mat <- NA
#        }
#      }
#      const.varse <- cbind(coef.path, const.path, ind.var.mat, ind.se.mat)
#      colnames(const.varse) <- c("coef.path", "const.path", "const.var", "const.se")
#
#
#
#    }#if lr closed
#
#    const.varse <- zoo(const.varse , order.by=x$aux$y.index)
#    return(const.varse)
#
#  } ##if (is null) closed
#
#} #end isatvar


##==================================================
## make matrix of impulse indicators:
iim <- function(x, which.ones=NULL)
{
  if(NROW(x)==1){
    n <- x
    mIIS <- matrix(0,n,n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", 1:n, sep="")
    if(!is.null(which.ones)){ mIIS <- mIIS[,which.ones] }
    mIIS <- as.zoo(mIIS)
  }else{
    n <- NROW(x)
    mIIS <- matrix(0,n,n)
    diag(mIIS) <- 1
    x <- as.zoo(x)
    x.index <- index(x)
    xIsRegular <- is.regular(x, strict=TRUE)
    if(xIsRegular && frequency(x)>1 ){
      xIndexObs <- floor(as.numeric(x.index))
      xCycle <- as.numeric(cycle(x))
      xIndexAsChar <- paste(xIndexObs, "(", xCycle, ")", sep="")
      xFrequency <- frequency(x)
    }else{
      xIndexAsChar <- as.character(x.index)
    }
    colnames(mIIS) <- paste("iis",
      xIndexAsChar, sep="")
    mIIS <- zoo(mIIS, order.by=x.index)
    if(xIsRegular){ mIIS <- as.zooreg(mIIS) }
    if(!is.null(which.ones)){
      where.indicators <- which(index(mIIS) %in% which.ones)
      if(length(where.indicators > 0)){
        mIIS <- cbind(mIIS[,where.indicators])
      }else{
        stop("'which.ones' not in index")
      }
    }
  } #end if(NROW(x)==1)else..
  return(mIIS)
} #close iim

##==================================================
##make matrix of step indicators:
sim <- function(x, which.ones=NULL)
{
  if(NROW(x)==1){
    x.is.scalar <- TRUE
    n <- x
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which.ones
    }
  }else{
    x.is.scalar <- FALSE
    n <- NROW(x)
    x <- as.zoo(x)
    x.index <- index(x)
    xIsRegular <- is.regular(x, strict=TRUE)
    if(xIsRegular && frequency(x)>1 ){
      xIndexObs <- floor(as.numeric(x.index))
      xCycle <- as.numeric(cycle(x))
      xIndexAsChar <- paste(xIndexObs, "(", xCycle, ")", sep="")
      xFrequency <- frequency(x)
    }else{
      xIndexAsChar <- as.character(x.index)
    }
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(x.index %in% which.ones)
      if(length(where.indicators)==0) stop("'which.ones' not in index")
    } #end if(is.null(which.ones))
  } #end if(NROW(x)==1)else(..)

  n.where.indicators <- length(where.indicators)
  loop.indx <- 1:n.where.indicators
  mSIS <-matrix(0,n,n.where.indicators)
  tmp <- function(i){ mSIS[ c(where.indicators[i]:n) ,i] <<- 1 }
  tmp <- sapply(loop.indx,tmp)
  if(x.is.scalar){
    colnames(mSIS) <- paste("sis", where.indicators, sep="")
    mSIS <- as.zoo(mSIS)
  }else{
    colnames(mSIS) <- paste("sis",
      xIndexAsChar[where.indicators], sep="")
    if(xIsRegular && frequency(x)>1 ){
      mSIS <- zooreg(mSIS, frequency=xFrequency,
        start=c(xIndexObs[1],xCycle[1]))
    }else{
      mSIS <- zoo(mSIS, order.by=x.index)
    }
  }
  return(mSIS)
} #close sim

##==================================================
## make matrix of trend indicators:
tim <- function(x, which.ones=NULL, log.trend=FALSE)
{
  if(NROW(x)==1){
    x.is.scalar <- TRUE
    n <- x
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(1:n %in% which.ones)
    }
  }else{
    x.is.scalar <- FALSE
    n <- NROW(x)
    x <- as.zoo(x)
    x.index <- index(x)
    xIsRegular <- is.regular(x, strict=TRUE)
    if(xIsRegular && frequency(x)>1 ){
      xIndexObs <- floor(as.numeric(x.index))
      xCycle <- as.numeric(cycle(x))
      xIndexAsChar <- paste(xIndexObs, "(", xCycle, ")", sep="")
      xFrequency <- frequency(x)
    }else{
      xIndexAsChar <- as.character(x.index)
    }
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(x.index %in% which.ones)
      if(length(where.indicators)==0) stop("'which.ones' not in index")
    }
  }
  n.where.indicators <- length(where.indicators)
  mTIS <-matrix(0,n,n.where.indicators)
  v1n <- seq(1,n)
  loop.indx <- 1:n.where.indicators
  tmp <- function(i){
    t.trend <- v1n[1:c(n-where.indicators[i]+1)]
    if(log.trend) t.trend <- log(t.trend)
    mTIS[c(where.indicators[i]:n),i] <<- t.trend
  }
  tmp <- sapply(loop.indx,tmp)
  if(x.is.scalar){
    colnames(mTIS) <- paste("tis", where.indicators, sep="")
    mTIS <- as.zoo(mTIS)
  }else{
    colnames(mTIS) <- paste("tis",
      xIndexAsChar[where.indicators], sep="")
    mTIS <- zoo(mTIS, order.by=x.index)
    if(xIsRegular){ mTIS <- as.zooreg(mTIS) }
  }
  return(mTIS)
} #close tim


################################
################# New Functions July 2019
##############################

##############
### isatdates
#############

## Extract breakdates from isat object (based on first break index):

isatdates <- function(x){
  
  mxbreak_iis <- "iis"
  mxbreak_sis <- "sis"
  mxbreak_tis <- "tis"
  
  iis.names <- c(x$ISnames[grep(mxbreak_iis, x$ISnames)])
  sis.names <- c(x$ISnames[grep(mxbreak_sis, x$ISnames)])
  tis.names <- c(x$ISnames[grep(mxbreak_tis, x$ISnames)])
  
  ##### iis
  if(length(iis.names) != 0){
    iis.breaks <- data.frame(matrix(NA, nrow=NROW(iis.names), ncol=1))
    names(iis.breaks) <- c("breaks")  
    is.m <- as.matrix(x$aux$mX[,iis.names])
    is.index <- which(x$aux$mXnames %in% iis.names)
    colnames(is.m) <- iis.names
    iis.date <- apply(is.m,2, function(x) (which(x>0))[1])  
    iis.date.index <- iis.date
    iis.date <- x$aux$y.index[iis.date]
    
    iis.breaks$breaks <- iis.names
    iis.breaks$date <- iis.date
    iis.breaks$index <- iis.date.index
    iis.breaks$coef <- x$mean.results$coef[is.index]
    iis.breaks$coef.se <- x$mean.results$std.error[is.index]
    iis.breaks$coef.t <- x$mean.results$`t-stat`[is.index]
    iis.breaks$coef.p <- x$mean.results$`p-value`[is.index]
    
  } else {
    iis.breaks <- NULL 
  }
  
  ##### sis
  if(length(sis.names) != 0){
    sis.breaks <- data.frame(matrix(NA, nrow=NROW(sis.names), ncol=1))
    names(sis.breaks) <- c("breaks")  
    sis.m <- as.matrix(x$aux$mX[,sis.names])
    sis.index <- which(x$aux$mXnames %in% sis.names)
    colnames(sis.m) <- sis.names
    sis.date <- apply(sis.m,2, function(x) (which(x>0))[1])  
    sis.date.index <- sis.date
    sis.date <- x$aux$y.index[sis.date]
    
    sis.breaks$breaks <- sis.names
    sis.breaks$date <- sis.date
    sis.breaks$index <- sis.date.index
    sis.breaks$coef <- x$mean.results$coef[sis.index]
    sis.breaks$coef.se <- x$mean.results$std.error[sis.index]
    sis.breaks$coef.t <- x$mean.results$`t-stat`[sis.index]
    sis.breaks$coef.p <- x$mean.results$`p-value`[sis.index]
    
  } else {
    sis.breaks <- NULL 
  }
  
  ##### tis
  if(length(tis.names) != 0){
    tis.breaks <- data.frame(matrix(NA, nrow=NROW(tis.names), ncol=1))
    names(tis.breaks) <- c("breaks")  
    tis.m <- as.matrix(x$aux$mX[,tis.names])
    tis.index <- which(x$aux$mXnames %in% tis.names)
    colnames(tis.m) <- tis.names
    tis.date <- apply(tis.m,2, function(x) (which(x>0))[1])  
    tis.date.index <- tis.date
    tis.date <- x$aux$y.index[tis.date]
    
    tis.breaks$breaks <- tis.names
    tis.breaks$date <- tis.date
    tis.breaks$index <- tis.date.index
    tis.breaks$coef <- x$mean.results$coef[tis.index]
    tis.breaks$coef.se <- x$mean.results$std.error[tis.index]
    tis.breaks$coef.t <- x$mean.results$`t-stat`[tis.index]
    tis.breaks$coef.p <- x$mean.results$`p-value`[tis.index]
    
  } else {
    tis.breaks <- NULL 
  }
  
  out <- list(iis.breaks, sis.breaks, tis.breaks)
  names(out) <- c("iis", "sis", "tis")
  return(out)
  
} #end function


####################
### isatvarcorrect
#####################
#### Function to correct variance estimates when using IIS

isatvarcorrect <- function(x,   mcor = 1){
  
  if (class(x)=="isat"){
    if (x$call$iis==TRUE) 
    {
      x$vcov.mean <- x$vcov.mean * as.numeric(isvarcor(x$aux$t.pval, 1)[2]^2)
      x$vcov.mean[x$keep, x$keep] <- x$vcov.mean[x$keep, x$keep] * as.numeric(isvareffcor(x$aux$t.pval, 1, mcor)[2]^2)
      
      x$mean.results$std.error <- sqrt(diag(x$vcov.mean))
      x$mean.results$`t-stat` <- x$mean.results$coef/x$mean.results$std.error
      x$mean.results$`p-value` <-  pt(abs(x$mean.results$`t-stat`), x$df, lower.tail=FALSE)*2   
      
      x$sigma2 <- x$sigma2*as.numeric(isvarcor(x$aux$t.pval, 1)[2]^2)
      x$logl <- -x$n*log(2*x$sigma2*pi)/2 - x$rss/(2*x$sigma2)
      
      return(x)
    } else {
      stop("iis not TRUE")
    }
  } else {
    stop("x must be an isat object")
  }  
  
} #function closed


###################
####### vargaugeiis
###################

##### Function to Compute Variance of the gauge (for use in outliertest)

vargaugeiis <- function(t.pval, T, infty=FALSE, m=1){
  
  alpha <- t.pval
  c <- abs(qnorm(alpha/2))
  fc <- dnorm(c)
  psi <- pnorm(c) - pnorm(-c) 
  tau <- psi - 2*c*dnorm(c)
  chi <- 3*psi - 2*c*(c^2+3)*fc
  rho_sig <- (c^2 - tau/psi)*c*fc/tau
  k <- 3
  
  ###m iterations
  m_min1 <- m-1
  eta_sigma_m_min1 <- ( ((1-rho_sig^m_min1)/((1-rho_sig)*tau))^2 + 2* ((1-rho_sig^m_min1)/((1-rho_sig)*tau)  )*rho_sig^m_min1  ) *(chi-tau^2/psi)/(k-1)+rho_sig^(2*m_min1)
  eta_gamma_m <- (c*fc)^2*eta_sigma_m_min1*(k-1) + 2*c*fc*rho_sig^(m-1)*(tau-psi)
  
  ###infinite iterations
  eta_infty <- (chi-tau^2/psi)*(c*fc)^2/((1-rho_sig)*tau)^2
  
  if (infty==TRUE){
    viis <- psi*(1-psi) + eta_infty  
  } else {
    viis <- psi*(1-psi) + eta_gamma_m
  }
  
  sdiis <- sqrt(viis)
  viis_T <- viis/T
  sdiis_T <- sqrt(viis/T)
  
  out <- data.frame(cbind(viis_T, sdiis_T, viis, sdiis))
  names(out) <- c("var_iisgauge", "sd_iisgauge", "asy_var_iisgauge", "asy_sd_iisgauge")
  return(out)
}  

################
#### outliertest
################

#### Outlier proportion and count tests from Jiao and Pretis (2019)


outliertest <- function(x=NULL, noutl=NULL, t.pval=NULL, T=NULL,  m=1, infty=FALSE, alternative="two.sided"){
  
  # noutl=x$obs.gauge$obs.num[i]
  # t.pval=x$obs.gauge$t.pval[i]
  # T=x$n
  # 
  
  if (!is.null(x)){
    
    if (class(x)=="isat"){
      
      if (any(x$call$sis, x$call$tis, x$call$uis)==TRUE){
        stop("Test only valid for iis")
      } else {
        if (  x$call$iis == TRUE){
          ISnames <- c(x$ISnames[grep("iis", x$ISnames)])
          noutl= length(ISnames) 
          t.pval <- x$aux$t.pval
          T <- x$n
        } else {
          stop("iis must be TRUE")
        }
      }
      
      
    } else {
      stop("x must be an isat object")
    }
  } #x null closed
  
  gauge.null <- t.pval
  gauge <- noutl/T
  
  #### Standard Normal Test
  gauge.sd <- vargaugeiis(t.pval=gauge.null, T=T, infty=infty, m=m)$sd_iisgauge 
  gaugetestval <- (gauge - gauge.null)/gauge.sd
  
  if (alternative=="two.sided"){
    pval <- 2*(pnorm(-abs(gaugetestval)))
  }
  if (alternative=="less"){
    pval <- pnorm(gaugetestval)  
  }
  if (alternative=="greater"){
    pval <- 1-  pnorm(gaugetestval) 
  }
  
  ### Poisson Test
  gauge_n <- noutl
  gauge_null_n <- gauge.null*T
  poistestval <- poisson.test(gauge_n, r=gauge_null_n, alternative=alternative)
  
  rval_norm <- list(statistic = gaugetestval, p.value = pval, estimate=gauge, null.value = gauge.null, alternative = alternative, method="Jiao-Pretis Outlier Proportion Test", data.name="Proportion of detected outliers")
  attr(rval_norm, "class") <- "htest"
  rval_pois <- list(statistic = gauge_n, p.value = poistestval$p.value, estimate=gauge_n, null.value = gauge_null_n, alternative = alternative, method="Jiao-Pretis Outlier Count Test", data.name="Number of detected outliers")
  attr(rval_pois, "class") <- "htest"
  out <- list(rval_norm, rval_pois)
  names(out) <- c("proportion", "count")
  
  return(out)
}

#######################
###mvrnormsim: reproduced from MASS
#######################

mvrnormsim <- function(n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE){
  
  p <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) {
    nm <- dn[[1L]]  
  }
  dimnames(X) <- list(nm, NULL)
  if (n == 1) {
    drop(X) } else {
      t(X)   
    } 
  #return(X)
}


####################################
############# isatloop
####################################

##### Looping over isat iterations for different p-values of selection to be used in outlierscaletest
isatloop <- function(num=c(seq(from=20, to=1, by=-1)), t.pval.spec = FALSE, print=FALSE, y, ar=NULL, iis=TRUE,  sis=FALSE,  ...){
  
  
  initialm <- arx(y=y, ar=ar) 
  n <- initialm$n
  num <- num[rev(order(num))]
  
  if (t.pval.spec == FALSE){
    p.num <- num/n # scaling significance levels
  } else { #if p-values are pre-specified
    p.num <- num
  }
  
  K <- length(p.num)
  obs.gauge <- data.frame(matrix(NA, nrow=K, ncol=4)) # record significance p value, sample gauge, expected number and sample number of outliers
  names(obs.gauge) <- c("t.pval", "obs.prop", "p.num", "obs.num")
  for (k in 1:K) 
  {
    pval <- p.num[k]
    if (print==TRUE)
    {
      print(paste("k:", k, "/", K, ", p:", pval, sep=""))
    }
    x <- isat(y=y, iis=iis, sis=sis, t.pval=pval, ar=ar, ...) 
    
    #mxreg=mxreg, mc=mc, 
    if (!is.null(x$ISnames)) {
      ISnames <- c(x$ISnames[grep("iis", x$ISnames)])
      noutl= length(ISnames) 
      is.gauge <- noutl/n
      is.num <- noutl
    }  else   {
      is.gauge <- 0
      is.num <- 0
    }
    if (print==TRUE)
    {
      print(paste("gauge k:", is.gauge, ", number k:", "is.num", sep=""))     
    }
    obs.gauge$t.pval[k] <- pval 
    obs.gauge$obs.prop[k] <- is.gauge
    obs.gauge$p.num[k] <- pval*n
    obs.gauge$obs.num[k] <- is.num    
  }
  
  out <- list(n, obs.gauge)
  names(out) <- c("n", "obs.gauge")
  return(out)
}

#######################
##### outlierscaletest
########################

#### Scaling outlier tests from Jiao and Pretis (2019)

outlierscaletest <- function(x, nsim = 10000){
  ###################
  
  obs.gauge <- x$obs.gauge
  n <- x$n ##need to extract n
  p.num <- obs.gauge$t.pval
  c <- qnorm(1 - p.num/2) # corresponding cut-offs # assume standardised error follows standard normal
  fc <- dnorm(c) 
  psi <- 1 - p.num
  tau_2_c <- psi - 2*c*fc
  tau_4 <- 3
  K <- length(p.num)
  
  #######################
  ##### Scale Sum Test
  ######################
  
  stand.obs.gauge <- n^(1/2)*(obs.gauge$obs.prop - obs.gauge$t.pval) # standardise gauge
  scalesum.stat <- sum(stand.obs.gauge) # scaling sum statistic
  
  if (K == 1) # covariance structure between test statistics of different significance levels
  {
    covar <- 0
  }  else  {
    covar <- 0 
    for (k in 1:(K - 1)) 
    {
      for (l in (k + 1):K)
      {
        covar <- covar + obs.gauge$t.pval[l] - obs.gauge$t.pval[k]*obs.gauge$t.pval[l]
      }
    }
  }
  var.scalesum.stat <- sum(obs.gauge$t.pval*(1 - obs.gauge$t.pval)) + 2*covar + sum(c*fc)^(2)*(tau_4 - 1) + 2*sum(c*fc)*sum(tau_2_c + obs.gauge$t.pval - 1) # variance for scaling sum test statistic
  sd.scalesum.stat <- sqrt(var.scalesum.stat)
  
  #### Sum Test Output
  stand.scalesum.stat <- scalesum.stat/sd.scalesum.stat
  scalesum.pval <- 2*(pnorm(-abs(stand.scalesum.stat)))
  
  #######################
  ##### Scale Sup Test
  ######################
  
  scalesup.stat <- max(abs(stand.obs.gauge)) # scaling sup statistic
  
  N <- nsim # sample size used to simulate the limiting distribution (could also be added as the argument of function)
  mu <- matrix(0, K, 1) # mean of GP
  Sigma <- matrix(NA, K, K) # covariance of GP
  for (s in 1:K)
  {
    for (t in 1:K)
    {
      if (s <= t)
      {
        Sigma[s, t] <- obs.gauge$t.pval[t]*(1 - obs.gauge$t.pval[s]) + c[s]*fc[s]*(tau_2_c[t] + obs.gauge$t.pval[t] - 1) + c[t]*fc[t]*(tau_2_c[s] + obs.gauge$t.pval[s] - 1) + c[s]*c[t]*fc[s]*fc[t]*(tau_4 - 1)
      }
      else
      {
        Sigma[s, t] <- Sigma[t, s]
      }
    }
  }
  GPsample <- mvrnormsim(N, mu, Sigma) # generate multivariate normal and dim of object is N by K
  Limitsample <- apply(abs(GPsample), 1, max) # find largest value over rows of absolute of GPsample
  
  pvalsample <- Limitsample # use empirical distribution to draw p value
  pvalsample[Limitsample <= scalesup.stat] <- 0
  pvalsample[Limitsample > scalesup.stat] <- 1
  scalesup.pval <- mean(pvalsample)
  

  ##################################
  ##### Output of Sum and Sup Tests
  
  rval_sum <- list(statistic = stand.scalesum.stat, p.value = scalesum.pval, estimate=NULL, null.value = NULL, alternative = NULL, method="Jiao-Pretis Outlier Scaling Sum Test", data.name="Scaling proportion of detected outliers (Sum)")
  attr(rval_sum, "class") <- "htest"
  rval_sup <- list(statistic = scalesup.stat, p.value = scalesup.pval, estimate=NULL, null.value = NULL, alternative = NULL, method="Jiao-Pretis Outlier Scaling Sup Test", data.name="Scaling proportion of detected outliers (Sup)")
  attr(rval_sup, "class") <- "htest"
  
  out <- list(rval_sum, rval_sup)
  names(out) <- c("sum", "sup")
  
  return(out)
  
} ###function closed






