##==================================================
##do block-based gets with full flexibility (for advanced users)
blocksFun <- function(y, x, blocks=NULL, no.of.blocks=NULL, 
  max.block.size=30, ratio.threshold=0.8, user.estimator=list(name="ols"),
  t.pval=0.001, wald.pval=t.pval, do.pet=FALSE, ar.LjungB=NULL,
  arch.LjungB=NULL, normality.JarqueB=NULL, user.diagnostics=NULL,
  gof.function=list(name = "infocrit", method = "sc"),
  gof.method=c("min","max"), keep=NULL, include.gum=FALSE,
  include.1cut=FALSE, include.empty=FALSE, max.paths=NULL,
  parallel.options=NULL, turbo=FALSE, force.invertibility=TRUE,
  tol=1e-07, LAPACK=FALSE, max.regs=NULL, print.searchinfo=TRUE,
  plot=NULL, alarm=FALSE)
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
  result$time.started <- date()
  result$time.finished <- NA
  result$call <- sys.call()
  result$messages <- NULL

  ##regressand:
  y.n <- NROW(y) #needed to determine no.of.blocks

#OLD (isat):
#  ##include.gum argument:
#  include.gum <- TRUE

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
  ## 2 x and blocks arguments
  ##-------------------------------
  
  ##if x is a matrix:
  if( is.matrix(x) ){

    ##matrix name:
    xMatrixName <- deparse(substitute(x))
    
    ##handle colnames:
    xColNames <- colnames(x)
    if( is.null(xColNames) ){
      xColNames <- paste0("xreg", 1:NCOL(x))
    }
    if( any(xColNames=="") ){
      missing.colnames <- which(xColNames == "")
      for(i in 1:length(missing.colnames)){
        #fixed by Jonas: 
        xColNames[ missing.colnames[i] ] <-
          paste0("xreg", missing.colnames[i])
      }
    }
    xColNames <- make.unique(xColNames)
    colnames(x) <- xColNames

    ##convert to list:
    x <- list(x=x)
    names(x) <- xMatrixName
  } #end if( is.matrix(x) )

  ##if x is a list of matrices:
  if( is.list(x) && length(x)>1 ){

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

    ##to do: check that colnames() are unique across matrices
    
  } #end if( is.list(x) )

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


  ##-------------------------------
  ## 3 loop on x matrices
  ##-------------------------------

  xGETSofUnionModels <- list()
  for(i in 1:length(x)){

    ##xGETSofUnionModels[[i]]
    xGETSofUnionModels[[i]] <- integer(0)
    
    ##blocks:
    if( !blocks.is.list ){

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

      ##
      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for(j in 1:no.of.blocks){
        if( blocksize*j <= ncol.adj ){
          partitions.t2[j] <- blocksize*j
        }
      }
      #check if last block contains last regressor:
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

    ##keep:
    ##do the mxkeep thing here!

    ##make blocks function for lapply/parLapply:
    XblocksFun <- function(j, i, x, blocks,
      parallel.options, y, user.estimator, t.pval, wald.pval, do.pet,
      ar.LjungB, arch.LjungB, normality.JarqueB, user.diagnostics,
      gof.function, gof.method, keep, include.gum, include.1cut,
      include.empty, max.paths, turbo, force.invertibility, tol,
      LAPACK, max.regs, print.searchinfo){

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
      getsx <- getsFun(y, mX, untransformed.residuals=NULL,
        user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, gof.function=gof.function,
        gof.method=gof.method, keep=keep, include.gum=include.gum,
        include.1cut=include.1cut, include.empty=include.empty,
        max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
        max.regs=max.regs, print.searchinfo=print.searchinfo,
        alarm=FALSE)

      if( is.null(getsx$specific.spec) ){
        xSpecificmodels <- NULL
      }else{
        xSpecificmodels <- names(getsx$specific.spec)
        ##For the future: add more info to return? E.g.
        ##  Xpaths[[j]] <- getsx$paths
        ##  Xterminals.results[[j]] <- getsx$terminals.results
      }

      ##return
      return(xSpecificmodels)

    } #close XblocksFun

    ##call XblocksFun/do gets on each block: NO parallel computing
    if( is.null(parallel.options) ){
      xSpecificmodels <- lapply(1:length(blocks[[i]]),
        XblocksFun, i, x, blocks, parallel.options,
        y, user.estimator, t.pval, wald.pval, do.pet, ar.LjungB,
        arch.LjungB, normality.JarqueB, user.diagnostics, gof.function,
        gof.method, keep, include.gum, include.1cut,
        include.empty, max.paths, turbo, force.invertibility, tol,
        LAPACK, max.regs, print.searchinfo)
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
        blocks, parallel.options, y, user.estimator, t.pval,
        wald.pval, do.pet, ar.LjungB, arch.LjungB, normality.JarqueB,
        user.diagnostics, gof.function, gof.method, keep,
        include.gum, include.1cut, include.empty, max.paths, turbo,
        force.invertibility, tol, LAPACK, max.regs, print.searchinfo)
      stopCluster(blocksClust)

    } #end if( parallel computing )

    ##gets of union of retained variables:
    ##------------------------------------
    if( print.searchinfo ){
      message("\n", appendLF=FALSE)
      message("GETS of union of retained ",
        names(x)[i], " variables... ",
        appendLF=TRUE)
      message("\n", appendLF=FALSE)
    }

    ##if no variables retained from the blocks:
    if( length(xSpecificmodels)==0 ){
      xNames <- NULL
#OLD:
#      xGETSofUnionModels[[i]] <- NULL
    }

    ##when indicators/variables(x) retained from the blocks:
    if( length(xSpecificmodels)>0 ){

      xNames <- NULL

      #which variables retained?:
      for(j in 1:length(xSpecificmodels)){
        #check if non-empty:
        if( !is.null(xSpecificmodels[[j]]) ){
          xNames <- union(xNames, xSpecificmodels[[j]])
        }
      } #end for(j) loop

      #do gets with union of retained x's:
#OLD:
#      if( length(xNames)==0 ){
#        xGETSofUnionModels[[i]] <- NULL
#      }
      if( length(xNames)>0 ){
        mX <- cbind(x[[i]][,xNames])
        colnames(mX) <- xNames
        if( force.invertibility ){
          mX <- dropvar(mX, tol=tol, LAPACK=LAPACK, silent=TRUE)
        }
        ##note: include.gum=TRUE is needed for the situation
        ##where all regressors are significant in the gum
        getsx <- getsFun(y, mX, untransformed.residuals=NULL,
          user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
          wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          user.diagnostics=user.diagnostics, gof.function=gof.function,
          gof.method=gof.method, keep=keep, include.gum=TRUE,
          include.1cut=include.1cut, include.empty=include.empty,
          max.paths=max.paths, turbo=turbo, tol=tol, LAPACK=LAPACK,
          max.regs=max.regs, print.searchinfo=print.searchinfo,
          alarm=FALSE)
        if( !is.null(names(getsx$specific.spec)) ){
          xGETSofUnionModels[[i]] <- names(getsx$specific.spec)
        }
      }

    } #end if(length(xSpecificmodels > 0)

  } #end for(i) loop (on x matrices)

  ##add names to blocks:
  names(blocks) <- names(x)
  names(xGETSofUnionModels) <- names(x)


  ##-------------------------------
  ## 4 make return object:
  ##-------------------------------

  result$y <- y
  result$x <- x
  result$blocks <- blocks
  result$best.specific <- xGETSofUnionModels
  result$time.finished <- date()
  if(alarm){ alarm() }
  return(result)
  
} #close blocksFun()
