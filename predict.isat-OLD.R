##==================================================
## forecast up to n.ahead
predict.isat <- function(object, n.ahead=12,
  newmxreg=NULL, newindex=NULL, return=TRUE, plot=NULL,
  plot.options=list(), ...)
{

  ##create new object to add stuff to in order to use predict.arx()
  object.new <- object

  ##------------
  ## arguments:
  ##------------

  if("mX" %in% names(object$aux)) {

    ##check if constant is retained:
    if("mconst" %in% object$aux$mXnames){
      object.new$call$mc <- TRUE
    }else{
      object.new$call$mc <- NULL
    }

    ##what dynamics specified in gum?
    gum.ar <- eval(object$call$ar)
    ##what dynamics remain in specific?
    spec.ar <- as.numeric(gsub("ar(\\d+)","\\1",object$aux$mXnames[grep("^ar\\d+$",object$aux$mXnames)]))
    if(NROW(spec.ar)==0) {
      object.new$call$ar <- NULL
    } else { ##check that dynamics in specific are subset of those in gum
      object.new$call$ar <- spec.ar[spec.ar %in% gum.ar]
    }

    ##"mxreg" argument:
    mc.and.ar.length <- length(object.new$call$mc)+length(object.new$call$ar)
    if(NCOL(object$aux$mX) > mc.and.ar.length){
      object.new$call$mxreg <- "mxreg"
    }else{
      object.new$call$mxreg <- NULL
    }

    ##if sis and tis terms retained need to adapt mxreg call ...
    if(!is.null(object$ISnames)) {

      ##J-dog, add your code here??
      if(is.null(object$call$mxreg)) { #need to ensure predict.arx knows there are mx variables
        object.new$call$mxreg <- "mXis"
      }

      ##... and need to specify newmxreg of right dimension:
      ##if we're here it means the isat call did not specify any mxregs
      ##hence what is in mX in the object are terms retained by isat
      ##we can automatically create isat terms into sample period (exception uis)
      if(is.null(newmxreg)) {
        ##if no newmxreg specified we add iis/sis/tis from scratch

        ##first check that there shouldn't be something in newmxreg...
        if(!is.null(object$call$mxreg)){ stop("'newmxreg' is NULL") }

        ##assuming not, then we start from scratch adding the indicators...
        newmxreg <- c()
      }

      if(any(regexpr("^iis",object$ISnames)>-1)){##isat retained some iis terms
        for(i in object$ISnames[grep("^iis",object$ISnames)]) {
          newmxreg <- cbind(newmxreg,rep(0,n.ahead))
        }
      }
      if(any(regexpr("^sis",object$ISnames)>-1)){##isat retained some sis terms
        for(i in object$ISnames[grep("^sis",object$ISnames)]) {
          newmxreg <- cbind(newmxreg,rep(1,n.ahead))
        }
      }
      if(any(regexpr("^tis",object$ISnames)>-1)){##isat retained some tis terms
        for(i in object$ISnames[grep("^tis",object$ISnames)]) {
          newmxreg <- cbind(newmxreg,
                            seq(1,n.ahead)+object$aux$mX[NROW(object$aux$mX),i])
        }
      }

    }

  } else {

    object.new$call$mc <- NULL
    object.new$call$ar <- NULL
    ##to do: log.ewma
    object.new$call$mxreg <- NULL

  }


  ##-----------------------------------
  ## pass on arguments to predict.arx:
  ##-----------------------------------

  out <- predict.arx(object.new, spec="mean", n.ahead=n.ahead,
    newmxreg=newmxreg, newvxreg=NULL, newindex=newindex,
    return=return, plot=plot, plot.options=plot.options)

  ##-------------------
  ## return forecasts:
  ##-------------------

  if(return){ return(out) }

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
