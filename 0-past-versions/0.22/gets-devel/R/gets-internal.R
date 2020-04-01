.onAttach <- function(libname, pkgname)
{
  ##set start-up message:
  txt <- c("\n",
    paste(sQuote("gets"), "version 0.22\n"),
    "\n",
    paste0("General-to-Specific (GETS) and Indicator Saturation (ISAT) methods, type help(", sQuote("gets-package"), ") for details"),
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
} #close .onAttach
