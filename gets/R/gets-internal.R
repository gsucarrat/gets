.onAttach <- function(libname, pkgname)
{
  ##set start-up message:
  txt <- c("\n",
    paste(sQuote("gets"), "version 0.34\n"),
    "\n",
    paste("General-to-Specific (GETS) and Indicator Saturation (ISAT) methods"),
    "\n",
    paste("CRAN website: https://CRAN.R-project.org/package=gets"),
    paste("Github (issues, discussions): https://github.com/gsucarrat/gets"),
    "\n",
    paste("* Type help(", dQuote("gets-package", q=FALSE), ") for details and browseVignettes(", dQuote("gets", q=FALSE), ") for an introduction", sep=""),
    "\n",
    paste("* For automatic plotting, set plot = TRUE in options: options(plot = TRUE)"),
    "\n",
    paste("* WARNING: New default 'mc = TRUE' in arx() as of version 0.28"),
    "\n")
  
  ##print message:
  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }
} #close .onAttach
