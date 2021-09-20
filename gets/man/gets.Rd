\name{gets}
\alias{gets}
\alias{gets.arx}
\alias{gets.lm}
\alias{gets.isat}

\title{General-to-Specific (GETS) Modelling}

\description{Generic function that enables new GETS and ISAT methods for new classes}

\usage{
gets(x, ...)

##wrapper to getsm() and getsv():
\method{gets}{arx}(x, spec=NULL, ...)

##wrapper to lm():
\method{gets}{lm}(x, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
  gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
  include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
  turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE,...)

##wrapper to isat():
\method{gets}{isat}(x, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
  gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
  include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
  turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE,...)

}



\arguments{
  \item{x}{an object of class 'arx' used to select a method}
  \item{spec}{\code{NULL} (default), \code{"mean"} or \code{"variance"}. If \code{"mean"}, then \code{\link{getsm}} is called. If \code{"variance"}, then \code{\link{getsv}} is called. If \code{NULL}, then it is automatically determined whether GETS-modelling of the mean or log-variance specification should be undertaken.}
  \item{...}{further arguments passed to or from other methods}
}

\details{
\code{gets.arx} is a convenience wrapper to \code{\link{getsm}} and \code{\link{getsv}}. 
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{getsm}}, \code{\link{getsv}}, \code{\link{getsFun}}
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
