\name{gets}
\alias{gets}
\alias{gets.arx}

\title{General-to-Specific (GETS) Modelling}

\description{
For an overview of the \bold{gets} package, see \code{\link{gets-package}}. Here, documentation of generic functions for GETS modelling is provided. Note that \code{gets.arx} is a convenience wrapper to \code{\link{getsm}} and \code{\link{getsv}}. For specific GETS methods for \code{lm}, \code{logitx} and \code{isat} models, see \code{\link{gets.lm}}, \code{\link{gets.logitx}} and \code{\link{gets.isat}}, respectively. 
}

\usage{
gets(x, ...)

\method{gets}{arx}(x, spec=NULL, ...)

}

\arguments{
  \item{x}{an object to be subjected to GETS modelling}
  \item{spec}{\code{NULL} (default), \code{"mean"} or \code{"variance"}. If \code{"mean"}, then \code{\link{getsm}} is called. If \code{"variance"}, then \code{\link{getsv}} is called. If \code{NULL}, then it is automatically determined whether GETS-modelling of the mean or log-variance specification should be undertaken.}
  \item{...}{further arguments passed to or from other methods}
}

\details{
\code{gets.arx} is a convenience wrapper to \code{\link{getsm}} and \code{\link{getsv}}. 
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{getsm}}, \code{\link{getsv}}, \code{\link{getsFun}}, \code{\link{blocksFun}}
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
