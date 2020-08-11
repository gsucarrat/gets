\name{gets}
\alias{gets}
\alias{gets.arx}

\title{General-to-Specific (GETS) Modelling}

\description{Generic function that enables new GETS and ISAT methods for new classes}

\usage{
gets(x, ...)

##wrapper to getsm() and getsv():
\method{gets}{arx}(x, spec=NULL, ...)

}

\arguments{
  \item{x}{an object of class 'arx' used to select a method}
  \item{...}{further arguments passed to or from other methods}
  \item{spec}{\code{NULL} (default), \code{"mean"} or \code{"variance"}. If \code{"mean"}, then \code{\link{getsm}} is called. If \code{"variance"}, then \code{\link{getsv}} is called. If \code{NULL}, then it is automatically determined whether GETS-modelling of the mean or log-variance specification should be undertaken.}
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
