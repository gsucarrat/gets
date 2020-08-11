\name{iim}
\alias{iim}
\alias{sim}
\alias{tim}

\title{Make Indicator Matrices (Impulses, Steps, Trends)}

\description{
Auxiliary functions to make, respectively, matrices of impulse indicators (\code{iim}), step indicators (\code{sim}) and trend indicators (\code{tim})
}

\usage{
##make matrix of impulse indicators:
iim(x, which.ones = NULL)

##make matrix of step indicators:
sim(x, which.ones = NULL)

##make matrix of trend indicators:
tim(x, which.ones = NULL, log.trend = FALSE)
}

\arguments{
  \item{x}{either an integer (the length of the series in question) or a series (a vector or matrix) from which to use the time-series index to make indicators of}
  \item{which.ones}{the locations of the impulses. If NULL (the default), then all impulses are returned}
  \item{log.trend}{logical. If TRUE, then the natural log is applied on the trends}
}

\details{
If \code{x} is a series or vector of observations, then the index of \code{x} will be used for the labelling of the impulses, and in the returned \code{\link{zoo}} object.

Note: For \code{sim} and \code{tim} the first indicator is removed, since it is exactly colinear with the others.
}

\value{
A \code{\link{zoo}} matrix containing the impulses
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{
\code{\link{isat}}, \code{\link{zoo}}
}
\examples{
##generate series:
y <- rnorm(40)

##make matrix of impulse indicators:
mIIM <- iim(40)

##make matrix of step-indicators, but only every third:
mSIM <- sim(y, which.ones=seq(1,40,3))

##give quarterly time-series attributes to y-series:
y <- zooreg(y, frequency=4, end=c(2015,4))

##make matrix of trend-indicators with quarterly labels:
mTIM <- tim(y)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
