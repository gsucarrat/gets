\name{eviews}
\alias{eviews}
\alias{stata}
\title{Exporting results to EViews and STATA}

\description{
Functions that facilitate the export of results to the commercial econometric softwares EViews and STATA, respectively.
}

\usage{
eviews(object, file=NULL, print=TRUE, return=FALSE)
stata(object, file=NULL, print=TRUE, return=FALSE)
}
\arguments{
  \item{object}{an \code{arx}, \code{gets} or \code{isat} object}
  \item{file}{filename, i.e. the destination of the exported data}
  \item{print}{logical. If TRUE, then the estimation code in EViews (or STATA) is printed}
  \item{return}{logical. If TRUE, then a list is returned}
}

\value{Either printed text or a \code{list} (if return=TRUE)}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}}

\examples{
##simulate random variates, estimate model:
y <- rnorm(30)
mX <- matrix(rnorm(30*2), 30, 2)
mymod <- arx(y, mc=TRUE, mxreg=mX)

##print EViews code:
eviews(mymod)

##print Stata code:
stata(mymod)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
