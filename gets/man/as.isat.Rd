\name{as.isat}
\alias{as.isat}

\title{Convert to 'isat' object}

\description{Experimental function to convert 'arx'/'gets' object to 'isat' object. Note that this conversion will naturally differ from a genuine 'isat' object, especially with regard to the list elements of the object that describe the saturation algorithm (paths, best.terminal, no.of.estimations, etc.). So please use with care.}

\usage{as.isat(object, indicator_regex = list(iis = "^iis", sis = "^sis", tis = "^tis", uis = "^uis"), ...)}

\arguments{
  \item{object}{object of class \code{\link{arx}} or \code{\link{gets}}}
  \item{indicator_regex}{List of regular expressions to identify the indicators in the object. The default is \code{list(iis = "^iis", sis = "^sis", tis = "^tis", uis = "^uis")}.}
  \item{...}{Additional arguments passed to \code{\link{isat}}.}
}

\value{Object of class \code{\link{isat}} }

\author{
Moritz Schwarz, \url{https://www.moritzschwarz.org/}\cr\cr
}

\seealso{
\code{\link{arx}}, \code{\link{gets}}, \code{\link{isat}}
}

\examples{
##generate data, estimate model of class 'arx':
set.seed(123)
y <- rnorm(30)
arxmod <- arx(y, mc=TRUE, ar=1:3)
as.lm(arxmod)

##from 'gets' to 'lm':
getsmod <- getsm(arxmod, keep=1)
as.lm(getsmod)

##from 'isat' to 'lm':
isatmod <- isat(y)
as.lm(isatmod)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}