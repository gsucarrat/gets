\name{as.lm}
\alias{as.lm}

\title{Convert to 'lm' object}

\description{Convert 'arx'/'gets'/'isat' object to 'lm' object}

\usage{as.lm(object)}

\arguments{
  \item{object}{object of class \code{\link{arx}}, \code{\link{gets}} or \code{\link{isat}} }
}

\value{Object of class \code{\link{lm}} }

\author{
Moritz Schwarz
Genaro Sucarrat
}

\seealso{
\code{\link{arx}}, \code{\link{gets}}, \code{\link{isat}}, \code{\link{lm}}
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