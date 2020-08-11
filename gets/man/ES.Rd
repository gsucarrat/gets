\name{ES}
\alias{ES}
\alias{VaR}

\title{Conditional Value-at-Risk (VaR) and Expected Shortfall (ES)}

\description{
Extract the in-sample conditional Value-at-Risk, or the in-sample conditional Expected Shortfall for the chosen risk level(s).
}

\usage{
ES(object, level=0.99, type=7, ...)
VaR(object, level=0.99, type=7, ...)
}

\arguments{
  \item{object}{an \code{\link{arx}} or \code{gets} object}
  \item{level}{the risk level(s), must be between 0 and 1}
  \item{type}{the method used to compute the empirical quantiles of the standardised residuals}
  \item{...}{arguments passed on (currently not used)}
}

\value{
A vector or matrix containing either the conditional Value-at-Risk (VaR) or the conditional Expected Shortfall (ES) for the chosen risk level(s).
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}} }

\examples{
##generate random variates, estimate model:
y <- rnorm(50)
mymodel <- arx(y, arch=1)

##extract 99% expected shortfall:
ES(mymodel)

##extract 99%, 95% and 90% expected shortfalls:
ES(mymodel, level=c(0.99, 0.95, 0.9))

##extract 99% value-at-risk:
VaR(mymodel)

##extract 99%, 95% and 90% values-at-risk:
VaR(mymodel, level=c(0.99, 0.95, 0.9))
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}