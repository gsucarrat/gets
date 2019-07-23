\name{eqwma}
\alias{eqwma}
\alias{leqwma}

\title{Equally Weighted Moving Average (EqWMA) of the pth. exponentiated values}

\description{
The function \code{eqwma} returns an Equally Weighted Moving Average (EqWMA) of the pth. exponentiated values lagged. Optionally, the absolute values are computed before averaging, and the log of is returned. The function \code{leqwma} is essentially a wrapper to \code{eqwma} in which the absolute values are used and the logarithm is applied.

If x is financial return (possibly mean-corrected) and p=2, then this gives the socalled 'historical' model, also known as an integrated ARCH model where the ARCH coefficients all have the same value with sum equal to one. In the log-variance specification the lag of log(EqWMA) is thus a financial volatility proxy. It may be an imperfect proxy compared with high-frequency data (which can also be included as regressors), but - in contrast to high-frequency data - is always available and easy to compute.
}
\usage{
eqwma(x, length = 5, lag = 1, start = 1, p = 1, log = FALSE, abs = FALSE,
  as.vector = TRUE)
leqwma(x, length = 5, lag = 1, start = 1, p = 2, as.vector=FALSE)
}

\arguments{
  \item{x}{numeric vector, time-series or \code{zoo} object. Missing values in the beginning and/or at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{length}{integer or vector of integers each equal to or greater than 1. The length or lengths of the moving window or windows of averages}
  \item{lag}{integer equal to or greater than 0. If 0, then the moving averages are not lagged}
  \item{start}{integer equal to or greater than 1 (default: start=1, i.e. the first observation). Where to start the moving windows of averages}
  \item{p}{numeric value greater than zero. The exponent p in x^p for \code{eqwma} and in abs(x)^p for \code{leqwma} }
  \item{log}{logical. If TRUE, then the logarithm of the moving average is returned. If FALSE (default), then the logarithm is not applied}
  \item{abs}{logical. If TRUE, then x is transformed to absolute values before x is exponentiated}
  \item{as.vector}{logical. If TRUE, then a univariate series is returned as a vector. If FALSE, then a univariate series is returnes as a matrix. Note: multivariate series are always returned as a matrix}
}

\details{
The intended primary use of \code{eqwma} is to construct mixed frequency regressors for the mean specification.

The intended primary use of \code{leqwma} is to construct volatility proxies for the log-variance specification. The default is the lagged log of an equally weighted moving average of the squared residuals, where each average is made up of m observations. This is equivalent to an integrated ARCH(p) model where the p coefficients are all equal. For further details on the use of log(EqWMA) as a volatility proxy, see Sucarrat and Escribano (2012).
}

\value{numeric vector, time series or \code{\link{zoo}} object}

\references{
Genaro Sucarrat and Alvaro Escribano (2012): 'Automated Financial
  Model Selection: General-to-Specific Modelling of the Mean and
  Volatility Specifications', Oxford Bulletin of Economics and
  Statistics 74, Issue no. 5 (October), pp. 716-735

Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{
\code{\link{zoo}}, \code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}
}

\examples{
##generate an iid normal series:
set.seed(123)
x <- rnorm(100)

##compute lag of EqWMA(20) for x^2:
eqwma(x, p=2)

##compute lag of EqWMA(5) and lag of EqWMA(10) for x:
eqwma(x, length=c(5,10))

##compute lag of log(EqWMA(20)) for x^2:
leqwma(x)

#compute lag of log(EqWMA(5)) and lag of log(EqWMA(8))
#for abs(x)^2:
leqwma(x, length=c(4,8))

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
