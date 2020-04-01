\name{eqwma}
\alias{eqwma}
\alias{leqwma}

\title{Equally Weighted Moving Average (EqWMA) of the pth. exponentiated values}

\description{
The function \code{eqwma} returns an Equally Weighted Moving Average (EqWMA) of the pth. exponentiated values lagged \code{k} times (the default of \code{k} is 1). Optionally, the absolute values are computed before averaging if \code{abs=TRUE}, and the natural log of the values is returned if \code{log=TRUE}. The function \code{leqwma} is a wrapper to \code{eqwma} with \code{abs=TRUE} and \code{log=TRUE}.

If x is financial return (possibly mean-corrected) and p=2, then this gives the socalled 'historical' model, also known as an integrated ARCH model where the ARCH coefficients all have the same value with sum equal to one. In the log-variance specification the lag of log(EqWMA) is thus a financial volatility proxy. It may be an imperfect proxy compared with high-frequency data (which can also be included as regressors), but - in contrast to high-frequency data - is always available and easy to compute.
}
\usage{
eqwma(x, length=5, k=1, p=1, abs=FALSE, log=FALSE,  as.vector=FALSE,
  lag=NULL, start=NULL)
leqwma(x, length=5, k=1, p=2, as.vector=FALSE, lag=NULL, start=NULL)
}

\arguments{
  \item{x}{numeric vector, time-series or \code{zoo} object}
  \item{length}{integer or vector of integers each equal to or greater than 1. The length or lengths of the moving window or windows of averages}
  \item{k}{integer that determines how many periods the term(s) should be lagged. If 0 (or smaller), then the moving averages are not lagged}
  \item{p}{numeric value. The exponent p in x^p when \code{abs=FALSE}, and in abs(x)^p when \code{abs=TRUE} }
  \item{log}{logical with default \code{FALSE}. If \code{TRUE}, then the logarithm of the moving average is returned}
  \item{abs}{logical with default \code{FALSE}. If \code{TRUE}, then x is transformed to absolute values before x is exponentiated}
  \item{as.vector}{logical with default \code{FALSE}. If \code{TRUE}, and if \code{length(length)==1}, then the result is returned as a vector. Otherwise the returned value is always a matrix}
  \item{lag}{deprecated}
  \item{start}{deprecated}
}

\details{
The intended primary use of \code{eqwma} is to construct mixed frequency regressors for the mean specification of an \code{\link{arx}} model.

The intended primary use of \code{leqwma} is to construct volatility proxies for the log-variance specification in an \code{\link{arx}} model. In the latter case, the default is the lagged log of an equally weighted moving average of the squared residuals, where each average is made up of m observations. This is equivalent to an integrated ARCH(p) model where the p coefficients are all equal. For further details on the use of log(EqWMA) as a volatility proxy, see Sucarrat and Escribano (2012).
}

\value{numeric matrix, vector or \code{\link{zoo}} object}

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
