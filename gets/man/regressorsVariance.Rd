\name{regressorsVariance}
\alias{regressorsVariance}

\title{Create the regressors of the variance equation}

\description{The function generates the regressors of the log-variance equation in an \code{\link{arx}} model. The returned value is a \code{matrix} with the regressors and, by default, the regressand in column one. By default, observations (rows) with missing values are removed in the beginning and the end with \code{\link{na.trim}}, and the returned matrix is a \code{\link{zoo}} object. }

\usage{
regressorsVariance(e, vc = TRUE, arch = NULL, asym = NULL,
  log.ewma = NULL, vxreg = NULL, zero.adj = 0.1, vc.adj = TRUE,
  return.regressand = TRUE, return.as.zoo = TRUE, na.trim = TRUE,
  na.omit = FALSE)
}

\arguments{
  \item{e}{numeric vector, time-series or \code{\link{zoo}} object.}
  \item{vc}{logical. \code{TRUE} includes an intercept in the log-variance specification, whereas \code{FALSE} (default) does not. If the log-variance specification contains any other item but the log-variance intercept, then vc is set to \code{TRUE} }
  \item{arch}{either \code{NULL} (default) or an integer vector, say, \code{c(1,3)} or \code{2:5}. The log-ARCH lags to include in the log-variance specification}
  \item{asym}{either \code{NULL} (default) or an integer vector, say, \code{c(1)} or \code{1:3}. The asymmetry (i.e. 'leverage') terms to include in the log-variance specification}
  \item{log.ewma}{either \code{NULL} (default) or a vector of the lengths of the volatility proxies, see \code{\link{leqwma}} }
  \item{vxreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link{zoo}} object, of conditioning variables. If both \code{y} and \code{mxreg} are \code{zoo} objects, then their samples are chosen to match}
  \item{zero.adj}{numeric value between 0 and 1. The quantile adjustment for zero values. The default 0.1 means the zero residuals are replaced by the 10 percent quantile of the absolute residuals before taking the logarithm}
  \item{vc.adj}{logical. If \code{TRUE} (default), then the log-variance intercept is adjusted by the estimate of E[ln(z^2)], where z is the standardised error. This adjustment is needed for the conditional scale to be equal to the conditional standard deviation. If \code{FALSE}, then the log-variance intercept is not adjusted}
  \item{return.regressand}{logical. \code{TRUE}, the default, includes the regressand as column one in the returned matrix. }
  \item{return.as.zoo}{\code{TRUE}, the default, returns the matrix as a \code{\link{zoo}} object. }
  \item{na.trim}{\code{TRUE}, the default, removes observations with \code{NA}-values in the beginning and the end with \code{\link{na.trim}}. }
  \item{na.omit}{\code{TRUE}, the non-default, removes observations with \code{NA}-values, not necessarily in the beginning or in the end, with \code{\link{na.omit}}. }
}

\value{A matrix, by default of class \code{\link{zoo}}, with the regressand as column one (the default). }

\references{
Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. DOI: https://www.jstatsoft.org/article/view/v086i03

Sucarrat, Genaro and Escribano, Alvaro (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{regressorsMean}}, \code{\link{arx}}, \code{\link{zoo}}, \code{\link{leqwma}}, \code{\link{na.trim}} and \code{\link{na.omit}}. }

\examples{

##generate some data:
eps <- rnorm(10) #error term
x <- matrix(rnorm(10*5), 10, 5) #regressors

##create regressors (examples):
regressorsVariance(eps, vxreg=x)
regressorsVariance(eps, vxreg=x, return.regressand=FALSE)
regressorsVariance(eps, arch=1:3, vxreg=x)
regressorsVariance(eps, arch=1:2, asym=1, vxreg=x)
regressorsVariance(eps, arch=1:2, asym=1, log.ewma=5)

##let eps and x be time-series:
eps <- ts(eps, frequency=4, end=c(2018,4))
x <- ts(x, frequency=4, end=c(2018,4))
regressorsVariance(eps, vxreg=x)
regressorsVariance(eps, arch=1:3, vxreg=x)
regressorsVariance(eps, arch=1:2, asym=1, vxreg=x)
regressorsVariance(eps, arch=1:2, asym=1, log.ewma=5)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}