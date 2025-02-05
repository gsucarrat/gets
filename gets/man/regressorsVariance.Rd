\name{regressorsVariance}
\alias{regressorsVariance}

\title{Create regressors for a log-variance model}

\description{The function creates the regressors of a log-variance model, e.g. in a \code{\link{arx}} model. The returned value is a \code{matrix} with the regressors and, by default, the regressand in the first column. By default, observations (rows) with missing values are removed in the beginning and the end with \code{\link[zoo]{na.trim}}, and the returned matrix is a \code{\link[zoo]{zoo}} object. }

\usage{
regressorsVariance(e, vc = TRUE, arch = NULL, harch = NULL, asym = NULL,
  asymind = NULL, log.ewma = NULL, vxreg = NULL, prefix = "v", zero.adj = NULL,
  vc.adj = TRUE, return.regressand = TRUE, return.as.zoo = TRUE, na.trim = TRUE,
  na.omit = FALSE) 
}

\arguments{
  \item{e}{numeric vector, time-series or \code{\link[zoo]{zoo}} object. }
  \item{vc}{logical. \code{TRUE} includes an intercept in the log-variance specification, whereas \code{FALSE} (default) does not. If the log-variance specification contains any other item but the log-variance intercept, then vc is set to \code{TRUE}. }
  \item{arch}{either \code{NULL} (default) or an integer vector, say, \code{c(1,3)} or \code{2:5}. The log-ARCH lags to include in the log-variance specification. }
  \item{harch}{either \code{NULL} (default) or an integer vector, say, \code{c(5,20)}. The log of heterogenous ARCH-terms as proposed by Muller et al. (1997). }
  \item{asym}{either \code{NULL} (default) or an integer vector, say, \code{c(1)} or \code{1:3}. The asymmetry (i.e. 'leverage') terms to include in the log-variance specification. }
  \item{asymind}{either \code{NULL} (default) or an integer vector, say, \code{c(1)} or \code{1:3}. The indicator ('binary') asymmetry terms to include in the log-variance specification. }
  \item{log.ewma}{either \code{NULL} (default) or a vector of the lengths of the volatility proxies, see \code{\link{leqwma}}. The log of heterogenous volatility proxies similar to those of Corsi (2009). }
  \item{vxreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link[zoo]{zoo}} object, of conditioning variables. If both \code{y} and \code{mxreg} are \code{\link[zoo]{zoo}} objects, then their samples are chosen to match. }
  \item{prefix}{a \code{character} used as prefix in the labelling of the variables in \code{vxreg} and of the intercept. } 
  \item{zero.adj}{\code{NULL} (default) or a strictly positive \code{numeric} scalar. If \code{NULL}, the zeros in the squared e's are replaced by the 10 percent quantile of the non-zero squared e's. If \code{zero.adj} is a strictly positive \code{numeric} scalar, then this value is used to replace the zeros of the squared e's. }
  \item{vc.adj}{deprecated and ignored. }
  \item{return.regressand}{\code{logical}. \code{TRUE} (default) includes the regressand as column one in the returned matrix. }
  \item{return.as.zoo}{\code{logical}. \code{TRUE} (default) returns the matrix as a \code{\link[zoo]{zoo}} object. }
  \item{na.trim}{\code{logical}. \code{TRUE} (default) removes observations with \code{NA}-values in the beginning and the end with \code{\link[zoo]{na.trim}}. }
  \item{na.omit}{\code{logical}. \code{FALSE} (default) means \code{NA}-observations that are not in the beginning or at the end are kept (i.e. not omitted). \code{TRUE} removes with \code{\link{na.omit}}. }
}

\value{
  A \code{matrix}, by default of class \code{\link[zoo]{zoo}}, with the regressand as column one (the default).
}

\references{
Corsi, Fulvio (2009): 'A Simple Approximate Long-Memory Model of Realized Volatility', Journal of Financial Econometrics 7, pp. 174-196

Muller, Ulrich A., Dacorogna, Michel M., Dave, Rakhal D., Olsen, Richard B, Pictet, Olivier, Weizsaker, Jacob E. (1997): 'Volatilities of different time resolutions - Analyzing the dynamics of market components'. Journal of Empirical Finance 4, pp. 213-239

Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. DOI: https://www.jstatsoft.org/article/view/v086i03

Sucarrat, Genaro and Escribano, Alvaro (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{regressorsMean}}, \code{\link{arx}}, \code{\link[zoo]{zoo}}, \code{\link{leqwma}}, \code{\link[zoo]{na.trim}} and \code{\link{na.omit}}. }

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

##example where eps and x are time-series:
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