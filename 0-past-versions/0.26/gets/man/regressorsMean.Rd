\name{regressorsMean}
\alias{regressorsMean}

\title{Create the regressors of the mean equation}

\description{The function generates the regressors of the mean equation in an \code{\link{arx}} model. The returned value is a \code{matrix} with the regressors and, by default, the regressand in column one. By default, observations (rows) with missing values are removed in the beginning and the end with \code{\link{na.trim}}, and the returned matrix is a \code{\link{zoo}} object. }

\usage{
regressorsMean(y, mc = FALSE, ar = NULL, ewma = NULL, mxreg = NULL,
  return.regressand = TRUE, return.as.zoo = TRUE, na.trim = TRUE,
  na.omit=FALSE)
}

\arguments{
  \item{y}{numeric vector, time-series or \code{\link{zoo}} object. }
  \item{mc}{logical. \code{TRUE} includes an intercept, whereas \code{FALSE} (default) does not. }
  \item{ar}{either \code{NULL} (default) or an integer vector, say, \code{c(2,4)} or \code{1:4} with the AR-lags to include in the mean specification. If \code{NULL}, then no lags are included. }
  \item{ewma}{either \code{NULL} (default) or a \code{\link{list}} with arguments sent to the \code{\link{eqwma}} function. In the latter case a lagged moving average of \code{y} is included as a regressor. }
  \item{mxreg}{either \code{NULL} (default), numeric vector or matrix, say, a \code{\link{zoo}} object, or \code{\link{data.frame}} containing conditioning variables (covariates). Note that, if both \code{y} and \code{mxreg} are \code{zoo} objects, then their samples are matched. }
  \item{return.regressand}{logical. \code{TRUE}, the default, includes the regressand as column one in the returned matrix. }
  \item{return.as.zoo}{\code{TRUE}, the default, returns the matrix as a \code{\link{zoo}} object. }
  \item{na.trim}{\code{TRUE}, the default, removes observations with \code{NA}-values in the beginning and the end with \code{\link{na.trim}}. }
  \item{na.omit}{\code{TRUE}, the non-default, removes observations with \code{NA}-values, not necessarily in the beginning or in the end, with \code{\link{na.omit}}. }
}

\value{A matrix, by default of class \code{\link{zoo}}, with the regressand as column one (the default). }

\references{
Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. DOI: https://www.jstatsoft.org/article/view/v086i03
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{arx}}, \code{\link{isat}}, \code{\link{regressorsVariance}}, \code{\link{zoo}}, \code{\link{eqwma}}, \code{\link{na.trim}} and \code{\link{na.trim}}. }

\examples{

##generate some data:
y <- rnorm(10) #regressand
x <- matrix(rnorm(10*5), 10, 5) #regressors

##create regressors (examples):
regressorsMean(y, mxreg=x)
regressorsMean(y, mxreg=x, return.regressand=FALSE)
regressorsMean(y, mc=TRUE, ar=1:3, mxreg=x)
regressorsMean(log(y^2), mc=TRUE, ar=c(2,4))

##let y and x be time-series:
y <- ts(y, frequency=4, end=c(2018,4))
x <- ts(x, frequency=4, end=c(2018,4))
regressorsMean(y, mxreg=x)
regressorsMean(y, mc=TRUE, ar=1:3, mxreg=x)
regressorsMean(log(y^2), mc=TRUE, ar=c(2,4))

##missing values (NA):
y[1] <- NA
x[10,3] <- NA
regressorsMean(y, mxreg=x)
regressorsMean(y, mxreg=x, na.trim=FALSE)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}