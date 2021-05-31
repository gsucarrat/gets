\name{as.arx}
\alias{as.arx}

\title{Convert to 'arx' object}

\description{Convert 'gets'/'isat'/'lm' object to 'arx' object}

\usage{as.arx(object, ar=NULL, ewma=NULL, vc=FALSE,
  arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL, zero.adj=0.1,
  vc.adj=TRUE, vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
  user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, plot=NULL)}

\arguments{
  \item{object}{object of class \code{\link{arx}}, \code{\link{gets}} or \code{\link{isat}} }
  \item{ar}{either \code{NULL} (default) or an integer vector, say, \code{c(2,4)} or \code{1:4}. The AR-lags   to include in the mean specification. If \code{NULL}, then no lags are included}
  \item{ewma}{either \code{NULL} (default) or a \code{\link{list}} with arguments sent to the \code{\link{eqwma}} function. In the latter case a lagged moving average of \code{y} is included as a regressor}
  \item{vc}{logical. \code{TRUE} includes an intercept in the log-variance specification, whereas \code{FALSE} (default) does not. If the log-variance specification contains any other item but the log-variance intercept, then vc is set to \code{TRUE} }
  \item{arch}{either \code{NULL} (default) or an integer vector, say, \code{c(1,3)} or \code{2:5}. The log-ARCH lags to include in the log-variance specification}
  \item{asym}{either \code{NULL} (default) or an integer vector, say, \code{c(1)} or \code{1:3}. The asymmetry (i.e. 'leverage') terms to include in the log-variance specification}
  \item{log.ewma}{either \code{NULL} (default) or a vector of the lengths of the volatility proxies, see \code{\link{leqwma}} }
  \item{vxreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link{zoo}} object, of conditioning variables. If both \code{y} and \code{mxreg} are \code{zoo} objects, then their samples are chosen to match}
  \item{zero.adj}{numeric value between 0 and 1. The quantile adjustment for zero values. The default 0.1 means the zero residuals are replaced by the 10 percent quantile of the absolute residuals before taking the logarithm}
  \item{vc.adj}{logical. If \code{TRUE} (default), then the log-variance intercept is adjusted by the estimate of E[ln(z^2)], where z is the standardised error. This adjustment is needed for the conditional scale to be equal to the conditional standard deviation. If \code{FALSE}, then the log-variance intercept is not adjusted}
  \item{vcov.type}{character vector, "ordinary" (default), "white" or "newey-west". If "ordinary", then the ordinary variance-covariance matrix is used for inference. If "white", then the White (1980) heteroscedasticity-robust matrix is used. If "newey-west", then the Newey and West (1987) heteroscedasticity and autocorrelation-robust matrix is used}
  \item{qstat.options}{\code{NULL} (default) or an integer vector of length two, say, \code{c(1,1)}. The first value sets the lag-order of the AR diagnostic test, whereas the second value sets the lag-order of the ARCH diagnostic test. If \code{NULL}, then the two values of the vector are set automatically}
  \item{normality.JarqueB}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, then the results of the Jarque and Bera (1980) test for non-normality in the residuals are included in the estimation results. }
  \item{user.estimator}{\code{NULL} (default) or a \code{\link{list}} with one entry, \code{name}, containing the name of the user-defined estimator. Additional items, if any, are passed on as arguments to the estimator in question }
  \item{user.diagnostics}{\code{NULL} (default) or a \code{\link{list}} with two entries, \code{name} and \code{pval}, see the \code{user.fun} argument in \code{\link{diagnostics}} }
  \item{tol}{numeric value (\code{default = 1e-07}). The tolerance for detecting linear dependencies in the columns of the regressors (see \code{\link{qr}} function). Only used if \code{LAPACK} is \code{FALSE} (default) and \code{user.estimator} is \code{NULL}.}
  \item{LAPACK}{logical. If \code{TRUE}, then use LAPACK. If \code{FALSE} (default), then use LINPACK (see \code{\link{qr}} function). Only used if \code{user.estimator} is \code{NULL}.}
  \item{plot}{\code{NULL} or logical. If \code{TRUE}, then the fitted values and the residuals are plotted. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
}

\value{Object of class \code{\link{arx}} }

\author{
Moritz Schwarz, \url{https://www.inet.ox.ac.uk/people/moritz-schwarz/}\cr
}

\seealso{
\code{\link{arx}}, \code{\link{gets}}, \code{\link{isat}}, \code{\link{lm}}
}

\examples{

lm_object <- lm(mpg ~ cyl + hp, data = mtcars)
as.arx(lm_object)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}