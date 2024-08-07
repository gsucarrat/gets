\name{ols}
\alias{ols}

\title{OLS estimation}

\description{
OLS estimation with the QR decomposition and, for some options, computation of variance-covariance matrices
}

\usage{
ols(y, x, untransformed.residuals=NULL, tol=1e-07, LAPACK=FALSE, method=3, 
  variance.spec=NULL, ...)
}

\arguments{
  \item{y}{numeric vector, the regressand}
  \item{x}{numeric matrix, the regressors}
  \item{untransformed.residuals}{\code{NULL} (default) or, when \code{\link{ols}} is used with \code{method=6}, a numeric vector containing the untransformed residuals}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the regressors, see the \code{\link{.lm.fit}} function}
  \item{LAPACK}{deprecated and ignored}
  \item{method}{an integer, 1 to 6, that determines the estimation method}
  \item{variance.spec}{\code{NULL} or a \code{\link{list}} with items that specifies the log-variance model to be estimated, see \code{\link{arx}} }
  \item{...}{further arguments (currently ignored)}
}

\details{
\code{method = 1} or \code{method = 2} only returns the OLS coefficient estimates together with the QR- information, the former being slightly faster. \code{method=3} returns, in addition, the ordinary variance-covariance matrix of the OLS estimator. \code{method=4} returns the White (1980) heteroscedasticity robust variance-covariance matrix in addition to the information returned by \code{method=3}, whereas \code{method=5} does the same except that the variance-covariance matrix now is that of Newey and West (1987). \code{method=6} undertakes OLS estimation of a log-variance model, see Pretis, Reade and Sucarrat (2018, Section 4). Alternatively, for \code{method} 1 to 5, a log-variance model is also estimated if \code{variance.spec} is not \code{NULL}.
}

\value{A list with items depending on \code{method}}

\references{
W. Newey and K. West (1987): 'A Simple Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708.

F. Pretis, J. Reade and G. Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks', Journal of Statistical Software 86, Issue 3, pp. 1-44, DOI: https://doi.org/10.18637/jss.v086.i03

H. White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{.lm.fit}}, \code{\link{qr}}, \code{\link{solve.qr}}, \code{\link{arx}} }

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
