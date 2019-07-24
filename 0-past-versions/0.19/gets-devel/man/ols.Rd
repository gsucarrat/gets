\name{ols}
\alias{ols}

\title{OLS estimation}

\description{
OLS estimation with the QR decomposition and, for some options, computation of variance-
covariance matrices
}

\usage{
ols(y, x, untransformed.residuals=NULL, tol=1e-07, LAPACK=FALSE, method=3,
  user.fun=NULL, user.options=NULL)
}

\arguments{
  \item{y}{numeric vector, the regressand}
  \item{x}{numeric matrix, the regressors}
  \item{untransformed.residuals}{\code{NULL}, the default, or - when \code{method=6} - a numeric vector containing the untransformed residuals}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the regressors, see the \code{\link{qr}} function. Only used if \code{LAPACK} is \code{FALSE}}
  \item{LAPACK}{logical, \code{TRUE} or \code{FALSE} (default). If true use LAPACK otherwise use LINPACK, see the \code{\link{qr}} function}
  \item{method}{an integer, 0 to 6, that determines the estimation method}
  \item{user.fun}{\code{NULL} (default) or the name (a character) of the user-defined estimator}
  \item{user.options}{\code{NULL} (default) or a list with arguments (entries) that are passed on to the user-defined function}
}

\details{
\code{method = 1} or 2 only returns the OLS coefficient estimates together with the QR-information. \code{method = 1} is slightly faster than \code{method=2}. \code{method=3} returns in addition the ordinary variance-covariance matrix of the OLS estimator. \code{method=4} returns the White (1980) heteroscedasticity robust variance-covariance matrix in addition to the information returned by \code{method=3}, whereas \code{method=5} does the same except that the variance-covariance matrix now is that of Newey and West (1987). \code{method=6} undertakes OLS estimation of a log-variance model.
}

\value{A list with items depending on the \code{method}}

\references{
H. White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838.

W. Newey and K. West (1987): 'A Simple Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{qr}}, \code{\link{solve.qr}}}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
