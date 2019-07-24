\name{recursive}
\alias{recursive}

\title{Recursive estimation}

\description{
Recursive estimation of coefficients and standard errors
}

\usage{
recursive(object, spec="mean", std.errors=TRUE, from=40, tol=1e-07,
  LAPACK=FALSE, plot=TRUE, return=TRUE)
}

\arguments{
  \item{object}{an \code{\link{arx}}, \code{gets} or \code{isat} object}
  \item{spec}{'mean' or 'variance'. If 'mean' (default), the the recursive estimates of the mean-equation are estimated}
  \item{std.errors}{logical. If TRUE (default), then the coefficient standard errors are also computed}
  \item{from}{integer. The starting point of the recursion}
  \item{tol}{numeric. The tolerance for linear dependency among regressors}
  \item{LAPACK}{logical, TRUE or FALSE (default). If true use LAPACK otherwise use LINPACK, see \code{\link{qr}} function}
  \item{plot}{NULL or logical. If TRUE, then the recursive coefficient estimates are plotted. If NULL (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
  \item{return}{logical. If TRUE (default), then the recursive estimates are returned in a list}
}

\value{If \code{return=TRUE}, then a \code{\link{list}} is returned with the following components:
  \item{estimates}{a \code{\link{zoo}} matrix with the recursive estimates}
  \item{standard.errors}{a \code{\link{zoo}} matrix with the standard errors}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{ols}}, \code{\link{qr}}, \code{\link{solve.qr}}}

\examples{
##generate random variates, estimate model:
y <- rnorm(100)
mX <- matrix(rnorm(4*100), 100, 4)
mymodel <- arx(y, mc=TRUE, mxreg=mX)

##compute recursive estimates and plot them:
recursive(mymodel)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}