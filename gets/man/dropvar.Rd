\name{dropvar}
\alias{dropvar}

\title{Drop variable}

\description{Drops columns in a matrix to avoid perfect multicollinearity.}

\usage{
dropvar(x, tol=1e-07, LAPACK=FALSE, silent=FALSE)
}

\arguments{
  \item{x}{a matrix, possibly less than full column rank.}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies among regressors, see \code{\link{qr}} function. Only used if LAPACK is FALSE}
  \item{LAPACK}{logical, TRUE or FALSE (default). If true use LAPACK otherwise use LINPACK, see \code{\link{qr}} function}
  \item{silent}{logical, TRUE (default) or FALSE. Whether to print a notification whenever a regressor is removed}
}

\value{a matrix whose regressors linearly independent}

\details{
Original function \code{drop.coef} developed by Rune Haubo B. Christensen in package \code{ordinal}, \url{https://cran.r-project.org/package=ordinal}.
}

\references{
Rune H.B. Christensen (2014): 'ordinal: Regression Models for Ordinal Data'. \url{https://cran.r-project.org/package=ordinal}
}

\author{
  Rune Haubo B. Christensen, with modifications by Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{isat}}}

\examples{
set.seed(1)
x <- matrix(rnorm(20), 5)
dropvar(x) #full rank, none are dropped

x[,4] <- x[,1]*2
dropvar(x) #less than full rank, last column dropped
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
