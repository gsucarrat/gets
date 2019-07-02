\name{periodicdummies}
\alias{periodicdummies}
\title{Make matrix of periodicity (e.g. seasonal) dummies}

\description{
Auxiliary function that generates periodicity dummies (e.g. seasonal dummies) for regular time series. The function is similar to, but more general than, the \code{seasonaldummy} function in the forecast package.
}

\usage{
periodicdummies(x, values=1)
}
\arguments{
  \item{x}{a regular time series (vector or matrix)}
  \item{values}{numeric of length 1 (default) or numeric vector of length equal to frequency(x) }
}

\value{A matrix of class zoo with the periodicity dummies}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{is.regular}}, \code{\link{zooreg}}, \code{\link{zoo}}, \code{\link{ts}} }

\examples{
##quarterly dummies:
x <- zooreg(rnorm(30), start=2000, frequency=4)
periodicdummies(x)

##monthly dummies:
y <- zooreg(rnorm(30), start=c(2000,1), frequency=12)
periodicdummies(y)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
