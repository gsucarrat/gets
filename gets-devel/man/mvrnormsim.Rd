\name{mvrnormsim}
\alias{mvrnormsim}

\title{Simulate from a Multivariate Normal Distribution}

\description{
  Produces one or more samples from the specified multivariate normal distribution. Used
  in \code{\link{outlierscaletest}}.
}

\usage{
mvrnormsim(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE)
}

\arguments{
  \item{n}{the number of samples required.}
  \item{mu}{a vector giving the means of the variables.}
  \item{Sigma}{a positive-definite symmetric matrix specifying the covariance matrix of the variables.}
  \item{tol}{tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.}
  \item{empirical}{logical. If true, mu and Sigma specify the empirical not population mean and covariance matrix.}
}

\value{If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.}

\details{
Original function \code{mvrnorm} developed by Venables, W. N. & Ripley. in package \code{MASS}, \url{https://CRAN.R-project.org/package=MASS}.
}

\references{
 Venables, W. N. & Ripley, B. D. (2019): 'MASS: Support Functions and Datasets for Venables and Ripley's MASS'. \url{https://CRAN.R-project.org/package=MASS}
 
 Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
 
}

\author{
  Venables, W. N. & Ripley, with modifications by Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{\code{\link{outlierscaletest}}}

\examples{
Sigma <- matrix(c(3,2,1,7),2,2)
mvrnormsim(n=2, mu=c(1,2), Sigma)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
