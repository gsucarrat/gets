\name{dlogitxSim}
\alias{dlogitxSim}

\title{Simulate from a dynamic logit-x model}

\description{
Simulate from a dynamic Autoregressive (AR) logit model with covariates ('X').
}

\usage{
dlogitxSim(n, intercept = 0, ar = NULL, xreg = NULL, verbose = FALSE, 
    as.zoo = TRUE)
}

\arguments{
  \item{n}{integer, the number of observations to generate}
  \item{intercept}{numeric, the value of the intercept in the logit specification}
  \item{ar}{\code{NULL} or a numeric vector with the autoregressive parameters}
  \item{xreg}{\code{NULL} or numeric vector with the values of the X-term}
  \item{verbose}{\code{logical}. If \code{FALSE}, then only the binary process (a vector) is returned. If \code{TRUE}, then a matrix with all the simulated information is returned (binary process, probabilities, etc.)}
  \item{as.zoo}{\code{logical}. If \code{TRUE}, then the returned object - a vector or matrix - will be of class \code{\link{zoo}} }
}

\details{No details, for the moment.}

\value{
A vector or matrix, depending on whether \code{verbose} is \code{FALSE} or \code{TRUE}, of class \code{\link{zoo}}, depending on whether \code{as.zoo} is \code{TRUE} or \code{FALSE}
}

\references{No references for the moment.}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{dlogitx}}
}

\examples{
##no examples for the moment
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
