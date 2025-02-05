\name{logitxSim}
\alias{logitxSim}
\alias{dlogitxSim}

\title{Simulate from a dynamic logit-x model}

\description{
Simulate from a dynamic Autoregressive (AR) logit model with covariates ('X'). This model is essentially a logit-version of the model of Kauppi and Saikkonen (2008).
}

\usage{
logitxSim(n, intercept = 0, ar = NULL, xreg = NULL, verbose = FALSE, 
    as.zoo = TRUE)

dlogitxSim(n, ...)
}

\arguments{
  \item{n}{integer, the number of observations to generate}
  \item{intercept}{numeric, the value of the intercept in the logit specification}
  \item{ar}{\code{NULL} or a numeric vector with the autoregressive parameters}
  \item{xreg}{\code{NULL} or numeric vector with the values of the X-term}
  \item{verbose}{\code{logical}. If \code{FALSE}, then only the binary process (a vector) is returned. If \code{TRUE}, then a matrix with all the simulated information is returned (binary process, probabilities, etc.)}
  \item{as.zoo}{\code{logical}. If \code{TRUE}, then the returned object - a vector or matrix - will be of class \code{\link[zoo]{zoo}} }
  \item{...}{arguments passed on to \code{logitxSim}}
}

\details{No details, for the moment.}

\value{
A vector or matrix, depending on whether \code{verbose} is \code{FALSE} or \code{TRUE}, of class \code{\link[zoo]{zoo}}, depending on whether \code{as.zoo} is \code{TRUE} or \code{FALSE}
}

\references{
Heikki Kauppi and Penti Saikkonen (2008): 'Predicting U.S. Recessions with Dynamic Binary Response Models'. The Review of Economic Statistics 90, pp. 777-791
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{logitx}}
}

\examples{

##simulate from ar(1):
set.seed(123) #for reproducibility
y <- logitxSim(100, ar=0.3)

##more output (value, probability, logit):
set.seed(123) #for reproducibility
y <- logitxSim(100, ar=0.3, verbose=TRUE)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
