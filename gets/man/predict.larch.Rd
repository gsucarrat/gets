\name{predict.larch}
\alias{predict.larch}

\title{Variance forecasting with 'larch' models}

\description{Generate out-of-sample variance forecasts up to \code{n.ahead} steps ahead. Optionally, quantiles of the forecasts are also returned if the argument \code{probs} is specified. The forecasts, confidence intervals and quantiles are obtained via simulation. By default, 5000 simulations is used, but this can be changed via the \code{n.sim} argument. Also by default, the simulations uses a classical bootstrap to sample from the standardised residuals. To use an alternative set of standardised innovations, for example the standard normal, use the \code{innov} argument }

\usage{
  \method{predict}{larch}(object, n.ahead=12, newvxreg=NULL, newindex=NULL, 
    n.sim=NULL, innov=NULL, probs=NULL, quantile.type=7, verbose = FALSE, ...)
}

\arguments{
  \item{object}{an object of class 'larch'}
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated (the default is 12)}
  \item{newvxreg}{a \code{matrix} of \code{n.ahead} rows and \code{NCOL(vxreg)} columns with the out-of-sample values of the \code{vxreg} regressors }
  \item{newindex}{\code{NULL} (default) or the date-index for the \code{\link[zoo]{zoo}} object returned by \code{predict.larch}. If \code{NULL}, then the function uses the in-sample \code{index} to generate the out-of-sample index}
  \item{n.sim}{\code{NULL} (default) or an \code{integer}, the number of replications used for the generation of the forecasts. If \code{NULL}, the number of simulations is determined internally (usually 5000) }
  \item{innov}{\code{NULL} (default) or a vector of length \code{n.ahead * n.sim} containing the standardised errors (i.e. mean zero and unit variance) used for the forecast simulations. If \code{NULL}, then a classic bootstrap procedure is used to draw from the standardised in-sample residuals}
  \item{probs}{\code{NULL} (default) or a \code{vector} with the quantile-levels (values strictly between 0 and 1) of the forecast distribution. If \code{NULL}, then no quantiles are returned }
  \item{quantile.type}{an integer between 1 and 9 that selects which algorithm to be used in computing the quantiles, see the argument \code{type} in \code{\link{quantile}} }
  \item{verbose}{logical with default \code{FALSE}. If \code{TRUE}, then additional information (typically the quantiles and/or the simulated series) used in the generation of forecasts is returned. If \code{FALSE}, then only the forecasts are returned}
  \item{\dots}{additional arguments}
}

\value{
a \code{vector} of class \code{\link[zoo]{zoo}} containing the out-of-sample forecasts, or a \code{matrix} of class \code{\link[zoo]{zoo}} containing the out-of-sample forecasts together with additional information (e.g. the prediction-quantiles)
}

\details{No details for the moment.}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\seealso{\code{\link{larch}}}

\examples{
##Simulate some data:
set.seed(123)
e <- rnorm(40)

##estimate log-ARCH(1) model:
mymod <- larch(e, arch=1)

##generate out-of-sample forecasts:
predict(mymod)

##same, but return also selected quantiles:
predict(mymod, probs=c(0.10,0.90))

##same, but using standard normals (instead of bootstrap) in the simulations:
n.sim <- 2000
n.ahead <- 12 #the default on n.ahead
predict(mymod, probs=c(0.10,0.90), n.sim=n.sim, innov=rnorm(n.ahead*n.sim))

##make x-regressors:
x <- matrix(rnorm(40*2), 40, 2)

##estimate log-ARCH(1) model w/covariates:
mymod <- larch(e, arch=1, vxreg=x)

##predict up to 5 steps ahead, setting x's to 0 out-of-sample:
predict(mymod, n.ahead=5, newvxreg=matrix(0,5,NCOL(x)))

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
