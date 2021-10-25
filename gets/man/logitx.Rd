\name{logitx}
\alias{logitx}
\alias{dlogitx}

\title{Estimate an autoregressive logit model with covariates}

\description{
Estimate a dynamic Autoregressive (AR) logit model with covariates ('X') by maximising the logit likelihood.
}

\usage{
logitx(y, intercept = TRUE, ar = NULL, ewma = NULL, xreg = NULL, 
    vcov.type = c("ordinary", "robust"), lag.length = NULL, 
    initial.values = NULL, lower = -Inf, upper = Inf, control = list(), 
    eps.tol = .Machine$double.eps, solve.tol = .Machine$double.eps, 
    plot = NULL)

dlogitx(y, ...)
}

\arguments{
  \item{y}{a binary numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{intercept}{logical. \code{TRUE}, the default, includes an intercept in the logit specification, whereas \code{FALSE} does not}
  \item{ar}{either \code{NULL} (default) or an integer vector, say, \code{c(2,4)} or \code{1:4}. The AR-lags to include in the logit specification. If \code{NULL}, then no lags are included}
  \item{ewma}{either \code{NULL} (default) or a \code{\link{list}} with arguments sent to the \code{\link{eqwma}} function. In the latter case a lagged moving average of \code{y} is included as a regressor}
  \item{xreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link{zoo}} object, of covariates. Note that, if both \code{y} and \code{xreg} are \code{zoo} objects, then their samples are chosen to match}
  \item{vcov.type}{character vector of length 1, either "ordinary" (default) or "robust". Partial matching is allowed. If "ordinary", then the ordinary variance-covariance matrix is used for inference. If "robust", then a robust coefficient-covariance of the Newey and West (1987) type is used }
  \item{lag.length}{\code{NULL} or an integer that determines the lag-length used in the robust coefficient covariance. If \code{lag.length} is an integer, then it is ignored unless \code{method = 3} }
  \item{initial.values}{\code{NULL} or a numeric vector with the initial parameter values passed on to the optimisation routine, \code{\link{nlminb}}. If \code{NULL}, the default, then the values are chosen automatically }
  \item{lower}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}} }
  \item{upper}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}} }
  \item{control}{a \code{list} passed on to the control argument of \code{\link{nlminb}} }
  \item{eps.tol}{numeric, a small value that ensures the fitted zero-probabilities are not too small when the log-transformation is applied when computing the log-likelihood }
  \item{solve.tol}{numeric value passed on to the \code{tol} argument of \code{\link{solve}}, which is called whenever the coefficient-coariance matrix is computed. The value controls the toleranse for detecting linear dependence between columns when inverting a matrix }
  \item{plot}{\code{NULL} or logical. If \code{TRUE}, then a plot is produced. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
  \item{...}{arguments passed on to \code{logitx}}
}

\details{
The function estimates a dynamic Autoregressive (AR) logit model with (optionally) covariates ('X') by maximising the logit likelihood. The estimated model is that of Kauppi and Saikkonen (2008). However, there (in contrast to here) estimation is by maximisation of the probit likelihood. 
}

\value{A list of class 'logitx'.}

\references{
Heikki Kauppi and Penti Saikkonen (2008): 'Predicting U.S. Recessions with Dynamic Binary Response Models'. The Review of Economic Statistics 90, pp. 777-791

Whitney K. Newey and Kenned D. West (1987): 'A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
Methods: \code{\link{coef.logitx}}, \code{\link{fitted.logitx}}, \code{\link{gets.logitx}}, \cr
\code{\link{logLik.logitx}}, \code{\link{plot.logitx}}, \code{\link{print.logitx}}, \code{\link{summary.logitx}}, \code{\link{toLatex.logitx}} and \code{\link{vcov.logitx}}\cr

Related functions: \code{\link{logitxSim}}, \code{\link{logit}}, \code{\link{nlminb}}
}

\examples{

##simulate from ar(1):
set.seed(123) #for reproducibility
y <- logitxSim(100, ar=0.3)

##estimate ar(1) and store result:
mymod <- logitx(y, ar=1)

##estimate ar(4) and store result:
mymod <- logitx(y, ar=1:4)

##create some more data, estimate new model:
x <- matrix(rnorm(5*100), 100, 5)
mymod <- logitx(y, ar=1:4, xreg=x)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
