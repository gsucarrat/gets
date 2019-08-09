\name{coef.gets}
\alias{coef.gets}
\alias{fitted.gets}
\alias{logLik.gets}
\alias{plot.gets}
\alias{predict.gets}
\alias{print.gets}
\alias{residuals.gets}
\alias{sigma.gets}
\alias{summary.gets}
\alias{vcov.gets}

\title{Extraction functions for 'gets' objects}

\description{Extraction functions for objects of class 'gets'}

\usage{
  \method{coef}{gets}(object, spec=NULL, ...)
  \method{fitted}{gets}(object, spec=NULL, ...)
  \method{logLik}{gets}(object, ...)
  \method{plot}{gets}(x, spec=NULL, col=c("red","blue"),
    lty=c("solid","solid"), lwd=c(1,1), ...)
  \method{predict}{gets}(object, spec=NULL, n.ahead=12, newmxreg=NULL,
    newvxreg=NULL, newindex=NULL, n.sim=5000, innov=NULL, probs=NULL,
    ci.levels=NULL, quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
    plot.options=list(), ...)
  \method{print}{gets}(x, ...)
  \method{residuals}{gets}(object, std=NULL, ...)
  \method{sigma}{gets}(object, ...)
  \method{summary}{gets}(object, ...)
  \method{vcov}{gets}(object, spec=NULL, ...)
}

\arguments{
  \item{object}{an object of class 'gets'}
  \item{x}{an object of class 'gets'}
  \item{spec}{NULL, "mean", "variance" or, in some instances, "both". When \code{NULL} is a valid value, then it is automatically determined whether information pertaining to the mean or variance specification should be returned}
  \item{std}{logical. If FALSE (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated (the default is 12)}
  \item{newmxreg}{a \code{matrix} of \code{n.ahead} rows and \code{NCOL(mxreg)} columns with the out-of-sample values of the \code{mxreg} regressors}
  \item{newvxreg}{a \code{matrix} of \code{n.ahead} rows and \code{NCOL(vxreg)} columns with the out-of-sample values of the \code{vxreg} regressors}
  \item{newindex}{\code{NULL} (default) or the date-index for the \code{zoo} object returned by \code{predict.arx}. If \code{NULL}, then the function uses the in-sample \code{index} to generate the out-of-sample index}
  \item{n.sim}{\code{integer}, the number of replications used for the generation of the forecasts}
  \item{innov}{\code{NULL} (default) or a vector of length \code{n.ahead * n.sim} containing the standardised errors (that is, zero mean and unit variance) used for the forecast simulations. If \code{NULL}, then a classica bootstrap procedure is used to draw from the standardised in-sample residuals}
  \item{probs}{\code{NULL} (default) or a \code{vector} with the quantile-levels (values strictly between 0 and 1) of the forecast distribution. If \code{NULL}, then no quantiles are returned unless \code{ci.levels} is non-\code{NULL}}
  \item{ci.levels}{\code{NULL} (default) or a \code{vector} with the confidence levels (expressed as values strictly between 0 and 1) of the forecast distribution. The upper and lower values of the confidence interval(s) are returned as quantiles}
  \item{quantile.type}{an integer between 1 and 9 that selects which algorithm to be used in computing the quantiles, see the argument \code{type} in \code{\link{quantile}} }
  \item{return}{logical. If \code{TRUE} (default), then the out-of-sample predictions are returned. The value \code{FALSE}, which does not return the predictions, may be of interest if only a prediction plot is of interest}
  \item{verbose}{logical with default \code{FALSE}. If \code{TRUE}, then additional information (typically the quantiles and/or the simulated series) used in the generation of forecasts is returned. If \code{FALSE}, then only the forecasts are returned}
  \item{plot}{\code{NULL} (default) or logical. If \code{NULL}, then the value set by \code{options$plot} (see \code{\link{options}}) determines whether a plot is produced or not. If \code{TRUE}, then the out-of-sample forecasts are plotted.}
  \item{plot.options}{a \code{list} of options related to the plotting of forecasts, see 'Details'}
  \item{col}{colours of fitted (default=red) and actual (default=blue) lines}
  \item{lty}{types of fitted (default=solid) and actual (default=solid) lines}
  \item{lwd}{widths of fitted (default=1) and actual (default=1) lines}
  \item{\dots}{additional arguments}
}

\details{
The \code{plot.options} argument is a \code{list} that controls the prediction plot, see 'Details' in \code{\link{predict.arx}}
}

\value{
  \item{coef:}{a numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{a numeric, the log-likelihood (normal density)}
  \item{plot:}{a plot of the fitted values and the residuals}
  \item{predict:}{a \code{vector} of class \code{\link{zoo}} containing the out-of-sample forecasts, or a \code{matrix} of class \code{\link{zoo}} containing the out-of-sample forecasts together with prediction-quantiles, or - if \code{return=FALSE} - \code{NULL} }
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the residuals}
  \item{sigma:}{the regression standard error ('SE of regression')}
  \item{summary:}{a print of the items in the \code{\link{gets}} object}
  \item{vcov:}{a variance-covariance matrix}
}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{ \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}} }

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and a log-ARCH(3) in the variance:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3)

##General-to-Specific (GETS) model selection of the mean:
meanmod <- getsm(mymod)

##General-to-Specific (GETS) model selection of the variance:
varmod <- getsv(mymod)

##print results:
print(meanmod)
print(varmod)

##plot the fitted vs. actual values, and the residuals:
plot(meanmod)
plot(varmod)

##generate and plot predictions of the mean:
predict(meanmod, plot=TRUE)

##print the entries of object 'gets':
summary(meanmod)
summary(varmod)

##extract coefficients of the simplified (specific) model:
coef(meanmod) #mean spec
coef(varmod) #variance spec

##extract log-likelihood:
logLik(mymod)

##extract coefficient-covariance matrix of simplified
##(specific) model:
vcov(meanmod) #mean spec
vcov(varmod) #variance spec

##extract and plot the fitted values:
mfit <- fitted(meanmod) #mean fit
plot(mfit)
vfit <- fitted(varmod) #variance fit
plot(vfit)

##extract and plot residuals:
epshat <- residuals(meanmod)
plot(epshat)

##extract and plot standardised residuals:
zhat <- residuals(varmod)
plot(zhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
