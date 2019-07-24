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
    newvxreg=NULL, newindex=NULL, n.sim=1000, innov=NULL, return=TRUE,
    plot=NULL, plot.options=list(), ...)
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
  \item{n.ahead}{generate forecasts up to n steps ahead (the default is 12)}
  \item{newmxreg}{a matrix (\code{n.ahead} rows and NCOL(mxregs) columns) with the out-of-sample values of the \code{mxreg} regressors}
  \item{newvxreg}{a matrix (\code{n.ahead} rows and NCOL(vxregs) columns) with the out-of-sample values of the \code{vxreg} regressors}
  \item{newindex}{date-index for the \code{zoo} object returned by \code{predict.arx}}
  \item{n.sim}{integer, the number of bootstrap replications for the generation of the variance forecasts}
  \item{innov}{\code{NULL} (default) or a vector of length \code{n.ahead * n.sim} containing the standardised errors (i.e. zero mean, unit variance) to bootstrap from}
  \item{return}{logical. If \code{TRUE} (default), then the out-of-sample forecasts are returned}
  \item{plot}{NULL or logical. If \code{TRUE} (default), then the out-of-sample forecasts are plotted. If NULL (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
  \item{plot.options}{a \code{list} of options related to the plotting of forecasts, see 'Details'}
  \item{col}{colours of fitted (default=red) and actual (default=blue) lines}
  \item{lty}{types of fitted (default=solid) and actual (default=solid) lines}
  \item{lwd}{widths of fitted (default=1) and actual (default=1) lines}
  \item{\dots}{additional arguments}
}

\details{
The \code{plot.options} argument is a list that can contain any of the following arguments:\cr

\code{keep}: integer greater than zero (default is 12) that controls the number of in-sample actual values to plot\cr
\code{fitted}: If \code{TRUE}, then the fitted values as well as actual values are plotted in-sample. The default is \code{FALSE}\cr
\code{legend.loc}: character string (the default is "topleft"). Allows the location of the plot legend to be altered\cr
\code{newmactual}: numeric vector or \code{NULL} (default). Enables the plotting of actual values out-of-sample in the mean in addition to the forecasts\cr
\code{newvactual}: numeric vector or \code{NULL} (default). Enables the plotting of actual values out-of-sample in the variance addition to the forecasts\cr

}

\value{
  \item{coef:}{a numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{a numeric, the log-likelihood (normal density)}
  \item{plot:}{a plot of the fitted values and the residuals}
  \item{predict}{a vector containing the out-of-sample forecasts}
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

##print the entries of object 'gets':
summary(meanmod)
summary(varmod)

##extract coefficients of the simplified (specific) model:
coef(meanmod) #mean spec
coef(varmod) #variance spec

##extract log-likelihood:
logLik(mymod)

##extract variance-covariance matrix of simplified
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
