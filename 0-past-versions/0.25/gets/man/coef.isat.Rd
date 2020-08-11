\name{coef.isat}
\alias{coef.isat}
\alias{fitted.isat}
\alias{logLik.isat}
\alias{plot.isat}
\alias{predict.isat}
\alias{print.isat}
\alias{residuals.isat}
\alias{sigma.isat}
\alias{summary.isat}
\alias{vcov.isat}

\title{Extraction functions for 'isat' objects}

\description{Extraction functions for objects of class 'isat'}

\usage{
  \method{coef}{isat}(object, ...)
  \method{fitted}{isat}(object, ...)
  \method{logLik}{isat}(object, ...)
  \method{plot}{isat}(x, col=c("red","blue"), lty=c("solid","solid"),
    lwd=c(1,1), coef.path=TRUE, ...)
  \method{predict}{isat}(object, n.ahead=12, newmxreg=NULL, newindex=NULL,
    n.sim=2000, probs=NULL, ci.levels=NULL, quantile.type=7,
    return=TRUE, verbose=FALSE, plot=NULL, plot.options=list(), ...)
  \method{print}{isat}(x, signif.stars=TRUE, ...)
  \method{residuals}{isat}(object, std=FALSE, ...)
  \method{sigma}{isat}(object, ...)
  \method{summary}{isat}(object, ...)
  \method{vcov}{isat}(object, ...)
}

\arguments{
  \item{object}{an object of class 'isat'}
  \item{x}{an object of class 'isat'}
  \item{std}{logical. If \code{FALSE} (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated (the default is 12)}
  \item{newmxreg}{a \code{matrix} of \code{n.ahead} rows and \code{NCOL(mxreg)} columns with the out-of-sample values of the \code{mxreg} regressors}
  \item{newindex}{\code{NULL} (default) or the date-index for the \code{zoo} object returned by \code{predict.arx}. If \code{NULL}, then the function uses the in-sample \code{index} to generate the out-of-sample index}
  \item{n.sim}{\code{integer}, the number of replications used for the generation of the forecasts}
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
  \item{coef.path}{logical. Only applicable if there are retained indicators after the application of \code{isat}}
  \item{signif.stars}{\code{logical}. If \code{TRUE}, then p-values are additionally encoded visually, 
  see \code{\link{printCoefmat}} }
  \item{\dots}{additional arguments}
}

\details{
The \code{plot.options} argument is a \code{list} that controls the prediction plot, see 'Details' in \code{\link{predict.arx}}
}

\value{
  \item{coef:}{numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{a numeric, the log-likelihood (normal density)}
  \item{plot:}{plot of the fitted values and the residuals}
  \item{predict:}{a \code{vector} of class \code{\link{zoo}} containing the out-of-sample forecasts, or a \code{matrix} of class \code{\link{zoo}} containing the out-of-sample forecasts together with prediction-quantiles, or - if \code{return=FALSE} - \code{NULL} }
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the residuals}
  \item{sigma:}{the regression standard error ('SE of regression')}
  \item{summary:}{a print of the items in the \code{\link{isat}} object}
  \item{vcov:}{variance-covariance matrix}
}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{ \code{\link{paths}}, \code{\link{terminals}}, \code{\link{coef.gets}}, \code{\link{getsm}}, \code{\link{arx}} }

\examples{
##step indicator saturation:
set.seed(123)
y <- rnorm(30)
isatmod <- isat(y)

##print results:
print(isatmod)

##plot the fitted vs. actual values, and the residuals:
plot(isatmod)

##print the entries of object 'isatmod':
summary(isatmod)

##extract coefficients of the simplified (specific) model:
coef(isatmod)

##extract log-likelihood:
logLik(isatmod)

##extract the coefficient-covariance matrix of simplified
##(specific) model:
vcov(isatmod)

##extract and plot the fitted values:
mfit <- fitted(isatmod)
plot(mfit)

##extract and plot (mean) residuals:
epshat <- residuals(isatmod)
plot(epshat)

##extract and plot standardised residuals:
zhat <- residuals(isatmod, std=TRUE)
plot(zhat)

##generate forecasts of the simplified (specific) model:
predict(isatmod, newmxreg=matrix(1,12,1), plot=TRUE)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
