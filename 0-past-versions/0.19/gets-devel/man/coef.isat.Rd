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
  \method{fitted}{isat}(object, spec=NULL, ...)
  \method{logLik}{isat}(object, ...)
  \method{plot}{isat}(x, col=c("red","blue"), lty=c("solid","solid"),
    lwd=c(1,1), coef.path=TRUE, ...)
  \method{predict}{isat}(object, n.ahead=12, newmxreg=NULL, newindex=NULL,
    return=TRUE, plot=NULL, plot.options=list(), ...)
  \method{print}{isat}(x, ...)
  \method{residuals}{isat}(object, std=FALSE, ...)
  \method{sigma}{isat}(object, ...)
  \method{summary}{isat}(object, ...)
  \method{vcov}{isat}(object, ...)
}

\arguments{
  \item{object}{an object of class 'isat'}
  \item{x}{an object of class 'isat'}
  \item{spec}{'mean' or 'variance'}
  \item{std}{logical. If FALSE (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{n.ahead}{generate forecasts up to n steps ahead (the default is 12)}
  \item{newmxreg}{a matrix (\code{n.ahead} rows and NCOL(mxregs) columns) with the out-of-sample values of the \code{mxreg} regressors}
  \item{newindex}{date-index for the \code{zoo} object returned by \code{predict.arx}}
  \item{return}{logical. If \code{TRUE} (default), then the out-of-sample forecasts are returned}
  \item{plot}{NULL or logical. If \code{TRUE}, then the out-of-sample forecasts are plotted. If NULL (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
  \item{plot.options}{a \code{list} of options related to the plotting of forecasts, see 'Details'}
  \item{col}{colours of fitted (default=red) and actual (default=blue) lines}
  \item{lty}{types of fitted (default=solid) and actual (default=solid) lines}
  \item{lwd}{widths of fitted (default=1) and actual (default=1) lines}
  \item{coef.path}{logical. Only applicable if there are retained indicators after the application of \code{isat}}
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
  \item{coef:}{numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{a numeric, the log-likelihood (normal density)}
  \item{plot:}{plot of the fitted values and the residuals}
  \item{predict}{a vector containing the out-of-sample forecasts}
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

##extract variance-covariance matrix of simplified
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

##generate forecasts with the simplified
##(specific) model:
predict(isatmod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
