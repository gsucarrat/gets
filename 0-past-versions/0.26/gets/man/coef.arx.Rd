\name{coef.arx}
\alias{coef.arx}
\alias{fitted.arx}
\alias{logLik.arx}
\alias{plot.arx}
\alias{print.arx}
\alias{residuals.arx}
\alias{sigma.arx}
\alias{summary.arx}
\alias{vcov.arx}

\title{Extraction functions for 'arx' objects}

\description{Extraction functions for objects of class 'arx'}

\usage{
  \method{coef}{arx}(object, spec=NULL, ...)
  \method{fitted}{arx}(object, spec=NULL, ...)
  \method{logLik}{arx}(object, ...)
  \method{plot}{arx}(x, spec=NULL, col=c("red","blue"),
    lty=c("solid","solid"), lwd=c(1,1), ...)
  \method{print}{arx}(x, signif.stars=TRUE, ...)
  \method{residuals}{arx}(object, std=FALSE, ...)
  \method{sigma}{arx}(object, ...)
  \method{summary}{arx}(object, ...)
  \method{vcov}{arx}(object, spec=NULL, ...)
}

\arguments{
  \item{object}{an object of class 'arx'}
  \item{x}{an object of class 'arx'}
  \item{spec}{\code{NULL}, \code{"mean"}, \code{"variance"} or, in some instances, "both". When \code{NULL} is a valid value, then it is automatically determined whether information pertaining to the mean or variance specification should be returned}
  \item{signif.stars}{\code{logical}. If \code{TRUE}, then p-values are additionally encoded visually, 
  see \code{\link{printCoefmat}} }
  \item{std}{\code{logical}. If \code{FALSE} (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{col}{colours of actual (default=blue) and fitted (default=red) lines}
  \item{lty}{types of actual (default=solid) and fitted (default=solid) lines}
  \item{lwd}{widths of actual (default=1) and fitted (default=1) lines}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{a numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{log-likelihood (normal density)}
  \item{plot:}{a plot of the fitted values and the residuals}
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the residuals}
  \item{sigma:}{the regression standard error ('SE of regression')}
  \item{summary:}{a print of the items in the \code{\link{arx}} object}
  \item{vcov:}{variance-covariance matrix}
}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}\cr
Moritz Schwarz, \url{https://www.inet.ox.ac.uk/people/moritz-schwarz/}\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{arx}}}

\examples{
##simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 40)

##simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*40), 40, 4)

##estimate an 'arx' model: An AR(2) with intercept and four conditioning
##regressors in the mean, and log-ARCH(3) in the variance:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3)

##print results:
print(mymod)

##plot the fitted vs. actual values, and the residuals:
plot(mymod)

##print the entries of object 'mymod':
summary(mymod)

##extract coefficient estimates (automatically determined):
coef(mymod)

##extract mean coefficients only:
coef(mymod, spec="mean")

##extract log-variance coefficients only:
coef(mymod, spec="variance")

##extract all coefficient estimates:
coef(mymod, spec="both")

##extract regression standard error:
sigma(mymod)

##extract log-likelihood:
logLik(mymod)

##extract variance-covariance matrix of mean equation:
vcov(mymod)

##extract variance-covariance matrix of log-variance equation:
vcov(mymod, spec="variance")

##extract and plot the fitted mean values (automatically determined):
mfit <- fitted(mymod)
plot(mfit)

##extract and plot the fitted variance values:
vfit <- fitted(mymod, spec="variance")
plot(vfit)

##extract and plot both the fitted mean and variance values:
vfit <- fitted(mymod, spec="both")
plot(vfit)

##extract and plot the fitted mean values:
vfit <- fitted(mymod, spec="mean")
plot(vfit)

##extract and plot residuals:
epshat <- residuals(mymod)
plot(epshat)

##extract and plot standardised residuals:
zhat <- residuals(mymod, std=TRUE)
plot(zhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
