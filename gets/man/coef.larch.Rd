\name{coef.larch}
\alias{coef.larch}
\alias{fitted.larch}
\alias{logLik.larch}
\alias{model.matrix.larch}
\alias{nobs.larch}
\alias{plot.larch}
\alias{print.larch}
\alias{residuals.larch}
\alias{summary.larch}
\alias{toLatex.larch}
\alias{vcov.larch}

\title{Methods and extraction functions for 'larch' objects}

\description{Methods and extraction functions for 'larch' objects}

\usage{
  \method{coef}{larch}(object, ...)
  \method{fitted}{larch}(object, ...)
  \method{logLik}{larch}(object, ...)
  \method{model.matrix}{larch}(object, response=FALSE, as.zoo=TRUE, ...)
  \method{nobs}{larch}(object, ...)
  \method{plot}{larch}(x, col=c("red","blue"), lty=c("solid","solid"),
    lwd=c(1,1), ...)
  \method{print}{larch}(x, signif.stars=TRUE, verbose=FALSE, ...)
  \method{residuals}{larch}(object, ...)
  \method{summary}{larch}(object, ...)
  \method{toLatex}{larch}(object, ...)
  \method{vcov}{larch}(object, ...)
}

\arguments{
  \item{object}{an object of class 'larch'}
  \item{x}{an object of class 'larch'}
  \item{response}{logical. If \code{FALSE} (default), the returned object does not contain the response (log(e^2)) used in the estimation}
  \item{as.zoo}{logical. If \code{TRUE} (default), the returned object is of class \code{\link[zoo]{zoo}} }
  \item{col}{a character vector of length two with the colours of actual (default=blue) and fitted (default=red) lines}
  \item{lty}{types of actual (default=solid) and fitted (default=solid) lines}
  \item{lwd}{widths of actual (default=1) and fitted (default=1) lines}
  \item{signif.stars}{logical. If \code{TRUE}, then p-values are additionally encoded visually, see \code{\link{printCoefmat}} }
  \item{verbose}{logical. If \code{TRUE}, and if \code{x} is the results of GETS modelling, then additional information pertaining to the GETS modelling is printed }
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{a vector containing the parameter estimates}
  \item{fitted:}{a \code{\link[zoo]{zoo}} object with fitted values}
  \item{logLik:}{the log-likelihood (normal density)}
  \item{model.matrix:}{the model matrix (see \code{\link{model.matrix}}) of class \code{\link[zoo]{zoo}} (default) of the model}
  \item{nobs:}{the number of observations}
  \item{plot:}{a plot of the fitted values and the residuals}
  \item{print:}{a print of the estimation results and, if \code{verbose=TRUE} and \code{x} is the result of GETS modelling, additional information pertaining to the GETS modelling }
  \item{residuals:}{a \code{\link[zoo]{zoo}} object with the standardised residuals}
  \item{summary:}{a print of the items in the \code{\link{larch}} object}
  \item{toLatex:}{a LaTeX print of the estimation results (equation format)}
  \item{vcov:}{variance-covariance matrix}
}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/} }

\seealso{\code{\link{larch}}, \code{\link[zoo]{zoo}} }

\examples{
##simulate some data:
set.seed(123)
e <- rnorm(40)
x <- matrix(rnorm(40*2), 40, 2)

##estimate a log-ARCH(3)-X model:
mymod <- larch(e, arch=1:3, vxreg=x)

##print results:
print(mymod)

##LaTeX print of the estimation results (equation format):
toLatex(mymod)

##plot the fitted vs. actual values, and the standardised residuals:
plot(mymod)

##extract coefficient estimates (automatically determined):
coef(mymod)

##extract the fitted values:
fitted(mymod)

##extract the standardised residuals:
residuals(mymod)

##extract variance-covariance matrix:
vcov(mymod)

##extract log-likelihood (based on the normal density):
logLik(mymod)

##extract the model matrix of the model:
model.matrix(mymod)

##print the entries of object 'mymod':
summary(mymod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
