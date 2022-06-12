\name{printtex}
\alias{printtex}
\alias{toLatex.arx}
\alias{toLatex.gets}

\title{Generate LaTeX code of an estimation result}

\description{
Convenience functions that generates LaTeX-code of an estimation result in equation-form. \code{printtex} can, in principle, be applied to any object for which \code{coef}, \code{vcov} and \code{logLik} methods exist. Note: The generated LaTeX-code contains an \code{eqnarray} environment, which requires that the \code{amsmath} package is loaded in the preamble of the LaTeX document.
}

\usage{
  printtex(x, fitted.name=NULL, xreg.names=NULL, digits=4,
      intercept=TRUE, gof=TRUE, diagnostics=TRUE, nonumber=FALSE,
      nobs="T", index="t", dec=NULL, print.info=TRUE)
  \method{toLatex}{arx}(object, ...)
  \method{toLatex}{gets}(object, ...)
}
\arguments{
  \item{x}{an estimation result, e.g. \code{\link{arx}}, \code{gets} or \code{isat} object }
  \item{object}{an estimation result of class \code{\link{arx}} or \code{gets} }
  \item{fitted.name}{\code{NULL} or a user-specified name of left-hand side variable}
  \item{xreg.names}{\code{NULL} or a user-specified character vector with the names of regressors}
  \item{digits}{integer, the number of digits to be printed}
  \item{intercept}{logical or numeric. The argument determines whether one of the regressors is an intercept or not, or its location. If \code{TRUE}, then the intercept is assumed to be located at \code{coef(x)[1]}, and hence the regressor-name of location 1 is excluded from the print. If \code{FALSE}, then it is assumed that there is no intercept among the regressors. If numeric, then it is assumed that the regressors contain an intercept at the location equal to the numeric value}
  \item{gof}{logical, whether to include goodness-of-fit in the print}
  \item{diagnostics}{logical, whether to include diagnostics in the print}
  \item{nonumber}{logical, whether to remove or not (default) the equation-numbering}
  \item{nobs}{character, the notation to use to denote the number of observations}
  \item{index}{\code{NULL} or a \code{\link{character}}, only relevant if \code{fitted.name} is not \code{NULL}, and if the object in question is of class \code{\link{arx}}, \code{\link{gets}} or \code{\link{isat}} }
  \item{dec}{\code{NULL} or a \code{\link{character}} (for example \code{","}). In the latter case, an attempt is made to replace the dot separator \code{.} with the character in \code{dec} }
  \item{print.info}{\code{logical}, whether to print the info at the start or not}
  \item{...}{arguments passed on to \code{printtex} }
}

\details{\code{toLatex.arx} and \code{toLatex.gets} are simply wrappers to \code{printtex} }

\value{LaTeX code of an estimation result}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{arx}}, \code{\link{logitx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}}

\examples{
##simulate random variates, estimate model:
y <- rnorm(30)
mX <- matrix(rnorm(30*2), 30, 2)
mymod <- arx(y, ar=1:3, mxreg=mX)

##print latex code of estimation result:
printtex(mymod)

##add intercept, at the end, to regressor matrix:
mX <- cbind(mX,1)
colnames(mX) <- c("xreg1", "xreg2", "intercept")
mymod <- arx(y, mc=FALSE, mxreg=mX)

##set intercept location to 3:
printtex(mymod, intercept=3)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
