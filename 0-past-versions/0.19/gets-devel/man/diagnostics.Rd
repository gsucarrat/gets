\name{diagnostics}
\alias{diagnostics}

\title{Diagnostics tests}

\description{Auxiliary function (i.e. not intended for the average user) called by the \code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}} and \code{\link{getsFun}} functions. The \code{diagnostics} function undertakes tests for autocorrelation, ARCH and non-normality in a residual series, and user-defined diagnostics provided via the \code{user.fun} argument (see details). The autocorrelation and ARCH tests are conducted as Ljung and Box (1979) tests for autocorrelation in the residuals and squared residuals, respectively, whereas the test for non-normality is that of Jarque and Bera (1980).}

\usage{
diagnostics(x, ar.LjungB=c(1, 0.025), arch.LjungB=c(1, 0.025),
  normality.JarqueB=NULL, verbose=TRUE, user.fun=NULL, ...)
}

\arguments{
  \item{x}{a \code{\link{list}}, typically the estimation result of \code{\link{ols}}}
  \item{ar.LjungB}{a two element vector or \code{NULL}. In the former case, the first element contains the AR-order, the second element the significance level. If \code{NULL}, then a test for autocorrelation is not conducted}
  \item{arch.LjungB}{a two element vector or \code{NULL}. In the former case, the first element contains the ARCH-order, the second element the significance level. If \code{NULL}, then a test for ARCH is not conducted}
  \item{normality.JarqueB}{\code{NULL} or a value between 0 and 1. In the latter case, a test for non-normality is conducted using a significance level equal to \code{normality.JarqueB}. If \code{NULL}, then no test for non-normality is conducted}
  \item{verbose}{logical. If \code{TRUE}, then a \code{\link{data.frame}} with the results of the diagnostics is returned. If \code{FALSE}, then the return-value is a \code{logical} that indicates whether the model passes the diagnostics (\code{TRUE} if it does, otherwise \code{FALSE}) }
  \item{user.fun}{\code{NULL} or a \code{\link{list}} with two entries, \code{name} and \code{pval}. The first item (\code{name}) should contain the name of the user-defined function, and must be of class \code{character}. The the second item should contain the chosen significance level or levels, i.e. either a scalar or a vector of length equal to the number of p-values returned by the user-defined diagnostics function, see details }
  \item{...}{further arguments (ignored) to accommodate deleted arguments from past versions of the functions}
}

\details{
  The argument \code{user.fun} enables the user to specify additional diagnostics. The function must be defined in the global environment (i.e. \code{.GlobalEnv}), and the user-defined function should refer to the named items of the estimation result \code{x} (see examples).
  The value returned by the user-defined function should be a matrix of dimension m x 3, where m is the number of diagnostic tests performed inside the user-defined function. (If only a single test is performed, then the returned value can be a vector of length 3.) The three columns of the m x 3 matrix should contain, in the following order: 1) the value(s) of the test-statistic(s), 2) the degrees of freedom(s) (or \code{NA} if there are none) of the tests, and 3) the p-value(s) of the test(s). When checking whether the model passes the diagnostics or not, the p-value(s) is(are) checked against the value(s) in the second entry of \code{user.fun}.
}

\value{
  \item{If \code{verbose=TRUE}:}{a \code{\link{data.frame}} that contains the diagnostics results}
  \item{If \code{verbose=FALSE}:}{a \code{logical} indicating whether the residuals and/or model passes the diagnostics (\code{TRUE} if it does, \code{FALSE} otherwise)}
}

\references{
C. Jarque and A. Bera (1980): 'Efficient Tests for Normality,
  Homoscedasticity and Serial Independence'. Economics Letters 6,
  pp. 255-259

G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time
  Series Models'. Biometrika 66, pp. 265-270
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}, \code{\link{getsFun}} }

\examples{
##return a data-frame with autocorrelation and ARCH diagnostics (default),
##and check whether they pass (the default p-value is 0.025):
set.seed(123)
vY <- rnorm(20)
mX <- matrix(rnorm(3*20), 20, 3)
est <- ols(vY,mX)
diagnostics(est)
diagnostics(est, verbose=FALSE)

##add the Jarque-Bera normality test to the diagnostics (w/p-value=0.05):
diagnostics(est, normality.JarqueB=0.05)
diagnostics(est, normality.JarqueB=0.05, verbose=FALSE)

##user-defined Shapiro-Wilks test for non-normality of the residuals:
SWtest <- function(residuals, ...){
  tmp <- shapiro.test(residuals) #do test on est$residuals
  return( c(tmp$statistic, NA, tmp$p.value) )
}
diagnostics(est, user.fun=list(name="SWtest", pval=0.05))
diagnostics(est, user.fun=list(name="SWtest", pval=0.05), verbose=FALSE)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
