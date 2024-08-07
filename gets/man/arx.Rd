\name{arx}
\alias{arx}
\title{Estimate an AR-X model with log-ARCH-X errors}

\description{

Estimation by OLS, two-step OLS if a variance specification is specified: In the first the mean specification (AR-X) is estimated, whereas in the second step the log-variance specification (log-ARCH-X) is estimated.

The AR-X mean specification can contain an intercept, AR-terms, lagged moving averages of the regressand and other conditioning covariates ('X'). The log-variance specification can contain log-ARCH terms, asymmetry or 'leverage' terms, log(EqWMA) where EqWMA is a lagged equally weighted moving average of past squared residuals (a volatility proxy) and other conditioning covariates ('X').

}

\usage{
arx(y, mc=TRUE, ar=NULL, ewma=NULL, mxreg=NULL, vc=FALSE,
  arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL, zero.adj=NULL,
  vc.adj=TRUE, vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
  user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, singular.ok=TRUE,
  plot=NULL)
}

\arguments{
  \item{y}{\code{numeric} vector, time-series or \code{\link[zoo]{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link[zoo]{na.trim}} command}
  \item{mc}{\code{logical}. \code{TRUE} (default) includes an intercept in the mean specification, whereas \code{FALSE} does not}
  \item{ar}{either \code{NULL} (default) or an integer vector, say, \code{c(2,4)} or \code{1:4}. The AR-lags to include in the mean specification. If \code{NULL}, then no lags are included}
  \item{ewma}{either \code{NULL} (default) or a \code{\link{list}} with arguments sent to the \code{\link{eqwma}} function. In the latter case a lagged moving average of \code{y} is included as a regressor}
  \item{mxreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link[zoo]{zoo}} object, of conditioning variables. Note that, if both \code{y} and \code{mxreg} are \code{\link[zoo]{zoo}} objects, then their samples are chosen to match}
  \item{vc}{\code{logical}. \code{TRUE} includes an intercept in the log-variance specification, whereas \code{FALSE} (default) does not. If the log-variance specification contains any other item but the log-variance intercept, then vc is set to \code{TRUE} }
  \item{arch}{either \code{NULL} (default) or an integer vector, say, \code{c(1,3)} or \code{2:5}. The log-ARCH lags to include in the log-variance specification}
  \item{asym}{either \code{NULL} (default) or an integer vector, say, \code{c(1)} or \code{1:3}. The asymmetry (i.e. 'leverage') terms to include in the log-variance specification}
  \item{log.ewma}{either \code{NULL} (default) or a vector of the lengths of the volatility proxies, see \code{\link{leqwma}} }
  \item{vxreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link[zoo]{zoo}} object, of conditioning variables. If both \code{y} and \code{mxreg} are \code{\link[zoo]{zoo}} objects, then their samples are chosen to match. }
  \item{zero.adj}{\code{NULL} (default) or a strictly positive \code{numeric} scalar. If \code{NULL}, the zeros in the squared residuals are replaced by the 10 percent quantile of the non-zero squared residuals. If \code{zero.adj} is a strictly positive \code{numeric} scalar, then this value is used to replace the zeros of the squared residuals. }
  \item{vc.adj}{\code{logical}. If \code{TRUE} (default), then the log-variance intercept is adjusted by the estimate of E[ln(z^2)], where z is the standardised error. This adjustment is needed for the conditional scale to be equal to the conditional standard deviation. If \code{FALSE}, then the log-variance intercept is not adjusted}
  \item{vcov.type}{\code{character} vector, "ordinary" (default), "white" or "newey-west". If "ordinary", then the ordinary variance-covariance matrix is used for inference. If "white", then the White (1980) heteroscedasticity-robust matrix is used. If "newey-west", then the Newey and West (1987) heteroscedasticity and autocorrelation-robust matrix is used}
  \item{qstat.options}{\code{NULL} (default) or an integer vector of length two, say, \code{c(1,1)}. The first value sets the lag-order of the AR diagnostic test, whereas the second value sets the lag-order of the ARCH diagnostic test. If \code{NULL}, then the two values of the vector are set automatically}
  \item{normality.JarqueB}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, then the results of the Jarque and Bera (1980) test for non-normality in the residuals are included in the estimation results. }
  \item{user.estimator}{\code{NULL} (default) or a \code{\link{list}} with one entry, \code{name}, containing the name of the user-defined estimator. Additional items, if any, are passed on as arguments to the estimator in question }
  \item{user.diagnostics}{\code{NULL} (default) or a \code{\link{list}} with two entries, \code{name} and \code{pval}, see the \code{user.fun} argument in \code{\link{diagnostics}} }
  \item{tol}{\code{numeric} value (\code{default = 1e-07}). The tolerance for detecting linear dependencies in the columns of the regressors (see \code{\link{qr}} function). Only used if \code{LAPACK} is \code{FALSE} (default) and \code{user.estimator} is \code{NULL}.}
  \item{LAPACK}{\code{logical}. If \code{TRUE}, then use LAPACK. If \code{FALSE} (default), then use LINPACK (see \code{\link{qr}} function). Only used if \code{user.estimator} is \code{NULL}.}
  \item{singular.ok}{\code{logical}. If \code{TRUE} (default), the regressors are checked for singularity, and the ones causing it are automatically removed. }
  \item{plot}{\code{NULL} or \code{logical}. If \code{TRUE}, the fitted values and the residuals are plotted. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
}

\details{

For an overview of the AR-X model with log-ARCH-X errors, see Pretis, Reade and Sucarrat (2018): \doi{10.18637/jss.v086.i03}.\cr

The arguments \code{user.estimator} and \code{user.diagnostics} enables the specification of user-defined estimators and user-defined diagnostics. To this end, the principles of the same arguments in \code{\link{getsFun}} are followed, see its documentation under "Details", and Sucarrat (2020): \url{https://journal.r-project.org/archive/2021/RJ-2021-024/}.

}

\value{A list of class 'arx'}

\references{

C. Jarque and A. Bera (1980): 'Efficient Tests for Normality, Homoscedasticity and Serial Independence'. Economics Letters 6, pp. 255-259. \doi{10.1016/0165-1765(80)90024-5}

Felix Pretis, James Reade and Genaro Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. \doi{10.18637/jss.v086.i03}

Genaro Sucarrat (2020): 'User-Specified General-to-Specific and Indicator Saturation Methods'. The R Journal 12:2, pages 388-401. \url{https://journal.r-project.org/archive/2021/RJ-2021-024/}

Halbert White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838.

Whitney K. Newey and Kenned D. West (1987): 'A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708.

}

\author{

\tabular{lll}{
	Jonas Kurle: \tab \tab \url{https://www.jonaskurle.com/}\cr
	Moritz Schwarz: \tab \tab  \url{https://www.inet.ox.ac.uk/people/moritz-schwarz}\cr
	Genaro Sucarrat: \tab \tab \url{https://www.sucarrat.net/}\cr
}

}

\seealso{
Extraction functions (mostly S3 methods): \code{\link{coef.arx}}, \code{\link{ES}}, \code{\link{fitted.arx}}, \code{\link{plot.arx}}, \cr
\code{\link{print.arx}}, \code{\link{recursive}}, \code{\link{residuals.arx}}, \code{\link{sigma.arx}}, \code{\link{rsquared}},\cr
\code{\link{summary.arx}}, \code{\link{VaR}} and \code{\link{vcov.arx}}\cr

Related functions: \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}
}

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 70)

##estimate an AR(2) with intercept:
arx(y, mc=TRUE, ar=1:2)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*70), 70, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean:
arx(y, ar=1:2, mxreg=xregs)

##estimate a log-variance specification with a log-ARCH(4)
##structure:
arx(y, mc=FALSE, arch=1:4)

##estimate a log-variance specification with a log-ARCH(4)
##structure and an asymmetry/leverage term:
arx(y, mc=FALSE, arch=1:4, asym=1)

##estimate a log-variance specification with a log-ARCH(4)
##structure, an asymmetry or leverage term, a 10-period log(EWMA) as
##volatility proxy, and the log of the squareds of the conditioning
##regressors in the log-variance specification:
arx(y, mc=FALSE,
  arch=1:4, asym=1, log.ewma=list(length=10), vxreg=log(xregs^2))

##estimate an AR(2) with intercept and four conditioning regressors
##in the mean, and a log-variance specification with a log-ARCH(4)
##structure, an asymmetry or leverage term, a 10-period log(EWMA) as
##volatility proxy, and the log of the squareds of the conditioning
##regressors in the log-variance specification:
arx(y, ar=1:2, mxreg=xregs,
  arch=1:4, asym=1, log.ewma=list(length=10), vxreg=log(xregs^2))
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
