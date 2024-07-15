\name{larch}
\alias{larch}

\title{Estimate a heterogeneous log-ARCH-X model}

\description{
The function \code{larch()} estimates a heterogeneous log-ARCH-X model, which is a generalisation of the dynamic log-variance model in Pretis, Reade and Sucarrat (2018). Internally, estimation is undertaken by a call to \code{\link{larchEstfun}}. The log-variance specification can contain log-ARCH terms, log-HARCH terms, asymmetry terms ('leverage'), the log of volatility proxies made up of past returns and other covariates ('X'), for example Realised Volatility (RV), volume or the range.
}

\usage{
larch(e, vc=TRUE, arch = NULL, harch = NULL, asym = NULL, asymind = NULL,
  log.ewma = NULL, vxreg = NULL, zero.adj = NULL, 
  vcov.type = c("robust", "hac"), qstat.options = NULL,
  normality.JarqueB = FALSE, tol = 1e-07, singular.ok = TRUE,  plot = NULL)
}

\arguments{
  \item{e}{\code{numeric} vector, time-series or \code{\link[zoo]{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link[zoo]{na.trim}} command }
  \item{vc}{\code{logical}. \code{TRUE} includes an intercept in the log-variance specification. Currently, \code{vc} cannot be set to any other value than \code{TRUE} }
  \item{arch}{either \code{NULL} (default) or an integer vector, say, \code{c(1,3)} or \code{2:5}. The log-ARCH lags to include in the log-variance specification }
  \item{harch}{either \code{NULL} (default) or an integer vector, say, \code{c(5,10)}. The (log of) heterogeneous ARCH terms (Muller et al. 1997) to include }
  \item{asym}{either \code{NULL} (default) or an integer vector, say, \code{c(1)} or \code{1:3}. The asymmetry (i.e. 'leverage') terms to include in the log-variance specification}
  \item{asymind}{either \code{NULL} (default or an integer vector. The indicator asymmetry terms to include }
  \item{log.ewma}{either \code{NULL} (default) or a vector of the lengths of the volatility proxies, see \code{\link{leqwma}}. The terms serve as (log of) volatility proxies similar to RVs in the HAR-model of Corsi (2009). Here, the \code{log.ewma} terms are made up of past e's }
  \item{vxreg}{either \code{NULL} (default) or a numeric vector or matrix, say, a \code{\link[zoo]{zoo}} object. If both \code{e} and \code{vxreg} are \code{\link[zoo]{zoo}} objects, then their samples are chosen to match }
  \item{zero.adj}{\code{NULL} (default) or a strictly positive \code{numeric} scalar. If \code{NULL}, the zeros in the squared residuals are replaced by the 10 percent quantile of the non-zero squared residuals. If \code{zero.adj} is a strictly positive \code{numeric} scalar, then this value is used to replace the zeros of the squared e's }
  \item{vcov.type}{\code{character}. "robust" (default) or "hac" (partial matching is allowed). If "robust", the robust variance-covariance matrix of the White (1980) type is used. If "hac", the Newey and West (1987) heteroscedasticity and autocorrelation-robust matrix is used}
  \item{qstat.options}{\code{NULL} (default) or an integer vector of length two, say, \code{c(1,1)}. The first value sets the lag-order of the AR diagnostic test of the standardised residuals, whereas the second value sets the lag-order of the ARCH diagnostic test of the standardised residuals. If \code{NULL}, then the two values of the vector are set automatically }
  \item{normality.JarqueB}{\code{FALSE} (default) or \code{TRUE}. If \code{TRUE}, then the results of the Jarque and Bera (1980) test for non-normality in the residuals are included in the estimation results }
  \item{tol}{\code{numeric} value. The tolerance (the default is \code{1e-07}) for detecting linear dependencies in the columns of the regressors (see \code{\link{ols}} and \code{\link{qr}}). Only used if \code{LAPACK} is \code{FALSE} (default) }
  \item{singular.ok}{\code{logical}. If \code{TRUE} (default), the regressors are checked for singularity, and the ones causing it are automatically removed. If \code{FALSE}, then the function returns an error }
  \item{plot}{\code{NULL} (default) or \code{logical}. If \code{TRUE}, the fitted values and the residuals are plotted. If \code{NULL}, then the value set by \code{\link{options}} determines whether a plot is produced or not }
}

\details{No details for the moment}

\value{A list of class 'larch'}

\references{

G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time Series Models'. Biometrika 66, pp. 265-270

F. Corsi (2009): 'A Simple Approximate Long-Memory Model of Realized Volatility', Journal of Financial Econometrics 7, pp. 174-196

C. Jarque and A. Bera (1980): 'Efficient Tests for Normality, Homoscedasticity and Serial Independence'. Economics Letters 6, pp. 255-259. \doi{10.1016/0165-1765(80)90024-5}

U. Muller, M. Dacorogna, R. Dave, R. Olsen, O. Pictet and J. von Weizsacker (1997): 'Volatilities of different time resolutions - analyzing the dynamics of market components'. Journal of Empirical Finance 4, pp. 213-239

F. Pretis, J. Reade and G. Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. \doi{10.18637/jss.v086.i03}

H. White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838.

W.K. Newey and K.D. West (1987): 'A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708.

}

\author{Genaro Sucarrat: \url{https://www.sucarrat.net/}}

\seealso{
Methods and extraction functions (mostly S3 methods): \code{\link{coef.larch}}, \code{\link{ES}}, \code{\link{fitted.larch}}, \code{\link{gets.larch}}, \cr
\code{\link{logLik.larch}}, \code{\link{nobs.larch}}, \code{\link{plot.larch}}, \code{\link{predict.larch}}, \code{\link{print.larch}}, \cr
\code{\link{residuals.larch}}, \code{\link{summary.larch}}, \code{\link{VaR}}, \code{\link{toLatex.larch}} and \code{\link{vcov.arx}}\cr

\code{\link{regressorsVariance}}
}

\examples{
##Simulate some data:
set.seed(123)
e <- rnorm(40)
x <- matrix(rnorm(40*2), 40, 2)

##estimate a log-variance specification with a log-ARCH(4)
##structure:
larch(e, arch=1:4)

##estimate a log-variance specification with a log-ARCH(4)
##structure, a log-HARCH(5) term and a first-order asymmetry/leverage
##term:
larch(e, arch=1:4, harch=5, asym=1)

##estimate a log-variance specification with a log-ARCH(4)
##structure, an asymmetry/leverage term, a 10-period log(EWMA) as
##volatility proxy, and the log of the squareds of the conditioning
##regressors in the log-variance specification:
larch(e, arch=1:4, asym=1, log.ewma=list(length=10), vxreg=log(x^2))

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
