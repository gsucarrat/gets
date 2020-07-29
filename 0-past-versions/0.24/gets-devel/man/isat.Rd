\name{isat}
\alias{isat}

\title{Indicator Saturation}

\description{
The \code{\link{isat}} function undertakes multi-path indicator saturation to detect outliers and mean-shifts using impulses (IIS), step-shifts (SIS), or trend-indicators (TIS). Indicators are partitioned into blocks and selected over at a chosen level of significance (\code{t.pval}) using the \code{\link{getsm}} function.
}

\usage{
isat(y, mc=TRUE, ar=NULL, ewma=NULL, mxreg=NULL, iis=FALSE, sis=TRUE,
  tis=FALSE, uis=FALSE, blocks=NULL, ratio.threshold=0.8, max.block.size=30,
  t.pval=0.001, wald.pval=t.pval,
  vcov.type= c("ordinary","white","newey-west"), do.pet=FALSE, ar.LjungB=NULL,
  arch.LjungB=NULL, normality.JarqueB=NULL, info.method=c("sc","aic","hq"),
  user.diagnostics=NULL, user.estimator=NULL, gof.function=NULL,
  gof.method = c("min", "max"), include.gum=NULL, include.1cut=FALSE,
  include.empty=FALSE, max.paths=NULL, parallel.options=NULL, turbo=FALSE,
  tol=1e-07, LAPACK=FALSE, max.regs=NULL, print.searchinfo=TRUE, plot=NULL,
  alarm=FALSE)
}

\arguments{
  \item{y}{numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{mc}{logical. \code{TRUE} (default) includes an intercept in the mean specification, whereas FALSE does not}
  \item{ar}{integer vector, say, c(2,4) or 1:4. The AR-lags to include in the mean specification}
  \item{ewma}{either NULL (default) or a list with arguments sent to the \code{\link{eqwma}} function. In the latter case a lagged moving average of y is included as a regressor}
  \item{mxreg}{numeric vector or matrix, say, a \code{\link{zoo}} object, of conditioning variables. Note that missing values in the beginning or at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command. Note also that, if both \code{y} and \code{mxreg} are \code{zoo} objects, then their samples are chosen to match}
  \item{iis}{logical. If \code{TRUE}, impulse indicator saturation is performed.}
  \item{sis}{logical. If \code{TRUE}, step indicator saturation is performed.}
  \item{tis}{logical. If \code{TRUE}, trend indicator saturation is performed.}
  \item{uis}{a matrix of regressors, or a list of matrices.}
  \item{blocks}{\code{NULL} (default), an integer (the number of blocks) or a user-specified \code{list} that indicates how blocks should be put together. If \code{NULL}, then the number of blocks is determined automatically}
  \item{ratio.threshold}{Minimum ratio of variables in each block to total observations to determine the block size, default=0.8. Only relevant if blocks = \code{NULL} }
  \item{max.block.size}{Maximum size of block of variables to be selected over, default=30. Block size used is the maximum of given by either the ratio.threshold and max.block.size }
  \item{t.pval}{numeric value between 0 and 1. The significance level used for the two-sided regressor significance t-tests}
  \item{wald.pval}{numeric value between 0 and 1. The significance level used for the Parsimonious Encompassing Tests (PETs)}
  \item{vcov.type}{the type of variance-covariance matrix used. If NULL (default), then the type used is that of the 'arx' object. This can be overridden by either "ordinary" (i.e. the ordinary variance-covariance matrix) or "white" (i.e. the White (1980) heteroscedasticity robust variance-covariance matrix)}
  \item{do.pet}{logical. If \code{TRUE}, then a Parsimonious Encompassing Test (PET) against the GUM is undertaken at each regressor removal for the joint significance of all the deleted regressors along the current path. If FALSE (default), then a PET is not undertaken at each regressor removal. By default, the numeric value is the same as that of \code{t.pval}}
  \item{ar.LjungB}{a two-item list with names \code{lag} and \code{pval}, or NULL (default). In the former case \code{lag} contains the order of the Ljung and Box (1979) test for serial correlation in the standardised residuals, and \code{pval} contains the significance level. If \code{lag=NULL} (default), then the order used is that of the estimated 'arx' object. If \code{ar.Ljungb=NULL}, then the standardised residuals are not checked for serial correlation}
  \item{arch.LjungB}{a two-item list with names \code{lag} and \code{pval}, or NULL (default). In the former case, \code{lag} contains the order of the Ljung and Box (1979) test for serial correlation in the squared standardised residuals, and \code{pval} contains the significance level. If \code{lag=NULL} (default), then the order used is that of the estimated 'arx' object. If \code{arch.Ljungb=NULL}, then the standardised residuals are not checked for ARCH}
  \item{normality.JarqueB}{\code{NULL} (the default) or a value between 0 and 1. In the latter case, a test for non-normality is conducted using a significance level equal to \code{normality.JarqueB}. If \code{NULL}, then no test for non-normality is conducted}
  \item{info.method}{character string, "sc" (default), "aic" or "hq", which determines the information criterion to be used when selecting among terminal models. The abbreviations are short for the Schwarz or Bayesian information criterion (sc), the Akaike information criterion (aic) and the Hannan-Quinn (hq) information criterion}
  \item{user.diagnostics}{\code{NULL} or a \code{\link{list}} with two entries, \code{name} and \code{pval}, see the \code{user.fun} argument in \code{\link{diagnostics}} }
  \item{user.estimator}{\code{NULL} or a \code{\link{list}} with at least one entry, \code{name}, see the \code{user.estimator} argument in \code{\link{getsFun}} }
  \item{gof.function}{\code{NULL} or a \code{\link{list}} with at least one entry, \code{name}, see the \code{user.estimator} argument in \code{\link{getsFun}} }
  \item{gof.method}{\code{NULL} or a \code{character} that determines whether the best Goodness-of-Fit is a minimum or maximum }
  \item{include.gum}{ignored (temporarily deprecated)}
  \item{include.1cut}{logical. If \code{TRUE}, then the 1-cut model is included among the terminal models, if it passes the diagnostic tests, even if it is not equal to one of the terminals. If FALSE (default), then the 1-cut model is not included (unless it is one of the terminals)}
  \item{include.empty}{logical. If \code{TRUE}, then an empty model is included among the terminal models, if it passes the diagnostic tests, even if it is not equal to one of the terminals. If FALSE (default), then the empty model is not included (unless it is one of the terminals)}
  \item{max.paths}{\code{NULL} (default) or an integer indicating the maximum number of paths to search}
  \item{parallel.options}{\code{NULL} or an integer, i.e. the number of cores/threads to be used for parallel computing (implemented w/\code{makeCluster} and \code{parLapply})}
  \item{turbo}{logical. If \code{TRUE}, then (parts of) paths are not searched twice (or more) unnecessarily, thus yielding a significant potential for speed-gain. However, the checking of whether the search has arrived at a point it has already been comes with a slight computational overhead. Accordingly, if \code{turbo=TRUE}, then the total search time might in fact be higher than if \code{turbo=FALSE}. This happens if estimation is very fast, say, less than quarter of a second. Hence the default is \code{FALSE}}
  \item{tol}{numeric value (default = 1e-07). The tolerance for detecting linear dependencies in the columns of the regressors (see \code{\link{qr}} function). Only used if LAPACK is FALSE (default)}
  \item{LAPACK}{logical. If \code{TRUE}, then use LAPACK. If \code{FALSE} (default), then use LINPACK (see \code{\link{qr}} function)}
  \item{max.regs}{integer. The maximum number of regressions along a deletion path. It is not recommended that this is altered}
  \item{print.searchinfo}{logical. If \code{TRUE} (default), then a print is returned whenever simiplification along a new path is started, and whenever regressors are dropped due to exact multicolinearity}
  \item{plot}{NULL or logical. If \code{TRUE}, then the fitted values and the residuals of the final model are plotted after model selection. If NULL (default), then the value set by \code{\link{options}} determines whether a plot is produced or not.}
  \item{alarm}{logical. If \code{TRUE}, then a sound is emitted (in order to alert the user) when the model selection ends}
}

\details{Multi-path indicator saturation using impulses (IIS), step-shifts (SIS), or trend-indicators (TIS). Indicators are partitioned into sequential blocks (as of beta version 0.7) where the block intervals are defined by the ratio of variables to observations in each block and a specified maximum block size. Indicators are selected over using the \code{\link{getsm}} function. Retained indicators in each block are combined and re-selected over. Fixed covariates that are not selected over can be included in the regression model either in the mxreg matrix, or for auto-regressive terms through the ar specification. See Hendry, Johansen and Santos (2007) and Castle, Doornik, Hendry, and Pretis (2015)}

\value{A list of class 'gets'}

\references{
Castle, Jennifer, L., Doornik, Jurgen, A., Hendry, David F., and Pretis, Felix (2015): 'Detecting Location Shifts during Model Selection by Step-Indicator Saturation', Econometrics, vol 3:2, 240-264.

Hendry, David, F., Johansen, Soren, and Santos, Carlos (2007):  'Automatic selection of indicators in a fully saturated regression'. Computational Statistics, vol 23:1, pp.317-335.

Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44
}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{
Extraction functions for 'gets' objects: \code{\link{coef.gets}}, \code{\link{fitted.gets}}, \code{\link{paths}}, \code{\link{plot.gets}}, \code{\link{print.gets}},\cr
\code{\link{residuals.gets}}, \code{\link{summary.gets}}, \code{\link{terminals}}, \code{\link{vcov.gets}}\cr

Related functions: \code{\link{arx}}, \code{\link{eqwma}}, \code{\link{leqwma}}, \code{\link{zoo}}, \code{\link{getsFun}}
}

\examples{
##SIS using the Nile data
data(Nile)
isat(Nile, sis=TRUE, iis=FALSE, plot=TRUE, t.pval=0.005)

##SIS using the Nile data in an autoregressive model
#isat(Nile, ar=1:2, sis=TRUE, iis=FALSE, plot=TRUE, t.pval=0.005)

##HP Data
##load Hoover and Perez (1999) data:
#data(hpdata)

##make quarterly data-matrix of zoo type
##(GCQ = personal consumption expenditure):
#y <- zooreg(hpdata$GCQ, 1959, frequency=4)

##transform data to log-differences:
#dlogy <- diff(log(y))

##run isat with step impulse saturation on four
##lags and a constant 1 percent significance level:
#isat(dlogy, ar=1:4, sis=TRUE, t.pval =0.01)

##Example with additional covariates entering through mxreg:

##(GYDQ = disposable personal income):
#x <- zooreg(hpdata$GYDQ, 1959, frequency=4)

##transform data to log-differences:
#dlogx <- diff(log(x))

##run isat with step impulse saturation on four
##lags and a constant 1 percent significance level:
#isat(dlogy, mxreg=dlogx, ar=1:4, sis=TRUE, t.pval =0.01)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
