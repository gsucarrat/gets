\name{isattest}
\alias{isattest}

\title{Indicator Saturation Test}

\description{
  Takes an 'isat' object returned by the \code{isat} function as input and returns the results of a hypothesis test on the time-varying intercept or long-run equilibrium against a specified null-hypothesis for a chosen level of significance - see Pretis (2015).
}

\usage{
  isattest(x, hnull=0, lr=FALSE, ci.pval=0.99, plot=NULL, plot.turn=FALSE,
    conscorr=FALSE, effcorr=FALSE, mcor = 1, biascorr=FALSE, mxfull = NULL,
    mxbreak=NULL)
}

\arguments{
  \item{x}{a 'gets' object obtained with the \code{\link{isat}} function}
  \item{hnull}{numeric. the null-hypothesis value to be tested against.}
  \item{lr}{logical. If TRUE and 'x' contains autoregressive elements, then \code{isattest} tests on the long-run equilibrium path. See Pretis (2015).}
  \item{ci.pval}{numeric between 0 and 1. Default is 0.99, the level of significance for the confidence interval of the test against 'hnull'.}
  \item{plot}{logical. If TRUE, then a plot showing the coefficient path and bias relative to 'hnull' is shown.}
  \item{plot.turn}{logical. If TRUE, then the plot output adds the time of the breaks to the plot showing the bias relative to 'hnull'.}
  \item{biascorr}{logical. If TRUE, then the coefficient path is bias-corrected using \code{\link{biascorr}}. This is only valid for the non-dynamic test without additional covariates.}
  \item{conscorr}{logical. If TRUE then the Johansen and Nielsen (2016) impulse-indicator consistency correction is applied to estimated residual variance.}
  \item{effcorr}{logical. If TRUE then the Johansen and Nielsen (2016) m-step efficiency correction is applied to estimated standard errors of `fixed' regressors.}
  \item{mcor}{integer. The m-step efficiency correction factor, where m=mcor.}
  \item{mxfull}{string. The name of the full-sample variable when constructing the coefficient path of user-specified break variables.}
  \item{mxbreak}{string. The name of the break variables used to construct the coefficient path of user-specified break variables.}
  
}

\details{The function tests the coefficient path (or long-run equilibrium path) against a specified null hypothesis at a chosen level of significance. If conducted on an \code{\link{isat}} model of a forecast error or relative forecast differential, then this corresponds to the test of time-varying predictive accuracy of Pretis (2015). The resulting output plot shows the coefficient path in the top panel (where 'hnull' is plotted as dotted lines), with the bias (significant difference relative to 'hnull') in the lower panel. If \code{mxfull} and \code{mxbreak} are specified, then the function tests on the coefficient path of the user-specified variable, where \code{mxfull} denotes the ful-sample variable name, to which the \code{mxbreak} variables are added. To correct for the under-estimation of the residual variance, the argument \code{conscorr} implements the Johansen and Nielsen (2016) consistency correction, and \code{effcorr} adds the efficiency correction for standard errors on fixed regressors which are not selected over.}

\value{A Tx4 matrix (with T = number of observations) where the first two columns denote the confidence interval of the coefficient path (or the long-run equilibrium path if 'lr=TRUE'). The third and fourth column denote the bias of the coefficient path relative to the chosen null-hypothesis, where 'bias.high' denotes the bias when the series tested is above the hypothesized value, and 'bias.low' denotes the bias when the series tested is significantly below the hypothesized value.}

\references{
Johansen, S., & Nielsen, B. (2016): 'Asymptotic theory of outlier detection algorithms for linear time series regression models.' Scandinavian Journal of Statistics, 43(2), 321-348.

Pretis, F. (2015): 'Testing for time-varying predictive accuracy using bias-corrected indicator saturation'. Oxford Department of Economics Discussion Paper.

Hendry, David, F., Johansen, Soren, and Santos, Carlos (2007):  'Automatic selection of indicators in a fully saturated regression'. Computational Statistics, vol 23:1, pp.317-335.
}

\author{
  Felix Pretis, \url{http://www.felixpretis.org/}\cr
}

\seealso{
  \code{\link{isat}}, \code{\link{coef.gets}}, \code{\link{plot.gets}}, \code{\link{biascorr}}, \code{\link{isatvar}}
}

\examples{
##Using artificial data:
#set.seed(123)
#d <- matrix(0,100,1)
#d[35:55] <- 1
#e <- rnorm(100, 0, 1)
#y <- d*2  +e 
#plot(y, type="l")

##Static Test against hnull=0 using bias-correction:

#ys <- isat(y, sis=TRUE, iis=FALSE, tis=FALSE, t.pval=0.01)
#isattest(ys, hnull=0, lr=FALSE, ci.pval = 0.99, plot.turn = FALSE, biascorr=TRUE)

##Dynamic Test of the long-run equilibrium against hnull=2 with breakpoints
##labelled in the plot:

#ys <- isat(y, sis=TRUE, iis=FALSE, tis=FALSE, t.pval=0.01, ar=1:2)
#isattest(ys, hnull=2, lr=TRUE, ci.pval = 0.99, plot.turn = TRUE, biascorr=FALSE)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}