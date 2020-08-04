\name{predict.arx}
\alias{predict.arx}

\title{Forecasting with 'arx' models}

\description{Generate out-of-sample forecasts up to n steps ahead for objects of class \code{\link{arx}}. Optionally, quantiles of the forecasts are also returned if any of the arguments \code{ci.levels} or \code{probs} are specified. The forecasts, confidence intervals and quantiles are obtained via simulation. By default, 5000 simulations is used, but this can be changed via the \code{n.sim} argument. Also by default, the simulations uses a classical bootstrap to sample from the standardised residuals. To use an alternative set of standardised innovations, for example the standard normal, use the \code{innov} argument. If \code{plot=TRUE}, then a plot of the forecasts is created.}

\usage{
  \method{predict}{arx}(object, spec=NULL, n.ahead=12, newmxreg=NULL,
    newvxreg=NULL, newindex=NULL, n.sim=5000, innov=NULL, probs=NULL,
    ci.levels=NULL, quantile.type=7, return=TRUE, verbose=FALSE,
    plot=NULL, plot.options=list(), ...)
}

\arguments{
  \item{object}{an object of class 'arx'}
  \item{spec}{\code{NULL} (default), \code{"mean"}, \code{"variance"} or \code{"both"}. If \code{NULL}, then it is automatically determined whether information pertaining to the mean or variance specification should be returned}
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated (the default is 12)}
  \item{newmxreg}{a \code{matrix} of \code{n.ahead} rows and \code{NCOL(mxreg)} columns with the out-of-sample values of the \code{mxreg} regressors}
  \item{newvxreg}{a \code{matrix} of \code{n.ahead} rows and \code{NCOL(vxreg)} columns with the out-of-sample values of the \code{vxreg} regressors}
  \item{newindex}{\code{NULL} (default) or the date-index for the \code{zoo} object returned by \code{predict.arx}. If \code{NULL}, then the function uses the in-sample \code{index} to generate the out-of-sample index}
  \item{n.sim}{\code{integer}, the number of replications used for the generation of the forecasts}
  \item{innov}{\code{NULL} (default) or a vector of length \code{n.ahead * n.sim} containing the standardised errors (that is, zero mean and unit variance) used for the forecast simulations. If \code{NULL}, then a classica bootstrap procedure is used to draw from the standardised in-sample residuals}
  \item{probs}{\code{NULL} (default) or a \code{vector} with the quantile-levels (values strictly between 0 and 1) of the forecast distribution. If \code{NULL}, then no quantiles are returned unless \code{ci.levels} is non-\code{NULL}}
  \item{ci.levels}{\code{NULL} (default) or a \code{vector} with the confidence levels (expressed as values strictly between 0 and 1) of the forecast distribution. The upper and lower values of the confidence interval(s) are returned as quantiles}
  \item{quantile.type}{an integer between 1 and 9 that selects which algorithm to be used in computing the quantiles, see the argument \code{type} in \code{\link{quantile}} }
  \item{return}{logical. If \code{TRUE} (default), then the out-of-sample predictions are returned. The value \code{FALSE}, which does not return the predictions, may be of interest if only a prediction plot is of interest}
  \item{verbose}{logical with default \code{FALSE}. If \code{TRUE}, then additional information (typically the quantiles and/or the simulated series) used in the generation of forecasts is returned. If \code{FALSE}, then only the forecasts are returned}
  \item{plot}{\code{NULL} (default) or logical. If \code{NULL}, then the value set by \code{options$plot} (see \code{\link{options}}) determines whether a plot is produced or not. If \code{TRUE}, then the out-of-sample forecasts are plotted.}
  \item{plot.options}{a \code{list} of options related to the plotting of forecasts, see 'Details'}
  \item{\dots}{additional arguments}
}

\value{
a \code{vector} of class \code{\link{zoo}} containing the out-of-sample forecasts, or a \code{matrix} of class \code{\link{zoo}} containing the out-of-sample forecasts together with prediction-quantiles, or \code{NULL} if \code{return=FALSE}
}

\details{
The \code{plot.options} argument is a \code{list} that, optionally, can contain any of the following arguments:\cr

\itemize{
  \item \code{keep}: integer greater than zero (the default is 12) that controls the number of in-sample actual values to plot\cr
  \item \code{line.at.origin}: logical. If \code{TRUE}, then a vertical line is drawn at the forecast origin, that is, at the last in-sample observation\cr
  \item \code{start.at.origin}: logical. If \code{TRUE}, then the drawing of the forecast line starts at the actual value of the forecast origin\cr
  \item \code{dot.at.origin}: logical. If \code{TRUE}, then a dot is drawn at the forecast origin\cr
  \item \code{hlines}: numeric vector that indicates where to draw grey horisontal grid lines\cr
  \item \code{col}: numeric vector of length two that controls the colour of the plotted lines. The first value controls the colour of the forecasts and the fitted values, whereas the second controls the colour of the actual values\cr
  \item \code{lty}: numeric vector of length two that controls the line type. The first value controls the line type of the forecast, whereas the second controls the line type of the actual and fitted values\cr
  \item \code{lwd}: an integer that controls the width of the plotted lines (the default is 1)\cr
  \item \code{ylim}: numeric vector of length two that contains the limits of the y-axis of the prediction plot\cr
  \item \code{ylab}: a character that controls the text on the y-axis\cr
  \item \code{main}: a character that controls the text in the overall title\cr
  \item \code{legend.text}: a character vector of length two that controls how the forecast and actual lines should be named or referred to in the legend of the plot\cr
  \item \code{fitted}: If \code{TRUE}, then the fitted values as well as actual values are plotted in-sample\cr
  \item \code{newmactual}: numeric vector or \code{NULL} (default). Enables the plotting of actual values out-of-sample in the mean in addition to the forecasts\cr
  \item \code{newvactual}: numeric vector or \code{NULL} (default). Enables the plotting of squared residuals ('actual values') out-of-sample in addition to the forecasts\cr
  \item \code{shades}: numeric vector of length \code{length(ci.levels)} that contains the shades of grey associated with the confidence intervals in the prediction plot. The shades can range from 100 (white) to 0 (black)\cr
}

}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{arx}}}

\examples{
##simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 40)

##estimate AR(2) model with intercept:
mymod <- arx(y, mc=TRUE, ar=c(1,2))

##generate out-of-sample forecasts:
predict(mymod)

##same, but plot the predictions in addition:
#predict(mymod, plot=TRUE)

##same, but return also the quantiles of the confidence intervals:
#predict(mymod, ci.levels=c(0.50,0.90), plot=TRUE)

##same, but with non-default levels on the confidence intervals:
#predict(mymod, ci.levels=c(0.20,0.80, 0.99), plot=TRUE)

##same, but with more confidence intervals:
#predict(mymod, ci.levels=seq(0.20, 0.95, by=0.05), plot=TRUE)

##same, but with less rugged ci's (achieved by increasing n.sim):
##predict(mymod, ci.levels=seq(0.20, 0.95, by=0.05), n.sim=50000, plot=TRUE)

##same, but using standard normals (instead of bootstrap) in the simulations:
#n.sim <- 2000
#n.ahead <- 12 #the default on n.ahead
#predict(mymod, ci.levels=seq(0.20, 0.95, by=0.05), n.sim=n.sim,
#  innov=rnorm(n.ahead*n.sim), plot=TRUE)

##make x-regressors:
x <- matrix(rnorm(40*3), 40, 3)

##estimate AR(1) model with intercept and covariates:
mymod <- arx(y, mc=TRUE, ar=1, mxreg=x)

##predict up to 5 steps ahead, setting x's to 0 out-of-sample:
predict(mymod, n.ahead=5, newmxreg=matrix(0,5,NCOL(x)))

##same, but do also plot:
#predict(mymod, n.ahead=5, newmxreg=matrix(0,5,NCOL(x)),
#  plot=TRUE)

##estimate an AR(2) model w/intercept and a log-ARCH(1) specification
##on the variance:
#mymodel <- arx(y, mc=TRUE, ar=1:2, arch=1)

##generate forecasts of the conditional variances:
#predict(mymodel, spec="variance")

##same, but do also plot:
#predict(mymodel, spec="variance", plot=TRUE)

##illustrations of the usage of plot.options:
#mymodel <- arx(y, mc=TRUE, ar=1)
#predict(mymodel, plot=TRUE, plot.options=list(keep=1))
#predict(mymodel, plot=TRUE, plot.options=list(line.at.origin=TRUE))
#predict(mymodel, plot=TRUE, plot.options=list(start.at.origin=FALSE))
#predict(mymodel, plot=TRUE, 
#  plot.options=list(start.at.origin=FALSE, fitted=TRUE))
#predict(mymodel, plot=TRUE, plot.options=list(dot.at.origin=FALSE))
#predict(mymodel, plot=TRUE, plot.options=list(hlines=c(-2,-1,0,1,2)))
#predict(mymodel, plot=TRUE, plot.options=list(col=c("darkred","green")))
#predict(mymodel, plot=TRUE, plot.options=list(lty=c(3,2)))
#predict(mymodel, plot=TRUE, plot.options=list(lwd=3))
#predict(mymodel, plot=TRUE, plot.options=list(ylim=c(-8,8)))
#predict(mymodel, plot=TRUE, plot.options=list(ylab="User-specified y-axis"))
#predict(mymodel, plot=TRUE, 
#  plot.options=list(main="User-specified overall title"))
#predict(mymodel, plot=TRUE, 
#  plot.options=list(legend.text=c("User-specified 1","User-specified 2")))
#predict(mymodel, plot=TRUE, plot.options=list(fitted=TRUE))
#predict(mymodel, plot=TRUE, plot.options=list(newmactual=rep(0,6)))
#predict(mymodel, plot=TRUE, plot.options=list(shades.of.grey=c(95,50)))
#predict(mymodel, plot=TRUE, plot.options=list(shades.of.grey=c(50,95))) #invert shades

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
