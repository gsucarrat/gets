\name{isatvarcorrect}
\alias{isatvarcorrect}

\title{Consistency and Efficiency Correction for Impulse Indicator Saturation}

\description{
  Takes an \code{\link{isat}} object and corrects the estimates of the error variance and the estimated standard errors of 'forced' regressors.
}

\usage{
  isatvarcorrect(x, mcor=1)
}

\arguments{
  \item{x}{an \code{\link{isat}} object}
  \item{mcor}{integer, number of iterations in the correction. Default = 1.}
}

\details{Impulse indicator saturation results in an under-estimation of the error variance as well as the variance of regressors not selected over. The magnitude of the inconsistency increases with the p-value of selection (\code{t.pval}). The function takes an \code{\link{isat}} object and applies the impulse indicator consistency (\code{\link{isvarcor}}) and efficiency correction (\code{\link{isvareffcor}}) of the estimated error variance and the estimated variance of regressors not selected over. See Johansen and Nielsen (2016a) and (2016b).}

\value{Returns an \code{\link{isat}} object in which the estimated standard errors, t-statistics, p-values, standard error of the regression, and log-likelihood are consistency and efficiency corrected when using impulse indicator saturation (\code{iis=TRUE}).}

\references{

Johansen, S., & Nielsen, B. (2016a). Asymptotic theory of outlier detection algorithms for linear time series regression models. Scandinavian Journal of Statistics, 43(2), 321-348.

Johansen, S., & Nielsen, B. (2016b). Rejoinder: Asymptotic Theory of Outlier Detection Algorithms for Linear. Scandinavian Journal of Statistics, 43(2), 374-381.

Pretis, F., Reade, J., & Sucarrat, G. (2018). Automated General-to-Specific (GETS) regression modeling and indicator saturation methods for the detection of outliers and structural breaks. Journal of Statistical Software, 86(3).

}

\author{
  Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{
  \code{\link{isat}}, \code{\link{isvarcor}}, \code{\link{isvareffcor}}
}

\examples{
###Consistency and Efficiency Correction of Impulse Indicator Estimates
nile <- as.zoo(Nile)
isat.nile <- isat(nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.1)
isat.nile.corrected <- isatvarcorrect(isat.nile)

isat.nile$sigma2
isat.nile.corrected$sigma2
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}