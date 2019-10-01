\name{outliertest}
\alias{outliertest}

\title{Jiao and Pretis Outlier Proportion and Count Tests}

\description{
  Tests whether the proportion (or number) of outliers detected using impulse indicator saturation is different from the proportion (or number) of outliers expected under the null hypothesis of no outliers using the Jiao and Pretis (2019) proportion and count outlier tests.
}

\usage{
  outliertest(x, noutl=NULL, t.pval=NULL, T=NULL, 
  m=1, infty=FALSE, alternative="two.sided")
  
}

\arguments{
  \item{x}{an \code{\link{isat}} object}
  \item{noutl}{integer, number of detected outliers if no \code{\link{isat}} object is provided i.e. x=NULL}
  \item{t.pval}{numeric, between 0 and 1. Selection p-value used in indicator saturation if no \code{\link{isat}} object is provided i.e. x=NULL}
  \item{T}{integer, sample sized used in indicator saturation if no \code{\link{isat}} object is provided i.e. x=NULL}
  \item{m}{integer, number of iterations in variance computation, default=1}
  \item{infty}{logical, argument used for variance computation}
  \item{alternative}{"two-sided", "less", "greater-than", alternative hypothesis of outlier test.}
}

\details{The function computes the estimated proportion of outliers (gauge) based on impulse indicator saturation and constructs the proportion and count outlier test statistics from Jiao and Pretis (2019). The null hypothesis is that the proportion (or count) of outliers is not different than the proportion (or count) of outliers detected under the null hypothesis of no outliers. The first test compares the estimated proportion of outliers scaled by its estimated variance against a standard normal distribution. The second test compares the number of outliers against a Poisson distribution.
  If an \code{\link{isat}} object is provided in \code{x}, then the function automatically extracts the detected impulses and computes the estimated outlier proportion. If no \code{\link{isat}} object is provided and \code{x=NULL}, then the tests can be conducted manually by providing the number of detected outliers (\code{noutl}), the sample size (\code{T}), and the chosen level of signficance used to detect outliers (\code{t.pval}).}

\value{Returns a list of two \code{htest} objects. The first providing the results of the test on the proportion of outliers against a standard normal distribution. The second providing the results on the number of outliers against the Poisson distribution.}

\references{
  
  Jiao, X. & Pretis, F. (2019). Testing the Presence of Outliers in Regression Models. Discussion Paper.
  
  Pretis, F., Reade, J., & Sucarrat, G. (2018). Automated General-to-Specific (GETS) regression modeling and indicator saturation methods for the detection of outliers and structural breaks. Journal of Statistical Software, 86(3).
  
}

\author{
  Xiyu Jiao, & Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{
  \code{\link{isat}}
}

\examples{
  ###Testing the Presence of Outliers in the Nile Data
  nile <- as.zoo(Nile)
  isat.nile <- isat(nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.1)
  outliertest(isat.nile)
  
  ###Testing the number of outliers when the sample is T=200, 
  ### with 7 detected outliers at t.pval=0.05 if no isat object is provided:
  outliertest(x=NULL, noutl=7, t.pval=0.05, T=200)
  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}