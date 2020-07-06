\name{vargaugeiis}
\alias{vargaugeiis}

\title{Variance of the Impulse Indicator Saturation Gauge}

\description{
  Computes the variance of the gauge (false-positive rate of outliers under the null of no outliers) in impulse indicator saturation based on Jiao and Pretis (2019).
}

\usage{
  vargaugeiis(t.pval, T, infty=FALSE, m=1)
  
}

\arguments{
  \item{t.pval}{numeric, between 0 and 1. Selection p-value used in indicator saturation.}
  \item{T}{integer, sample sized used in indicator saturation.}
  \item{m}{integer, number of iterations in variance computation, default=1}
  \item{infty}{logical, argument used for variance computation}
}

\details{The function computes the variance of the Gauge (false-positive rate of outliers in impulse indicator saturation) for a given level of significance of selection (\code{t.pval}) and sample size ({\code{T}) based on Jiao and Pretis (2019). This is an auxilliary function used within the \code{\link{outliertest}} function.}}

\value{Returns a dataframe of the variance and standard deviation of the gauge, as well the asymptotic variance and standard deviation.}

\references{
  
  Jiao, X. & Pretis, F. (2019). Testing the Presence of Outliers in Regression Models. Discussion Paper.
  
  Pretis, F., Reade, J., & Sucarrat, G. (2018). Automated General-to-Specific (GETS) regression modeling and indicator saturation methods for the detection of outliers and structural breaks. Journal of Statistical Software, 86(3).
  
}

\author{
  Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{
  \code{\link{isat}}, \code{\link{outliertest}}
}

\examples{
  ###Computing the variance of the gauge under the null for a sample of T=200 observations:
  vargaugeiis(t.pval=0.05, T=200, infty=FALSE, m=1)
  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}