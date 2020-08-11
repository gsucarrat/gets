\name{outlierscaletest}
\alias{outlierscaletest}

\title{Sum and Sup Scaling Outlier Tests}

\description{
  Computes the Sum and Supremum Scaling Tests for the overall presence of outliers based on Jiao and Pretis (2019).
}

\usage{
  outlierscaletest(x, nsim = 10000)
  
}

\arguments{
  \item{x}{list, output of the \code{\link{isatloop}} function}
  \item{nsim}{integer, number of replications to simulate critical values for the Sup test}
}

\details{The function takes the output of the \code{\link{isatloop}} function and computes the Scaling Sum and Supremum Tests for the presence of outliers from Jiao and Pretis (2019). The test compares the expected and observed proportion of outliers over the range of different significance levels of selection specified in \code{\link{isatloop}}. The Sum test compares the sum of deviations against the standard normal distribution, the Sup test compares the supremum of deviations against critical values simulated with \code{nsim} replications. The null hypothesis is that the observed proportion of outliers scales with the proportion of outliers under the null of no outliers.}

\value{Returns a list of two \code{htest} objects. The first providing the results of the Sum test on the sum of the deviation of outliers against a standard normal distribution. The second providing the results on the supremum of the deviation of outliers against simulated critical values.}

\references{
  
  Jiao, X. & Pretis, F. (2019). Testing the Presence of Outliers in Regression Models. Discussion Paper.
  
  Pretis, F., Reade, J., & Sucarrat, G. (2018). Automated General-to-Specific (GETS) regression modeling and indicator saturation methods for the detection of outliers and structural breaks. Journal of Statistical Software, 86(3).
  
}

\author{
  Xiyu Jiao, & Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{
  \code{\link{isat}}, \code{\link{isatloop}}
}

\examples{
  ###Repeated isat models using the Nile dataset
  ### where p-values are chosen such that the expected number of outliers under the null
  ### corresponds to 1, 2, ..., 20. Then computing the Outlier Scaling Tests:
  
  #nile <- as.zoo(Nile)
  #isat.nile.loop <- isatloop(y=nile)
  #outlierscaletest(isat.nile.loop)
  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}