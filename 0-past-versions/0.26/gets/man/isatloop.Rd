\name{isatloop}
\alias{isatloop}

\title{Repeated Impulse Indicator Saturation}

\description{
  Runs \code{\link{isat}} repeatedly at pre-specified significance levels to yield multiple iterations used in\cr
  \code{\link{outlierscaletest}}.
}

\usage{
  isatloop(num=c(seq(from=20, to=1, by=-1)), t.pval.spec = FALSE,  
  print=FALSE, y, ar=NULL, iis=TRUE, sis=FALSE, ...)
  
}

\arguments{
  \item{num}{numeric, target expected number of outliers under the null hypothesis, or target proportion of outliers if \code{t.pval.spec==TRUE}}
  \item{t.pval.spec}{logical, if \code{TRUE}, then \code{num} specifies proportion rather than number of targeted outliers}
  \item{print}{logical, if \code{TRUE}, then iterations are printed}
  \item{y}{numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{ar}{integer vector, say, c(2,4) or 1:4. The AR-lags to include in the mean specification}
  \item{iis}{logical, whether to use \code{iis}}
  \item{sis}{logical, whether to use \code{sis}, default is \code{FALSE}}
  \item{...}{any argument from \code{\link{isat}} can also be used in \code{isatloop}}
}

\details{The function repeatedly runs \code{\link{isat}} detecting outliers in a model of \code{y} at different chosen target levels of significance speciefied in \code{num}.  The output of this function is used as the input for the \code{\link{outlierscaletest}} function. All additional arguments from \code{\link{isat}} can be passed to \code{isatloop}.}

\value{Returns a list of two items. The first item is the number of observations. The second item is a dataframe containing the expected and observed proportion (and number of outliers) for each specified significance level of selection.}

\references{
  
  Jiao, X. & Pretis, F. (2019). Testing the Presence of Outliers in Regression Models. Discussion Paper.
  
  Pretis, F., Reade, J., & Sucarrat, G. (2018). Automated General-to-Specific (GETS) regression modeling and indicator saturation methods for the detection of outliers and structural breaks. Journal of Statistical Software, 86(3).
  
}

\author{
  Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{
  \code{\link{isat}}, \code{\link{outlierscaletest}}
}

\examples{
  ###Repeated isat models using the Nile dataset
  ### where p-values are chosen such that the expected number of outliers under the null
  ### corresponds to 1, 2, 3, 4 and 5.
  nile <- as.zoo(Nile)
  isat.nile.loop <- isatloop(y=nile, iis=TRUE, num=c(1,2, 3, 4, 5))
  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}