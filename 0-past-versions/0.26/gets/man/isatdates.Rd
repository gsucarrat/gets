\name{isatdates}
\alias{isatdates}

\title{Extracting Indicator Saturation Breakdates}

\description{
  Takes an \code{\link{isat}} object and extracts the break dates together with their estimated coefficients.
}

\usage{
  isatdates(x)
}

\arguments{
  \item{x}{an \code{\link{isat}} object}
}

\details{The function extracts the breakdates determined by \code{\link{isat}} for \code{iis}, \code{sis}, and \code{tis}, together with their estimated coefficients and standard errors.}

\value{Returns a list of three elements (one for \code{iis}, \code{sis}, and \code{tis}). Each element lists the name of the break variable, the time index of the break (labelled 'date'), the index of the break date, the estimated coefficient, the standard error of the estimated coefficient, as well as the corresponding t-statistic and p-value.}

\references{

Pretis, F., Reade, J., & Sucarrat, G. (2018). Automated General-to-Specific (GETS) regression modeling and indicator saturation methods for the detection of outliers and structural breaks. Journal of Statistical Software, 86(3).

}

\author{
  Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{
  \code{\link{isat}}
}

\examples{
###Break date extraction of the Nile data
nile <- as.zoo(Nile)
isat.nile <- isat(nile, sis=TRUE, iis=FALSE, plot=TRUE, t.pval=0.005)
isatdates(isat.nile)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}