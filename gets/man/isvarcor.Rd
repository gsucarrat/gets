\name{isvarcor}
\alias{isvarcor}

\title{IIS Consistency Correction}

\description{Consistency correction for estimate of residual variance when using impulse indicator saturation.}

\usage{
isvarcor(t.pval, sigma)
}

\arguments{
  \item{t.pval}{numeric value. the p-value of selection in the impulse indicator saturation model.}
  \item{sigma}{numeric value. The estimated standard deviation of the residuals from the impulse indicator saturation model.}
}

\value{a data frame containing the corrected standard deviation \code{$sigma.cor} and the correction factor used \code{$corxi}}

\details{
The Johansen and Nielsen (2016) impulse-indicator consistency correction for the estimated residual standard deviation.
}

\references{
Johansen, S., & Nielsen, B. (2016): 'Asymptotic theory of outlier detection algorithms for linear time series regression models.' Scandinavian Journal of Statistics, 43(2), 321-348.

Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44
}

\author{
 Felix Pretis, \url{http://www.felixpretis.org/}
}

\seealso{\code{\link{isatvar}}}

\examples{

isvarcor(t.pval=0.05, sigma=2)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
