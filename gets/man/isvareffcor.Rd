\name{isvareffcor}
\alias{isvareffcor}

\title{IIS Efficiency Correction}

\description{Efficiency correction for the estimates of coefficient standard errors on fixed regressors.}

\usage{
isvareffcor(t.pval, se, m=1)
}

\arguments{
  \item{t.pval}{numeric value. the p-value of selection in the impulse indicator saturation model.}
  \item{se}{numeric value or vector. The estimated standard errors of the coefficients on fixed regressors in impulse indicator saturation model.}
  \item{m}{integer. The m-step correction factor.}
 }

\value{a data frame containing the corrected standard deviation \code{$se.cor} and the correction factor used \code{$eta.m}}

\details{
The Johansen and Nielsen (2016) impulse-indicator efficiency correction for the estimated standard errors on fixed regressors in impulse indicator models.
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

isvareffcor(t.pval=0.05, se=2, m=1)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
