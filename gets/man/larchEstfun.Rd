\name{larchEstfun}
\alias{larchEstfun}

\title{Estimation of a log-variance model}

\description{
Two-step estimation of a log-variance model: OLS in step 1, bias correction w/residuals in step 2 (see the code for details). The function \code{larchEstfun()} is not intended for the average user, but is called by \code{\link{larch}} and \code{\link{gets.larch}}.
}

\usage{
larchEstfun(loge2, x, e, vcov.type = c("robust", "hac"), tol = 1e-07)
}

\arguments{
  \item{loge2}{numeric vector, the log of the squared errors 'e' (adjusted for zeros on e, if any) }
  \item{x}{numeric matrix, the regressors}
  \item{e}{numeric vector, the errors}
  \item{vcov.type}{\code{character} vector, "robust" (default) or "hac". If "robust", then the White (1980) heteroscedasticity-robust variance-covariance matrix is used for inference. If "hac", then the Newey and West (1987) heteroscedasticity and autocorrelation-robust matrix is used}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the regressors in the first step estimation by OLS, see \code{\link{ols}}. Only used if \code{LAPACK} is \code{FALSE} }
}

\details{No details for the moment.}

\value{A \code{\link{list}}.}

\references{No references for the moment.}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{qr}}, \code{\link{larch}}, \code{\link{gets.larch}} }

\examples{
##no examples for the moment
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
