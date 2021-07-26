\name{logit}
\alias{logit}

\title{Estimation of a logit model}

\description{
Maximum Likelihood (ML) estimation of a logit model.
}

\usage{
logit(y, x, initial.values = NULL, lower = -Inf, upper = Inf, 
    method = 2, lag.length = NULL, control = list(), eps.tol = .Machine$double.eps, 
    solve.tol = .Machine$double.eps )
}

\arguments{
  \item{y}{numeric vector, the binary process}
  \item{x}{numeric matrix, the regressors}
  \item{initial.values}{\code{NULL} or a numeric vector with the initial parameter values passed on to the optimisation routine, \code{\link{nlminb}}. If \code{NULL}, the default, then the values are chosen automatically }
  \item{lower}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}} }
  \item{upper}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}} }
  \item{method}{an integer that determines the expression for the coefficient-covariance, see "details"}
  \item{lag.length}{\code{NULL} or an integer that determines the lag-length used in the robust coefficient covariance. If \code{lag.length} is an integer, then it is ignored unless \code{method = 3} }
  \item{control}{a \code{list} passed on to the control argument of \code{\link{nlminb}} }
  \item{eps.tol}{numeric, a small value that ensures the fitted zero-probabilities are not too small when the log-transformation is applied when computing the log-likelihood }
  \item{solve.tol}{numeric value passed on to the \code{tol} argument of \code{\link{solve}}, which is called whenever the coefficient-coariance matrix is computed. The value controls the toleranse for detecting linear dependence between columns when inverting a matrix }
}

\details{No details for the moment.}

\value{A \code{\link{list}}.}

\references{No references for the moment.}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{nlminb}}, \code{\link{solve}} }

\examples{
##no examples for the moment
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
