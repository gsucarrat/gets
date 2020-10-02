\name{gmm}
\alias{gmm}

\title{Generalised Method of Moment (GMM) estimation of linear models}

\description{
Generalised Method of Moment (GMM) estimation of linear models with either ordinary (homoscedastic error) or robust (heteroscedastic error) coefficient-covariance, see Hayashi (2000) chapter 3.
}

\usage{
gmm(y, x, z, tol = .Machine$double.eps,
  weighting.matrix = c("efficient", "2sls", "identity"),
  vcov.type = c("ordinary", "robust"))
}

\arguments{
  \item{y}{numeric vector, the regressand}
  \item{x}{numeric matrix, the regressors}
  \item{z}{numeric matrix, the instruments}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the matrices that are inverted, see the \code{\link{solve}} function}
  \item{weighting.matrix}{a character that determines the weighting matrix to bee used, see "details"}
  \item{vcov.type}{a character that determines the expression for the coefficient-covariance, see "details"}
}

\details{
\code{weighting.matrix = "identity"} corresponds to the Instrumental Variables (IV) estimator, \code{weighting.matrix = "2sls"} corresponds to the 2 Stage Least Squares (2SLS) estimator, whereas \code{weighting.matrix = "efficient"} corresponds to the efficient GMM estimator, see chapter 3 in Hayashi(2000).

\code{vcov.type = "ordinary"} returns the ordinary expression for the coefficient-covariance, which is valid under conditionally homoscedastic errors. \code{vcov.type = "robust"} returns an expression that is also valid under conditional heteroscedasticity, see chapter 3 in Hayashi (2000).
}

\value{A list with, amongst other, the following items:

\item{n}{number of observations}
\item{k}{number of regressors}
\item{df}{degrees of freedom, i.e. n-k}
\item{coefficients}{a vector with the coefficient estimates}
\item{fit}{a vector with the fitted values}
\item{residuals}{a vector with the residuals}
\item{residuals2}{a vector with the squared residuals}      
\item{rss}{the residual sum of squares}
\item{sigma2}{the regression variance}
\item{vcov}{the coefficient-covariance matrix}
\item{logl}{the normal log-likelihood}

}

\references{
F. Hayashi (2000): 'Econometrics'. Princeton: Princeton University Press
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{solve}}, \code{\link{ols}} }

\examples{

##generate data where regressor is correlated with error:
set.seed(123) #for reproducibility
n <- 100
z1 <- rnorm(n) #instrument
eps <- rnorm(n) #ensures cor(z,eps)=0
x1 <- 0.5*z1 + 0.5*eps #ensures cor(x,eps) is strong
y <- 0.4 + 0.8*x1 + eps #the dgp
cor(x1, eps) #check correlatedness of regressor
cor(z1, eps) #check uncorrelatedness of instrument

x <- cbind(1,x1) #regressor matrix
z <- cbind(1,z1) #matrix with instruments

##efficient gmm estimation:
mymod <- gmm(y, x, z)
mymod$coefficients

##ols (for comparison):
mymod <- ols(y,x)
mymod$coefficients

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
