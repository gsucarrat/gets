\name{boot.distorttest}
\alias{boot}
\alias{boot.distorttest}
\alias{print.boot.distorttest}

\title{Bootstrapped Jiao-Pretis-Schwarz Outlier Distortion Test}

\description{
  Implements the Jiao-Pretis-Schwarz bootstrap test for coefficient distortion due to outliers by comparing coefficient estimates obtained using OLS to estimates obtained using the robust IIS estimator implemented using \code{isat}. Three bootstrap schemes are available – using the original sample (not recommended), the clean (outlier-removed) data, and using the clean (outlier-removed) sample with scaled cut-offs used to detect outliers in IIS implemented using isat. See the referenced Jiao-Pretis-Schwarz Paper below for more information.
  %The function \code{as.arx} is a generic function and its methods returns an object of class \code{\link{arx}}.
}

\usage{
boot(x, ...)

##S3 method for objects of class 'distorttest':
\method{boot}{distorttest}(x, nboot, clean.sample = TRUE, 
  parametric = FALSE, scale.t.pval = 1, parallel.options = NULL, 
  quantiles = c(0.90, 0.95, 0.99), ...)

##S3 printing method for objects of class 'boot.distorttest':
\method{print}{boot.distorttest}(x, print.proportion = FALSE, ...)
}

\arguments{
  \item{x}{object of class \code{\link{isat}} }
  \item{nboot}{numeric. Number of bootstrap replications. A high number of replications are recommended for final estimation (more than 200 at least).}
  \item{clean.sample}{logical. Whether the outlier-removed sample should be used in resampling.}
  \item{parametric}{logical. Whether to use a parametric bootstrap. Default is non-parametric (FALSE). Parametric currently not implemented for autoregressive models.}
  \item{scale.t.pval}{numeric. Scaled target p-value (for selection) relative to the initial p-value used in isat. Default is 1. E.g. a value of 0.5 would scale an initial target p-value of 0.05 to 0.025.}
  \item{parallel.options}{NULL (Default) or an integer, i.e. the number of cores/threads to be used for parallel computing (implemented w/makeCluster and parLapply).}
  \item{print.proportion}{logical. Should the bootstraped Jiao-Pretis Outlier Proportion Test be printed. Default is FALSE.}
  \item{quantiles}{numeric vector. Quantiles to be shown based on the bootstrapped results. Default is c(0.90, 0.95, 0.99).}
  \item{...}{Further arguments passed to \code{\link{isat}}.}}

\value{A list including an object of class \code{h-test}.}

\references{
Xiyu Jiao, Felix Pretis,and Moritz Schwarz. Testing for Coefficient Distortion due to Outliers with an Application to the Economic Impacts of Climate Change. Available at SSRN: \url{https://ssrn.com/abstract=3915040 or http://dx.doi.org/10.2139/ssrn.3915040}
}

\author{
  Xiyu Jiao \url{https://sites.google.com/view/xiyujiao}\cr\cr
  Felix Pretis \url{https://www.felixpretis.org/}\cr\cr
  Moritz Schwarz \url{https://moritzschwarz.org}\cr\cr
}

\seealso{
  \code{\link{isat}, \link{distorttest}}
}


\examples{
  \dontrun{
  data(Nile)
  nile <- isat(Nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.01)
  
  distorttest(nile)
  # bootstrap (with nboot = 5 to save time. Higher replications are recommended)
  boot(distorttest(nile), nboot = 5)
  
  data("hpdata")
  # Another example with co-variates
  dat <- hpdata[,c("GD", "GNPQ", "FSDJ")]
  Y <- ts(dat$GD,start = 1959, frequency = 4)
  mxreg <- ts(dat[,c("GNPQ","FSDJ")],start = 1959, frequency = 4)
  m1 <- isat(y = Y, mc = TRUE, sis = FALSE, iis = TRUE)
  m2 <- isat(y = Y, mc = TRUE, sis = FALSE, iis = TRUE, ar = 1)
  m3 <- isat(y = Y, mxreg = mxreg, mc = TRUE, sis = FALSE, iis = TRUE)
  m4 <- isat(y = Y, mxreg = mxreg, mc = TRUE, sis = FALSE, iis = TRUE, ar = 1, t.pval = 0.01)
  distorttest(m1, coef = "all")
  distorttest(m2, coef = "all")
  distorttest(m3, coef = "GNPQ")
  distorttest(m4, coef = c("ar1", "FSDJ"))
  
  # bootstrap (with nboot = 5 to save time)
  boot(distorttest(m1), nboot = 5)
  boot(distorttest(m2), nboot = 5)
  boot(distorttest(m3), nboot = 5)
  boot(distorttest(m4), nboot = 5)
  boot(distorttest(m4), nboot = 5, parametric = TRUE, scale.t.pval = 0.5)

  }
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
