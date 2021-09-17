\name{boot.distorttest}
\alias{boot.distorttest}
\alias{boot}
\alias{print.boot.distorttest}

\title{Bootstrapped Jiao-Pretis-Schwarz Outlier Distortion Test}

\description{
  Test whether Outliers distort regression coefficients using Bootstrapping.
  %The function \code{as.arx} is a generic function and its methods returns an object of class \code{\link{arx}}.
}

\usage{
boot(x, ...)

##S3 method for objects of class 'distorttest':
\method{boot}{distorttest(x, nboot = 199, clean.sample = TRUE, 
parametric = FALSE, scale.t.pval = 1, parallel.options = NULL, ...)}
##S3 printing method for objects of class 'boot.distorttest':
\method{print}{boot.distorttest}(x, ...)
}

\arguments{
  \item{x}{object of class \code{\link{isat}} }
  \item{nboot}{numeric. Number of bootstrap replications. Default is 199 replications.}
  \item{clean.sample}{logical. Whether the outlier-removed sample should be used in resampling.}
  \item{parametric}{logical. Whether to use a parametric bootstrap. Default is non-parametric (FALSE). Parametric currently not implemented for autoregressive models.}
  \item{scale.t.pval}{numeric. Scaled target p-value. Default is 1.}
  \item{parallel.options}{NULL (Default) or an integer, i.e. the number of cores/threads to be used for parallel computing (implemented w/makeCluster and parLapply).}
  \item{...}{Further arguments passed to \code{\link{isat}}.}}

\value{A list including an object of class \code{h-test}.}

\author{
  Felix Pretis \url{https://felixpretis.climateeconometrics.org/}\cr\cr
  Moritz Schwarz \url{https://moritzschwarz.org}\cr\cr
}

\seealso{
  \code{\link{isat}, \link{distorttest}}
}

\examples{
  
  data(Nile)
  Nile[10] <- Nile[10]+1000 # to ensure there is an indicator retained
  nile <- isat(Nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.005)
  
  distorttest(nile)
  # bootstrap (with nboot = 2 to save time)
  boot(distorttest(nile), nboot = 2)
  
  data("hpdata")
  # Another example with co-variates
  dat <- hpdata[,c("GD", "GNPQ", "FSDJ")]
  Y <- ts(dat$GD,start = 1959, frequency = 4)
  mxreg <- ts(dat[,c("GNPQ","FSDJ")],start = 1959, frequency = 4)
  m1 <- isat(y = Y, mc = TRUE, sis = FALSE, iis = TRUE)
  m2 <- isat(y = Y, mc = TRUE, sis = FALSE, iis = TRUE, ar = 1)
  m3 <- isat(y = Y, mxreg = mxreg, mc = TRUE, sis = FALSE, iis = TRUE)
  m4 <- isat(y = Y, mxreg = mxreg, mc = TRUE, sis = FALSE, iis = TRUE, ar = 1)
  distorttest(m1, coef = "all")
  distorttest(m2, coef = "all")
  distorttest(m3, coef = "GNPQ")
  distorttest(m4, coef = c("ar1", "FSDJ"))
  
  # bootstrap (with nboot = 2 to save time)
  boot(distorttest(m1), nboot = 2)
  boot(distorttest(m2), nboot = 2)
  boot(distorttest(m3), nboot = 2)
  boot(distorttest(m4), nboot = 2)
  boot(distorttest(m4), nboot = 2, clean.sample = TRUE)

  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
