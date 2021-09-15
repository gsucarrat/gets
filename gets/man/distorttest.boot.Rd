\name{distorttest.boot}
%\alias{as.arx.lm}

\title{Bootstrapped Jiao-Pretis-Schwarz Outlier Distortion Test}

\description{
  Test whether Outliers distort regression coefficients using Bootstrapping.
  %The function \code{as.arx} is a generic function and its methods returns an object of class \code{\link{arx}}.
}

\usage{
distorttest(x, nboot = 199, clean.sample = TRUE, parametric = FALSE, scale.t.pval = 1, parallel.options = NULL, ...)}

\arguments{
  \item{x}{object of class \code{\link{isat}} }
  \item{nboot}{numeric. Number of bootstrap replications. Default is 199 replications.}
  \item{clean.sample}{logical. Whether the sample should be cleaned of outliers before resampling.}
  \item{parametric}{logical.}
  \item{scale.t.pval}{logical. Default is 1.}
  \item{parallel.options}{NULL (Default) or an integer, i.e. the number of cores/threads to be used for parallel computing (implemented w/makeCluster and parLapply).}
  \item{...}{Further arguments passed to \code{\link{isat}}.}}

\value{A list including an object of class \code{\link{h-test}}.}

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
  distorttest.boot(nile, nboot = 10)
  
  
  # Another example with co-variates
  dat <- hpdata[,c("GD", "GNPQ", "FSDJ")]
  Y <- ts(dat$GD,start = 1959, frequency = 4)
  mxreg <- ts(dat[,c("GNPQ","FSDJ")],start = 1959, frequency = 4)
  m1 <- isat(y = Y, mc = TRUE, sis = FALSE, iis = TRUE)
  m2 <- isat(y = Y, mc = TRUE, sis = FALSE, iis = TRUE, ar = 1)
  m3 <- isat(y = Y, mxreg = mxreg, mc = TRUE, sis = FALSE, iis = TRUE)
  m4 <- isat(y = Y, mxreg = mxreg, mc = TRUE, sis = FALSE, iis = TRUE, ar = 1)
  distorttest(m1, coef = "all")
  distorttest(m2, coef = "GNPQ")
  distorttest(m3, coef = "GNPQ")
  distorttest(m4, coef = c("AR1", "FSDJ"))
  
  # bootstrap
  db1 <- distorttest.boot(m1, nboot = 10)
  db2 <- distorttest.boot(m2, nboot = 10)
  db3 <- distorttest.boot(m3, nboot = 10)
  db4 <- distorttest.boot(m4, nboot = 10)
  db4.clean <- distorttest.boot(m4, nboot = 10, clean.sample = TRUE)

  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
