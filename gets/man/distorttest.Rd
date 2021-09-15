\name{distorttest}
\alias{distorttest}

\title{Jiao-Pretis-Schwarz Outlier Distortion Test}

\description{
  Test whether Outliers distort regression coefficients. 
  %The function \code{as.arx} is a generic function and its methods returns an object of class \code{\link{arx}}.
}

\usage{
  distorttest(x, coef = "all")
}

\arguments{
  \item{x}{object of class \code{\link{isat}} }
  \item{coef}{Either "all" (Default) to test the distortion on all coefficients or a character vector of explanatory variable names.}
}

\value{Object of class \code{\link{arx}} }

\author{
  Felix Pretis \url{https://felixpretis.climateeconometrics.org/}\cr\cr
  Moritz Schwarz \url{https://moritzschwarz.org}\cr\cr
}

\seealso{
  \code{\link{isat}, \link{distorttest.boot}}
}

\examples{
  
data(Nile)
Nile[10] <- Nile[10]+1000 # to ensure there is an indicator retained
nile <- isat(Nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.005)

distorttest(nile)


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
  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
