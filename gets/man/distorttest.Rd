\name{distorttest}
\alias{distorttest}

\title{Jiao-Pretis-Schwarz Outlier Distortion Test}

\description{
  Test whether Outliers distort regression coefficients. See the referenced Jiao-Pretis-Schwarz Paper below for more information.
  %The function \code{as.arx} is a generic function and its methods returns an object of class \code{\link{arx}}.
}

\usage{
  distorttest(x, coef = "all")
}

\arguments{
  \item{x}{object of class \code{\link{isat}} }
  \item{coef}{Either "all" (Default) to test the distortion on all coefficients or a character vector of explanatory variable names.}

}

\value{Object of class \code{\link{isat}} }

\references{
Xiyu Jiao, Felix Pretis,and Moritz Schwarz. Testing for Coefficient Distortion due to Outliers with an Application to the Economic Impacts of Climate Change. Available at SSRN: \url{https://ssrn.com/abstract=3915040 or http://dx.doi.org/10.2139/ssrn.3915040}
}

\author{
  Felix Pretis \url{https://felixpretis.climateeconometrics.org/}\cr\cr
  Moritz Schwarz \url{https://moritzschwarz.org}\cr\cr
}

\seealso{
  \code{\link{isat}, \link{boot.distorttest}}
}

\examples{
  
data(Nile)
Nile[10] <- Nile[10]+1000 # to ensure there is an indicator retained
nile <- isat(Nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.005)
distorttest(nile)

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
  
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
