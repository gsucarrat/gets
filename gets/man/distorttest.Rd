\name{distorttest}
\alias{distorttest}

\title{Jiao-Pretis-Schwarz Outlier Distortion Test}

\description{
Implements the Jiao-Pretis-Schwarz test for coefficient distortion due to outliers by comparing coefficient estimates obtained using OLS to estimates obtained using the robust IIS estimator implemented using \code{isat}. See the referenced Jiao-Pretis-Schwarz Paper below for more information.
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
Xiyu Jiao, Felix Pretis,and Moritz Schwarz. Testing for Coefficient Distortion due to Outliers with an Application to the Economic Impacts of Climate Change. Available at SSRN: \url{https://www.ssrn.com/abstract=3915040} or \doi{10.2139/ssrn.3915040}
}

\author{
  Xiyu Jiao \url{https://sites.google.com/view/xiyujiao}\cr\cr
  Felix Pretis \url{https://felixpretis.climateeconometrics.org/} \cr\cr
  Moritz Schwarz \url{https://moritzschwarz.org}\cr\cr
}

\seealso{
  \code{\link{isat}, \link{distorttestboot}}
}

\examples{
\dontrun{  
data(Nile)
nile <- isat(Nile, sis=FALSE, iis=TRUE, plot=TRUE, t.pval=0.01)
distorttest(nile)

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
 } 
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
