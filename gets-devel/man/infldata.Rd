\name{infldata}
\alias{infldata}
\docType{data}

\title{Quarterly Norwegian year-on-year CPI inflation}

\description{
Quarterly Norwegian year-on-year CPI inflation from 1989(1) to 2015(4).
}

\usage{data("infldata")}

\format{
  A data frame with 108 observations on the following 5 variables:
  \describe{
    \item{\code{date}}{a factor containing the dates}
    \item{\code{infl}}{year-on-year inflation}
    \item{\code{q2dum}}{a dummy variable equal to 1 in quarter 2 and 0 otherwise}
    \item{\code{q3dum}}{a dummy variable equal to 1 in quarter 3 and 0 otherwise}
    \item{\code{q4dum}}{a dummy variable equal to 1 in quarter 4 and 0 otherwise}
  }
}

\source{
Statistics Norway (SSB): \url{http://www.ssb.no/}. The raw data comprise monthly CPI data obtained via \url{https://www.ssb.no/statbank/table/08183}.
}

\references{
Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44
}

\examples{
data(infldata)
infldata <- zooreg(infldata[,-1], frequency=4, start=c(1989,1))
plot(infldata[,"infl"])
}

\keyword{datasets}
\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
