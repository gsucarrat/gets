\name{sp500data}
\alias{sp500data}
\docType{data}
\title{
Daily Standard and Poor's 500 index data
}
\description{
Daily Standard and Poor's 500 (SP500) index data from 3 January 1950 to 8 March 2016.
}
\usage{data("sp500data")}
\format{
  A data frame with 16652 observations on the following 7 variables:
  \describe{
    \item{\code{Date}}{the dates}
    \item{\code{Open}}{the opening values of the index}
    \item{\code{High}}{the daily maximum value of the index}
    \item{\code{Low}}{the daily minimum value of the index}
    \item{\code{Close}}{the closing values of the index}
    \item{\code{Volume}}{the traded volume}
    \item{\code{Adj.Close}}{the adjusted closing values of the index}
  }
}
\source{
Yahoo Finance, retrieved 9 March 2016
}
\references{
Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44
}
\examples{
data(sp500data)
sp500data <- zoo(sp500data[, -1], order.by = as.Date(sp500data[, "Date"]))
plot(window(sp500data, start = as.Date("2000-01-03")))
}
\keyword{datasets}
