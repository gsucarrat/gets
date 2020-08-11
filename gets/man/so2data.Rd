\name{so2data}
\alias{so2data}
\docType{data}
\title{
UK SO2 Data
}
\description{
UK Annual Total Anthropogenic Sulphur Dioxide (SO2) Emissions 1946-2005.
}
\usage{data("so2data")}
\format{
  A data frame with 60 observations on the following 4 variables.
  \describe{
    \item{\code{year}}{Year of observation}
    \item{\code{uk_tot_so2}}{UK annual total anthropogenic SO2 emissions in gigagrams}
    \item{\code{Luk_tot_so2}}{Log of UK annual total anthropogenic SO2 emissions}
    \item{\code{DLuk_tot_so2}}{First difference of Log UK annual total anthropogenic SO2 emissions}
  }
}
\details{
Data reports the total estimated anthropogenic SO2 emissions aggregated over coal, petroleum, biomass combustion, smelting, fuel processing, and other processes.
}
\source{
Smith, SJ, J van Aardenne, Z Klimont, RJ Andres, A Volke, and S Delgado Arias. (2011). Anthropogenic Sulfur Dioxide Emissions, 1850-2005: National and Regional Data Set by Source Category, Version 2.86. Data distributed by the NASA Socioeconomic Data and Applications Center (SEDAC), CIESIN, Columbia University, Palisades, New York. Available at

http://sedac.ciesin.columbia.edu/data/set/haso2-anthro-sulfur-dioxide-emissions-1850-2005-v2-86
}
\references{
Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44

Smith, SJ, J van Aardenne, Z Klimont, RJ Andres, A Volke, and S Delgado Arias. (2011). Anthropogenic Sulfur Dioxide Emissions: 1850-2005, Atmospheric Chemistry and Physics, 11:1101-1116.
}
\examples{
data(so2data)

##create annual zoo object:
newso2data<- zooreg(so2data[,-1], start=1946, frequency=1)

##plot UK annual total anthropogenic SO2 emissions:
plot(newso2data$uk_tot_so2)

}
\keyword{datasets}
\keyword{emissions}
