\name{gets-package}
\alias{gets-package}
\docType{package}
\title{General-to-Specific (GETS) Modelling and Indicator Saturation Methods}

\description{
Automated multi-path General-to-Specific (GETS) modelling of the mean and variance of a regression, and indicator saturation methods for detecting structural breaks in the mean. The mean can be specified as an autoregressive model with covariates (an 'AR-X' model), and the variance can be specified as a dynamic log-variance model with covariates (a 'log-ARCH-X' model). For the statistical details of the model, see Section 4 in Pretis, Reade and Sucarrat (2018).

The four main functions of the package are \code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}} and \code{\link{isat}}. The first function, \code{arx}, estimates an AR-X model with (optionally) a log-ARCH-X specification on the log-variance. The second function, \code{getsm}, undertakes GETS model selection of the mean specification of an \code{arx} object. The third function, \code{getsv}, undertakes GETS model selection of the log-variance specification of an \code{arx} object. The fourth function, \code{isat}, undertakes GETS model selection of an indicator saturated mean specification. Extraction functions (mainly S3 methods) are also available, together with additional auxiliary functions used by the main functions.

For the advanced user, the package provides facilities for user-defined GETS and ISAT methods. To this end, the main functions are \code{\link{getsFun}} and \code{\link{isat}}.
}

\details{
\tabular{ll}{
Package: \tab gets\cr
Type: \tab Package\cr
Version: \tab 0.20\cr
Date: \tab 2019-09-18\cr
License: \tab GPL-2\cr
}
The code originated in relation with G. Sucarrat and A. Escribano (2012). Subsequently, Felix Pretis and James Reade joined for the development of the \code{\link{isat}} code and related functions. For a recent introduction to the package, see Pretis, Reade and Sucarrat (2018).
}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade,  \url{https://sites.google.com/site/jjamesreade/}\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}\cr\cr
Maintainer: Genaro Sucarrat
}

\references{
Jurgen, A. Doornik, Hendry, David F., and Pretis, Felix (2013): 'Step Indicator Saturation', Oxford Economics Discussion Paper, 658.

Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. DOI: https://www.jstatsoft.org/article/view/v086i03

Carlos Santos, Hendry, David, F. and Johansen, Soren (2007):  'Automatic selection of indicators in a fully saturated regression'. Computational Statistics, vol 23:1, pp.317-335

Sucarrat, Genaro and Escribano, Alvaro (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735
}

\seealso{
\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}
}

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 60)

##Estimate an AR(2) with intercept as mean specification
##and a log-ARCH(4) as log-volatility specification:
myModel <- arx(y, mc=TRUE, ar=1:2, arch=1:4)

##GETS modelling of the mean of myModel:
simpleMean <- getsm(myModel)

##GETS modelling of the log-variance of myModel:
simpleVar <- getsv(myModel)

##results:
print(simpleMean)
print(simpleVar)

##step indicator saturation of an iid normal series:
set.seed(123)
y <- rnorm(30)
isat(y)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
