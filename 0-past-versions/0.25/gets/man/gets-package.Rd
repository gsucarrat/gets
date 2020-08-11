\name{gets-package}
\alias{gets-package}
\docType{package}
\title{General-to-Specific Modelling and Indicator Saturation Methods}

\description{
Multi-path General-to-Specific (GETS) modelling of the mean and/or variance of a regression, and Indicator Saturation (ISAT) methods for detecting structural breaks in the mean. The mean can be specified as an autoregressive model with covariates (an 'AR-X' model), and the variance can be specified as a dynamic log-variance model with covariates (a 'log-ARCH-X' model). For the statistical details of the model, see Section 4 in Pretis, Reade and Sucarrat (2018).

The main functions of the package are \code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}} and \code{\link{isat}}. The first function, \code{arx}, estimates an AR-X model with (optionally) a log-ARCH-X specification on the log-variance. The second function, \code{getsm}, undertakes GETS model selection of the mean specification of an \code{arx} object. The third function, \code{getsv}, undertakes GETS model selection of the log-variance specification of an \code{arx} object. The fourth function, \code{isat}, undertakes GETS model selection of an indicator saturated mean specification. Extraction functions (mainly S3 methods) are also available, together with additional auxiliary functions used by the main functions.

For an introduction to the package, see Pretis, Reade and Sucarrat (2018): \url{https://www.jstatsoft.org/article/view/v086i03}. The package also provides facilities for user-defined GETS and ISAT methods. While this is to some extent available via the main functions, full flexibility is provided by \code{\link{getsFun}} and \code{\link{blocksFun}}, see Sucarrat (2019): \url{https://mpra.ub.uni-muenchen.de/96653/}.

}

\details{

\tabular{ll}{
Package: \tab gets\cr
Type: \tab Package\cr
Version: \tab 0.25\cr
Date: \tab 2020-08-11\cr
License: \tab GPL-2\cr
}
The code originated in relation with G. Sucarrat and A. Escribano (2012). Felix Pretis and James Reade joined for the development of the \code{\link{isat}} code and related functions. Subsequently, Moritz Schwarz and Jonas Kurle have made various contributions.

}

\author{
Jonas Kurle, \url{https://www.jonaskurle.com/}\cr
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade,  \url{https://sites.google.com/site/jjamesreade/}\cr
Moritz Schwarz, \url{https://www.inet.ox.ac.uk/people/moritz-schwarz/}\cr
Genaro Sucarrat \url{http://www.sucarrat.net/}\cr\cr
Maintainer: Genaro Sucarrat            
}

\references{
Jurgen A. Doornik, David F. Hendry, and Felix Pretis (2013): 'Step Indicator Saturation', Oxford Economics Discussion Paper, 658.

Felix Pretis, James Reade and Genaro Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. \url{https://www.jstatsoft.org/article/view/v086i03}

Carlos Santos, David F. Hendry and Soren Johansen (2007):  'Automatic selection of indicators in a fully saturated regression'. Computational Statistics, vol 23:1, pp.317-335

Genaro Sucarrat (2019): 'User-Specified General-to-Specific and Indicator Saturation Methods'. \url{https://mpra.ub.uni-muenchen.de/96653/}

Genaro Sucarrat and Alvaro Escribano (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735. \url{https://doi.org/10.1111/j.1468-0084.2011.00669.x}

}

\seealso{
\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}, \code{\link{getsFun}}, \code{\link{blocksFun}}
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
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
