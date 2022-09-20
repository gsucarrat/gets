\name{gets-package}
\alias{gets-package}
\docType{package}
\title{General-to-Specific (GETS) and Indicator Saturation (ISAT) Modelling}

\description{

The \pkg{gets} package provides functions and methods for General-to-Specific (GETS) and Indicator Saturation (ISAT) modelling. GETS modelling is a powerful and flexible variable selection algorithm that returns a parsimonious and interpretable model. It is ideally suited for the development of models that can be used for counterfactual and predictive scenario analysis (e.g. conditional forecasting). ISAT modelling provides a comprehensive, flexible and powerful approach to the identification of structural breaks and outliers.\cr

The code of the package originated in relation with the research project G. Sucarrat and A. Escribano (2012). In 2014, Felix Pretis and James Reade joined for the development of the \code{\link{isat}} code and related functions. Moritz Schwarz and Jonas Kurle joined the development team in 2020.

}

\details{

\tabular{ll}{
	Version: \tab 0.37\cr
	Date: \tab 2022-09-20\cr
	Licence: \tab GPL-2\cr
}

}

\section{GETS modelling}{

	GETS methods are available for the following model classes:\cr

	\itemize{
	  \item Linear regression, both static and dynamic, see \code{\link{arx}}, \code{\link{gets.arx}} and \code{\link{gets.lm}}
	  
	  \item Variance models, both static and dynamic, see \code{\link{arx}}
	  
	  \item Logit models, both static and dynamic, see \code{\link{logitx}} and \code{\link{gets.logitx}}
	}

	The function \code{arx} estimates a static linear regression, or a dynamic AR-X model with (optionally) a log-variance specification. The log-variance specification can either be static or a dynamic log-variance model with covariates (a 'log-ARCH-X' model). For the statistical details of the model, see Section 4 in Pretis, Reade and Sucarrat (2018). The function \code{\link{logitx}} estimates a static logit model, or a dynamic logit model with covariates (optionally). For complete user-specified GETS modelling, see \code{\link{getsFun}}.\cr

}

\section{ISAT modelling}{

	ISAT methods are available for:\cr

	\itemize{
	  \item Linear regression, both static and dynamic, see \code{\link{isat}}
	}

	The \code{isat} function undertakes GETS model selection of an indicator saturated mean specification. Extraction functions (mainly S3 methods) are also available, together with additional auxiliary functions. For complete user-specified ISAT modelling, see \code{\link{blocksFun}}.\cr

}

\section{Vignettes}{

	Two vignettes are available in the package (type \code{browseVignettes("gets")} to access them):\cr

	\itemize{
	  \item An introduction to the \emph{gets} package
	  \item User-Specified General-to-Specific (GETS) and Indicator Saturation (ISAT) Methods
	}

	The former is a mildly modified version of Pretis, Reade and Sucarrat (2018), whereas the latter is an updated version of Sucarrat (2020).

}

\author{

\tabular{lll}{
	Jonas Kurle: \tab \tab \url{https://www.jonaskurle.com/}\cr
	Felix Pretis: \tab \tab \url{http://www.felixpretis.org/}\cr
	James Reade: \tab \tab \url{https://sites.google.com/site/jjamesreade/}\cr
	Moritz Schwarz: \tab \tab  \url{https://www.inet.ox.ac.uk/people/moritz-schwarz/}\cr
	Genaro Sucarrat: \tab \tab \url{http://www.sucarrat.net/}\cr
}

Maintainer: Genaro Sucarrat            

}

\references{
Jurgen A. Doornik, David F. Hendry, and Felix Pretis (2013): 'Step Indicator Saturation', Oxford Economics Discussion Paper, 658. \url{https://ideas.repec.org/p/oxf/wpaper/658.html}

Felix Pretis, James Reade and Genaro Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. \doi{10.18637/jss.v086.i03}

Carlos Santos, David F. Hendry and Soren Johansen (2007):  'Automatic selection of indicators in a fully saturated regression'. Computational Statistics, vol 23:1, pp.317-335. \doi{10.1007/s00180-007-0054-z}

Genaro Sucarrat (2020): 'User-Specified General-to-Specific and Indicator Saturation Methods'. The R Journal 12:2, pages 388-401. \url{https://journal.r-project.org/archive/2021/RJ-2021-024/}

Genaro Sucarrat and Alvaro Escribano (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735.

}

\seealso{
\code{\link{arx}}, \code{\link{gets.arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}, \code{\link{getsFun}}, \code{\link{blocksFun}}
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
