\name{gets.isat}
\alias{gets.isat}

\title{General-to-Specific (GETS) Modelling 'isat' objects}

\description{
General-to-Specific (GETS) Modelling of a objects of class \code{\link{isat}}.
}

\usage{
\method{gets}{isat}(x, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
  gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
  include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
  turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE,...) 
}

\arguments{
  \item{x}{an object of class 'isat'}
  \item{t.pval}{numeric value between 0 and 1. The significance level used for the two-sided regressor significance t-tests}
  \item{wald.pval}{numeric value between 0 and 1. The significance level used for the Parsimonious Encompassing Tests (PETs). By default, it is the same as \code{t.pval}}
  \item{vcov.type}{the type of variance-covariance matrix used. If \code{NULL} (default), then the type used in the estimation of the 'arx' object is used. This can be overridden by either "ordinary" (i.e. the ordinary variance-covariance matrix) or "white" (i.e. the White (1980) heteroscedasticity robust variance-covariance matrix)}
  \item{do.pet}{logical. If \code{TRUE} (default), then a Parsimonious Encompassing Test (PET) against the GUM is undertaken at each regressor removal for the joint significance of all the deleted regressors along the current path. If \code{FALSE}, then a PET is not undertaken at each regressor removal}
  \item{ar.LjungB}{a two-item list with names \code{lag} and \code{pval}, or \code{NULL}. In the former case \code{lag} contains the order of the Ljung and Box (1979) test for serial correlation in the standardised residuals, and \code{pval} contains the significance level. If \code{lag=NULL} (default), then the order used is that of the estimated 'arx' object. If \code{ar.Ljungb=NULL}, then the standardised residuals are not checked for serial correlation}
  \item{arch.LjungB}{a two-item list with names \code{lag} and \code{pval}, or \code{NULL}. In the former case, \code{lag} contains the order of the Ljung and Box (1979) test for serial correlation in the squared standardised residuals, and \code{pval} contains the significance level. If \code{lag=NULL} (default), then the order used is that of the estimated 'arx' object. If \code{arch.Ljungb=NULL}, then the standardised residuals are not checked for ARCH}
  \item{normality.JarqueB}{a value between 0 and 1, or \code{NULL}. In the former case, the Jarque and Bera (1980) test for non-normality is conducted using a significance level equal to the numeric value. If \code{NULL}, then no test for non-normality is undertaken}
  \item{user.diagnostics}{\code{NULL} or a \code{\link{list}} with two entries, \code{name} and \code{pval}, see the \code{user.fun} argument in \code{\link{diagnostics}} }
  \item{info.method}{character string, "sc" (default), "aic" or "hq", which determines the information criterion to be used when selecting among terminal models. The abbreviations are short for the Schwarz or Bayesian information criterion (sc), the Akaike information criterion (aic) and the Hannan-Quinn (hq) information criterion}
  \item{gof.function}{\code{NULL} (default) or a \code{list}, see \code{\link{getsFun}}. If \code{NULL}, then \code{\link{infocrit}} is used }
  \item{gof.method}{\code{NULL} (default) or a \code{character}, see \code{\link{getsFun}}. If \code{NULL} and \code{gof.function} is also \code{NULL}, then the best goodness-of-fit is characterised by a minimum value}
  \item{keep}{the regressors to be excluded from removal in the specification search. Note that \code{keep=c(1)} is obligatory when using \code{getsv}. This excludes the log-variance intercept from removal. The regressor numbering is contained in the \code{reg.no} column of the GUM}
  \item{include.gum}{logical. If \code{TRUE}, then the GUM (i.e. the starting model) is included among the terminal models. If \code{FALSE} (default), then the GUM is not included}
  \item{include.1cut}{logical. If \code{TRUE}, then the 1-cut model is added to the list of terminal models. If \code{FALSE} (default), then the 1-cut is not added, unless it is a terminal model in one of the paths }
  \item{include.empty}{logical. If \code{TRUE}, then an empty model is included among the terminal models, if it passes the diagnostic tests, even if it is not equal to one of the terminals. If \code{FALSE} (default), then the empty model is not included (unless it is one of the terminals)}
  \item{max.paths}{\code{NULL} (default) or an integer greater than 0. If \code{NULL}, then there is no limit to the number of paths. If an integer (e.g. 1), then this integer constitutes the maximum number of paths searched (e.g. a single path) }
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the variance-covariance matrix when computing the Wald-statistic used in the Parsimonious Encompassing Tests (PETs), see the \code{\link{qr.solve}} function}
  \item{turbo}{logical. If \code{TRUE}, then (parts of) paths are not searched twice (or more) unnecessarily, thus yielding a significant potential for speed-gain. However, the checking of whether the search has arrived at a point it has already been comes with a slight computational overhead. Accordingly, if \code{turbo=TRUE}, then the total search time might in fact be higher than if \code{turbo=FALSE}. This happens if estimation is very fast, say, less than quarter of a second. Hence the default is \code{FALSE}}
  \item{print.searchinfo}{logical. If \code{TRUE} (default), then a print is returned whenever simiplification along a new path is started}
  \item{plot}{\code{NULL} or logical. If \code{TRUE}, then the fitted values and the residuals of the final model are plotted after model selection. If \code{FALSE}, then they are not. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not}
  \item{alarm}{logical. If \code{TRUE}, then a sound or beep is emitted (in order to alert the user) when the model selection ends}
  \item{\dots}{further arguments passed on to and from methods }
}


\details{
Internally, \code{gets.isat} invokes \code{\link{getsm}} for the GETS-modelling.
}

\value{A list of class \code{\link{gets}}.  }

\author{
Moritz Schwarz, \url{https://www.inet.ox.ac.uk/people/moritz-schwarz}\cr
Genaro Sucarrat, \url{https://www.sucarrat.net/}
}

\seealso{
\code{\link{isat}}, \code{\link{getsm}}, \code{\link{getsFun}}, \code{\link{paths}} and \code{\link{terminals}}
}

\examples{

##generate some data:
#set.seed(123) #for reproducibility
#y <- rnorm(30) #generate Y
#isatmod <- isat(y)
#gets(isatmot)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
