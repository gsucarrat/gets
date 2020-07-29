\name{getsm}
\alias{getsm}
\alias{getsv}

\title{General-to-Specific (GETS) Modelling of an AR-X model (the mean specification) with log-ARCH-X errors (the log-variance specification).}

\description{The starting model, an object of the 'arx' class, is referred to as the General Unrestricted Model (GUM). The \code{getsm} function undertakes multi-path GETS modelling of the mean specification, whereas \code{getsv} does the same for the log-variance specification. The diagnostic tests are undertaken on the standardised residuals, and the \code{keep} option enables regressors to be excluded from possible removal.}

\usage{
##GETS-modelling of mean specification:
getsm(object, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL, 
    do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025), 
    arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL, 
    user.diagnostics=NULL, info.method=c("sc","aic","aicc", "hq"),
    gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
    include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, turbo=FALSE,
    print.searchinfo=TRUE, plot=NULL, alarm=FALSE)

##GETS modelling of log-variance specification:
getsv(object, t.pval=0.05, wald.pval=t.pval,
    do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
    arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
    user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
    gof.function=NULL, gof.method=NULL, keep=c(1), include.gum=FALSE,
    include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, turbo=FALSE,
    print.searchinfo=TRUE, plot=NULL, alarm=FALSE)
}

\arguments{
  \item{object}{an object of class 'arx'}
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
  \item{turbo}{logical. If \code{TRUE}, then (parts of) paths are not searched twice (or more) unnecessarily, thus yielding a significant potential for speed-gain. However, the checking of whether the search has arrived at a point it has already been comes with a slight computational overhead. Accordingly, if \code{turbo=TRUE}, then the total search time might in fact be higher than if \code{turbo=FALSE}. This happens if estimation is very fast, say, less than quarter of a second. Hence the default is \code{FALSE}}
  \item{print.searchinfo}{logical. If \code{TRUE} (default), then a print is returned whenever simiplification along a new path is started}
  \item{plot}{\code{NULL} or logical. If \code{TRUE}, then the fitted values and the residuals of the final model are plotted after model selection. If \code{FALSE}, then they are not. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not}
  \item{alarm}{logical. If \code{TRUE}, then a sound or beep is emitted (in order to alert the user) when the model selection ends}
}

\details{
For an overview, see Pretis, Reade and Sucarrat (2018): \url{https://www.jstatsoft.org/article/view/v086i03}.

The arguments \code{user.diagnostics} and \code{gof.function} enable the specification of user-defined diagnostics and a user-defined goodness-of-fit function. For the former, see the documentation of \code{\link{diagnostics}}. For the latter, the principles of the same arguments in \code{\link{getsFun}} are followed, see its documentation under "Details", and Sucarrat (2019): \url{https://mpra.ub.uni-muenchen.de/96653/}.
}

\value{A list of class 'gets'}

\references{
C. Jarque and A. Bera (1980): 'Efficient Tests for Normality, Homoscedasticity and Serial Independence'. Economics Letters 6, pp. 255-259

G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time Series Models'. Biometrika 66, pp. 265-270

Felix Pretis, James Reade and Genaro Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. \url{https://www.jstatsoft.org/article/view/v086i03}

Genaro Sucarrat (2019): 'User-Specified General-to-Specific and Indicator Saturation Methods'. \url{https://mpra.ub.uni-muenchen.de/96653/}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{
Extraction functions: \code{\link{coef.gets}}, \code{\link{fitted.gets}}, \code{\link{paths}}, \code{\link{plot.gets}}, \code{\link{print.gets}},\cr
\code{\link{residuals.gets}}, \code{\link{summary.gets}}, \code{\link{terminals}}, \code{\link{vcov.gets}}\cr

Related functions: \code{\link{arx}}, \code{\link{eqwma}}, \code{\link{leqwma}}, \code{\link{zoo}}, \code{\link{getsFun}}
}

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 80)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(2*80), 80, 2)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and a log-ARCH(3) with log(xregs^2) as
##regressors in the log-variance:
gum01 <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3,
  vxreg=log(xregs^2))

##GETS model selection of the mean:
meanmod01 <- getsm(gum01)

##GETS model selection of the log-variance:
varmod01 <- getsv(gum01)

##GETS model selection of the mean with the mean intercept
##excluded from removal:
meanmod02 <- getsm(gum01, keep=1)

##GETS model selection of the mean with non-default
#serial-correlation diagnostics settings:
meanmod03 <- getsm(gum01, ar.LjungB=list(pval=0.05))

##GETS model selection of the mean with very liberal
##(20 percent) significance levels:
meanmod04 <- getsm(gum01, t.pval=0.2)

##GETS model selection of log-variance with all the
##log-ARCH terms excluded from removal:
varmod03 <- getsv(gum01, keep=2:4)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
