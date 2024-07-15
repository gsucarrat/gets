\name{gets.larch}
\alias{gets.larch}

\title{General-to-Specific (GETS) Modelling of a heterogeneous log-ARCH-X model}

\description{The starting model, an object of the 'larch' class (see \code{\link{larch}}, is referred to as the General Unrestricted Model (GUM). The \code{gets.larch()} function undertakes multi-path GETS modelling of the log-variance specification. The diagnostic tests are undertaken on the standardised residuals, and the \code{keep} option enables regressors to be excluded from possible removal. }

\usage{
  \method{gets}{larch}(x, t.pval=0.05, wald.pval=t.pval, do.pet=TRUE, 
    ar.LjungB=NULL, arch.LjungB=NULL, normality.JarqueB=NULL, 
    user.diagnostics=NULL, info.method=c("sc", "aic", "aicc", "hq"),
    gof.function=NULL, gof.method=NULL, keep=c(1), include.gum=FALSE,
    include.1cut=TRUE, include.empty=FALSE, max.paths=NULL, tol=1e-07,
    turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE, ...)
}

\arguments{
  \item{x}{an object of class 'larch'}
  \item{t.pval}{numeric value between 0 and 1. The significance level used for the two-sided regressor significance t-tests}
  \item{wald.pval}{numeric value between 0 and 1. The significance level used for the Parsimonious Encompassing Tests (PETs). By default, \code{wald.pval} is equal to \code{t.pval}}
  \item{do.pet}{logical. If \code{TRUE} (default), then a Parsimonious Encompassing Test (PET) against the GUM is undertaken at each regressor removal for the joint significance of all the deleted regressors along the current path. If \code{FALSE}, then a PET is not undertaken at each regressor removal}
  \item{ar.LjungB}{\code{NULL} (default), or a \code{\link{list}} with named items \code{lag} and \code{pval}, or a two-element numeric vector where the first element contains the lag and the second the p-value. If \code{NULL}, then the standardised residuals are not checked for autocorrelation. If \code{ar.LjungB} is a \code{list}, then \code{lag} contains the order of the Ljung and Box (1979) test for serial correlation in the standardised residuals, and \code{pval} contains the significance level. If \code{lag=NULL}, then the order used is that of the estimated 'larch' object }
  \item{arch.LjungB}{\code{NULL} (default), or a \code{\link{list}} with named items \code{lag} and \code{pval}, or a two-element numeric vector where the first element contains the lag and the second the p-value. If \code{NULL}, then the standardised residuals are not checked for ARCH (autocorrelation in the squared standardised residuals). If \code{ar.LjungB} is a \code{list}, then \code{lag} contains the order of the test, and \code{pval} contains the significance level. If \code{lag=NULL}, then the order used is that of the estimated 'larch' object }
  \item{normality.JarqueB}{\code{NULL} (default) or a numeric value between 0 and 1. If \code{NULL}, then no test for non-normality is undertaken. If a numeric value between 0 and 1, then the Jarque and Bera (1980) test for non-normality is conducted using a significance level equal to the numeric value }
  \item{user.diagnostics}{\code{NULL} (default) or a \code{\link{list}} with two entries, \code{name} and \code{pval}, see the \code{user.fun} argument in \code{\link{diagnostics}} }
  \item{info.method}{character string, "sc" (default), "aic", "aicc" or "hq", which determines the information criterion to be used when selecting among terminal models. See \code{\link{infocrit}} for the details }
  \item{gof.function}{\code{NULL} (default) or a \code{list}, see \code{\link{getsFun}}. If \code{NULL}, then \code{\link{infocrit}} is used }
  \item{gof.method}{\code{NULL} (default) or a \code{character}, see \code{\link{getsFun}}. If \code{NULL} and \code{gof.function} is also \code{NULL}, then the best goodness-of-fit is characterised by a minimum value }
  \item{keep}{the regressors to be kept (i.e. excluded from removal) in the specification search. Currently, \code{keep=c(1)} is obligatory, which excludes the log-variance intercept from removal }
  \item{include.gum}{logical. If \code{TRUE}, the GUM (i.e. the starting model) is included among the terminal models. If \code{FALSE} (default), the GUM is not included }
  \item{include.1cut}{logical. If \code{TRUE} (default), then the 1-cut model is added to the list of terminal models. If \code{FALSE}, the 1-cut is not added, unless it is a terminal model in one of the paths }
  \item{include.empty}{logical. If \code{TRUE}, then an empty model is included among the terminal models, if it passes the diagnostic tests. If \code{FALSE} (default), then the empty model is not included }
  \item{max.paths}{\code{NULL} (default) or an integer equal to or greater than 0. If \code{NULL}, then there is no limit to the number of paths. If an integer (e.g. 1), then this integer constitutes the maximum number of paths searched (e.g. a single path) }
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the variance-covariance matrix when computing the Wald-statistic used in the Parsimonious Encompassing Tests (PETs), see the \code{\link{qr.solve}} function}
  \item{turbo}{logical. If \code{TRUE}, then paths are not searched twice (or more) unnecessarily, thus yielding a significant potential for speed-gain. However, the checking of whether the search has arrived at a point it has already been comes with a computational overhead. Accordingly, if \code{turbo=TRUE}, the total search time might in fact be higher than if \code{turbo=FALSE}. This is particularly likely to happen if estimation is very fast, say, less than a quarter of a second. Hence the default is \code{FALSE} }
  \item{print.searchinfo}{logical. If \code{TRUE} (default), then a print is returned whenever simiplification along a new path is started}
  \item{plot}{\code{NULL} or logical. If \code{TRUE}, then the fitted values and the standardised residuals of the final model are plotted after model selection. If \code{FALSE}, then they are not plotted. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not}
  \item{alarm}{logical. If \code{TRUE}, then a sound or beep is emitted (in order to alert the user) when the model selection ends, see \code{\link{alarm}} }
  \item{\dots}{additional arguments}
}

\details{

See Pretis, Reade and Sucarrat (2018): \doi{10.18637/jss.v086.i03}, and Sucarrat (2020): \url{https://journal.r-project.org/archive/2021/RJ-2021-024/}.\cr

The arguments \code{user.diagnostics} and \code{gof.function} enable the specification of user-defined diagnostics and a user-defined goodness-of-fit function. For the former, see the documentation of \code{\link{diagnostics}}. For the latter, the principles of the same arguments in \code{\link{getsFun}} are followed, see its documentation under "Details", and Sucarrat (2020): \url{https://journal.r-project.org/archive/2021/RJ-2021-024/}.

}

\value{A list of class 'larch', see \code{\link{larch}}, with additional information about the GETS modelling }

\references{

C. Jarque and A. Bera (1980): 'Efficient Tests for Normality, Homoscedasticity and Serial Independence'. Economics Letters 6, pp. 255-259. \doi{10.1016/0165-1765(80)90024-5}

G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time Series Models'. Biometrika 66, pp. 265-270

Felix Pretis, James Reade and Genaro Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44. \doi{10.18637/jss.v086.i03}

Genaro Sucarrat (2020): 'User-Specified General-to-Specific and Indicator Saturation Methods'. The R Journal 12:2, pages 388-401. \url{https://journal.r-project.org/archive/2021/RJ-2021-024/}

}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{
Methods and extraction functions (mostly S3 methods): \code{\link{coef.larch}}, \code{\link{ES}}, \code{\link{fitted.larch}}, \code{\link{gets.larch}}, \cr
\code{\link{logLik.larch}}, \code{\link{nobs.larch}}, \code{\link{plot.larch}}, \code{\link{predict.larch}}, \code{\link{print.larch}}, \cr
\code{\link{residuals.larch}}, \code{\link{summary.larch}}, \code{\link{VaR}}, \code{\link{toLatex.larch}} and \code{\link{vcov.arx}}\cr

Related functions: \code{\link{eqwma}}, \code{\link{leqwma}}, \code{\link{regressorsVariance}}, \code{\link[zoo]{zoo}}, \code{\link{getsFun}}, \code{\link{qr.solve}}
}

\examples{
##Simulate some data:
set.seed(123)
e <- rnorm(40)
x <- matrix(rnorm(4*40), 40, 4)

##estimate a log-ARCH(3) with asymmetry and log(x^2) as regressors:
gum <- larch(e, arch=1:3, asym=1, vxreg=log(x^2))

##GETS modelling of the log-variance:
simple <- gets(gum)

##GETS modelling with intercept and log-ARCH(1) terms
##excluded from removal:
simple <- gets(gum, keep=c(1,2))

##GETS modelling with non-default autocorrelation
##diagnostics settings:
simple <- gets(gum, ar.LjungB=list(pval=0.05))

##GETS modelling with very liberal (40%) significance level:
simple <- gets(gum, t.pval=0.4)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
