\name{blocksFun}
\alias{blocksFun}

\title{ Block-based General-to-Specific (GETS) modelling }

\description{
Auxiliary function (i.e. not intended for the average user) that enables block-based GETS-modelling with user-specified estimators, diagnostics and goodness-of-fit criteria.
}

\usage{
blocksFun(y, x, untransformed.residuals=NULL, blocks=NULL,
  no.of.blocks=NULL, max.block.size=30, ratio.threshold=0.8,
  gets.of.union=TRUE, force.invertibility=FALSE,
  user.estimator=list(name="ols"), t.pval=0.001, wald.pval=t.pval,
  do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL, normality.JarqueB=NULL,
  user.diagnostics=NULL, gof.function=list(name="infocrit"),
  gof.method=c("min", "max"), keep=NULL, include.gum=FALSE,
  include.1cut=FALSE, include.empty=FALSE, max.paths=NULL,
  turbo=FALSE, parallel.options=NULL, tol=1e-07, LAPACK=FALSE,
  max.regs=NULL, print.searchinfo=TRUE, alarm=FALSE)
}

\arguments{
  \item{y}{a numeric vector (with no missing values, i.e. no non-numeric 'holes') }
  \item{x}{a \code{matrix}, or a \code{list} of matrices }
  \item{untransformed.residuals}{\code{NULL} (default) or, when \code{\link{ols}} is used with \code{method=6} in \code{user.estimator}, a numeric vector containing the untransformed residuals }
  \item{blocks}{\code{NULL} (default) or a \code{list} of lists with vectors of integers that indicate how blocks should be put together. If \code{NULL}, then the block composition is undertaken automatically by an internal algorithm that depends on \code{no.of.blocks}, \code{max.block.size} and \code{ratio.threshold} }
  \item{no.of.blocks}{\code{NULL} (default) or \code{integer}. If \code{NULL}, then the number of blocks is determined automatically by an internal algorithm }
  \item{max.block.size}{\code{integer} that controls the size of blocks }
  \item{ratio.threshold}{\code{numeric} between 0 and 1 that controls the minimum ratio of variables in each block to total observations }
  \item{gets.of.union}{\code{logical}. If \code{TRUE}, then GETS modelling is undertaken of the union of retained variables. Otherwise it is not }
  \item{force.invertibility}{\code{logical}. If \code{TRUE}, then the x-matrix is ensured to have full row-rank before it is passed on to \code{\link{getsFun}} }
  \item{user.estimator}{\code{list}, see \code{\link{getsFun}} for the details }
  \item{t.pval}{\code{numeric} value between 0 and 1. The significance level used for the two-sided coefficient significance t-tests }
  \item{wald.pval}{\code{numeric} value between 0 and 1. The significance level used for the Parsimonious Encompassing Tests (PETs) }
  \item{do.pet}{\code{logical}. If \code{TRUE}, then a Parsimonious Encompassing Test (PET) against the GUM is undertaken at each variable removal for the joint significance of all the deleted regressors along the current GETS path. If \code{FALSE}, then a PET is not undertaken at each removal}
  \item{ar.LjungB}{a two element \code{vector}, or \code{NULL}. In the former case, the first element contains the AR-order, the second element the significance level. If \code{NULL}, then a test for autocorrelation in the residuals is not conducted }
  \item{arch.LjungB}{a two element \code{vector}, or \code{NULL}. In the former case, the first element contains the ARCH-order, the second element the significance level. If \code{NULL}, then a test for ARCH in the residuals is not conducted }
  \item{normality.JarqueB}{\code{NULL} or a \code{numeric} value between 0 and 1. In the latter case, a test for non-normality in the residuals is conducted using a significance level equal to\cr \code{normality.JarqueB}. If \code{NULL}, then no test for non-normality is conducted }
  \item{user.diagnostics}{\code{NULL} (default) or a \code{list} with two entries, \code{name} and \code{pval}. See \code{\link{getsFun}} for the details }
  \item{gof.function}{\code{list}. The first item should be named \code{name} and contain the name (a character) of the Goodness-of-Fit (GOF) function used. Additional items in the list \code{gof.function} are passed on as arguments to the GOF-function. . See \code{\link{getsFun}} for the details }
  \item{gof.method}{\code{character}. Determines whether the best Goodness-of-Fit is a minimum (default) or maximum }
  \item{keep}{\code{NULL} (default), \code{vector} of integers or a \code{list} of vectors of integers. In the latter case, the number of vectors should be equal to the number of matrices in \code{x} }
  \item{include.gum}{\code{logical}. If \code{TRUE}, then the GUM (i.e. the starting model) is included among the terminal models }
  \item{include.1cut}{\code{logical}. If \code{TRUE}, then the 1-cut model is added to the list of terminal models }
  \item{include.empty}{\code{logical}. If \code{TRUE}, then the empty model is added to the list of terminal models }
  \item{max.paths}{\code{NULL} (default) or \code{integer} greater than 0. If \code{NULL}, then there is no limit to the number of paths. If \code{integer} (e.g. 1), then this integer constitutes the maximum number of paths searched (e.g. a single path) }
  \item{turbo}{\code{logical}. If \code{TRUE}, then (parts of) paths are not searched twice (or more) unnecessarily in each GETS modelling. Setting \code{turbo} to \code{TRUE} entails a small additional computational costs, but may be outweighed substantially if estimation is slow, or if the number of variables to delete in each path is large }
  \item{parallel.options}{\code{NULL} or \code{integer} that indicates the number of cores/threads to use for parallel computing (implemented w/\code{makeCluster} and \code{parLapply})}
  \item{tol}{\code{numeric} value, the tolerance for detecting linear dependencies in the columns of the variance-covariance matrix when computing the Wald-statistic used in the Parsimonious Encompassing Tests (PETs), see the \code{\link{qr.solve}} function}
  \item{LAPACK}{currently not used}
  \item{max.regs}{\code{integer}. The maximum number of regressions along a deletion path. Do not alter unless you know what you are doing! }
  \item{print.searchinfo}{\code{logical}. If \code{TRUE} (default), then a print is returned whenever simiplification along a new path is started }
  \item{alarm}{\code{logical}. If \code{TRUE}, then a sound or beep is emitted (in order to alert the user) when the model selection ends }

}

\details{
\code{blocksFun} undertakes block-based GETS modelling by a repeated but structured call to \code{getsFun}. For the details of how to user-specify an estimator via \code{user.estimator}, diagnostics via
\code{user.diagnostics} and a goodness-of-fit function via \code{gof.function}, see documentation of \code{\link{getsFun}} under "Details".

The algorithm of \code{blocksFun} is similar to that of \code{\link{isat}}, but more flexible. The main use of \code{blocksFun} is the creation of user-specified methods that employs block-based GETS modelling, e.g. indicator saturation techniques.
}

\value{ A \code{\link{list}} with the results of the block-based GETS-modelling. }

\references{
F. Pretis, J. Reade and G. Sucarrat (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44

G. sucarrat (2019): 'User-Specified General-to-Specific and Indicator Saturation Methods', Munich Personal RePEc Archive: \url{https://mpra.ub.uni-muenchen.de/96653/}
}

\author{
Genaro Sucarrat, with contributions from Jonas kurle, Felix Pretis and James Reade \cr
}

\seealso{
\code{\link{ols}}, \code{\link{diagnostics}}, \code{\link{infocrit}} and \code{\link{isat}}
}

\examples{

## more variables than observations:
y <- rnorm(20)
x <- matrix(rnorm(length(y)*40), length(y), 40)
blocksFun(y, x)

## 'x' as list of matrices:
z <- matrix(rnorm(length(y)*40), length(y), 40)
blocksFun(y, list(x,z))

## ensure regressor no. 3 in matrix no. 2 is not removed:
blocksFun(y, list(x,z), keep=list(integer(0), 3))

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
