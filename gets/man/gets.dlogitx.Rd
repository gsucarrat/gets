\name{gets.dlogitx}
\alias{gets.dlogitx}

\title{General-to-Specific (GETS) Modelling of objects of class 'dlogitx'}

\description{
General-to-Specific (GETS) Modelling of autoregressive logit models of class 'dlogitx'
}

\usage{
\method{gets}{dlogitx}(x, t.pval = 0.05, wald.pval = t.pval, do.pet = TRUE, 
    keep = NULL, include.gum = FALSE, include.1cut = TRUE, include.empty = FALSE, 
    max.paths = NULL, turbo = TRUE, print.searchinfo = TRUE, 
    plot = NULL, alarm = FALSE, ...) 
}

\arguments{
  \item{x}{an object of class 'dlogitx', see \code{\link{dlogitx}} }
  \item{t.pval}{numeric value between 0 and 1. The significance level used for the two-sided regressor significance t-tests}
  \item{wald.pval}{numeric value between 0 and 1. The significance level used for the Parsimonious Encompassing Tests (PETs). By default, it is the same as \code{t.pval}}
  \item{do.pet}{\code{logical} that determines whether a Parsimonious Encompassing Test (PET) against the GUM should be undertaken at each regressor removal for the joint significance of all the deleted regressors along the current path. If \code{FALSE}, then a PET is not undertaken at each regressor removal }
  \item{keep}{\code{NULL} or a vector of integers that determines which regressors to be excluded from removal in the specification search }
  \item{include.gum}{\code{logical} that determines whether the GUM (i.e. the starting model) should be included among the terminal models. If \code{FALSE} (default), then the GUM is not included }
  \item{include.1cut}{\code{logical} that determines whether the 1-cut model should be added to the list of terminal models. If \code{FALSE} (default), then the 1-cut is not added, unless it is a terminal model in one of the paths }
  \item{include.empty}{\code{logical} that determines whether an empty model should be added to the list of terminal models, if it passes the diagnostic tests. If \code{FALSE} (default), then the empty model is not added, unless it is a terminal model in one of the paths }
  \item{max.paths}{\code{NULL} (default) or an integer greater than 0. If \code{NULL}, then there is no limit to the number of paths. If an integer (e.g. 1), then this integer constitutes the maximum number of paths searched (e.g. a single path) }
  \item{turbo}{\code{logical}. If \code{TRUE} (the default), then (parts of) paths are not searched twice (or more) unnecessarily, thus yielding a significant potential for speed-gain. The checking of whether the search has arrived at a point it has already been at comes with a slight computational overhead. So faster search is not guaranteed when \code{turbo=TRUE} }
  \item{print.searchinfo}{\code{logical}. If \code{TRUE} (default), then a print is returned whenever simiplification along a new path is started}
  \item{plot}{\code{NULL} or logical. If \code{TRUE}, then a plot is produced. If \code{NULL} (default), then the value set by \code{\link{options}} determines whether a plot is produced or not}
  \item{alarm}{\code{logical}. If \code{TRUE}, then a sound or beep is emitted (in order to alert the user) when the model selection ends }
  \item{\dots}{further arguments passed to or from other methods}
}

\details{No details form the moment}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{
\code{\link{dlogitx}}, \code{\link{dlogitxSim}}, \code{\link{coef.dlogitx}}, \code{\link{getsFun}}
}

\examples{
##no examples for the moment
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
