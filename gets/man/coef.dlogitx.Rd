\name{coef.dlogitx}
\alias{coef.dlogitx}
\alias{fitted.dlogitx}
\alias{logLik.dlogitx}
\alias{plot.dlogitx}
\alias{print.dlogitx}
\alias{summary.dlogitx}
\alias{toLatex.dlogitx}
\alias{vcov.dlogitx}

\title{Extraction functions for 'dlogitx' objects}

\description{Extraction functions (of type S3 methods) for objects of class 'dlogitx'}

\usage{
  \method{coef}{dlogitx}(object, ...)
  \method{fitted}{dlogitx}(object, zero.prob=FALSE, ...)
  \method{logLik}{dlogitx}(object, ...)
  \method{plot}{dlogitx}(x, ...)
  \method{print}{dlogitx}(x, signif.stars=TRUE, ...)
  \method{summary}{dlogitx}(object, ...)
  \method{toLatex}{dlogitx}(object, digits = 4, gof = TRUE, nonumber = FALSE, nobs = "T", ...)
  \method{vcov}{dlogitx}(object, ...)
}

\arguments{
  \item{object}{an object of class 'dlogitx'}
  \item{x}{an object of class 'dlogitx'}
  \item{zero.prob}{\code{logical}. If \code{FALSE} (default), then the probabilities of a one are returned as fitted values. If \code{TRUE}, then the zero probabilities are returned as fitted values }
  \item{signif.stars}{\code{logical}. If \code{TRUE}, then p-values are additionally encoded visually, see \code{\link{printCoefmat}} }
  \item{digits}{integer, the number of digits in the LaTeX print}
  \item{gof}{logical that determines whether goodness-of-fit information should be included in the LaTeX print}
  \item{nonumber}{logical that determines whether a "nonumber" tag should be added to each equation in the LaTeX print}
  \item{nobs}{character that determines the label for the number of observations in the LaTeX print}
  \item{\dots}{additional arguments}
}

\value{Various, depending on the method}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{
\code{\link{dlogitx}}
}

\examples{
##no examples for the moment
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
