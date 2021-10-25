\name{coef.logitx}
\alias{coef.logitx}
\alias{fitted.logitx}
\alias{logLik.logitx}
\alias{plot.logitx}
\alias{print.logitx}
\alias{summary.logitx}
\alias{toLatex.logitx}
\alias{vcov.logitx}

\title{Extraction functions for 'logitx' objects}

\description{Extraction functions (of type S3 methods) for objects of class 'logitx'}

\usage{
  \method{coef}{logitx}(object, ...)
  \method{fitted}{logitx}(object, zero.prob=FALSE, ...)
  \method{logLik}{logitx}(object, ...)
  \method{plot}{logitx}(x, ...)
  \method{print}{logitx}(x, signif.stars=TRUE, ...)
  \method{summary}{logitx}(object, ...)
  \method{toLatex}{logitx}(object, digits = 4, gof = TRUE, nonumber = FALSE, nobs = "T", ...)
  \method{vcov}{logitx}(object, ...)
}

\arguments{
  \item{object}{an object of class 'logitx'}
  \item{x}{an object of class 'logitx'}
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
\code{\link{logitx}}, \code{\link{logitxSim}}, \code{\link{gets.logitx}}
}

\examples{

##simulate from ar(1):
set.seed(123) #for reproducibility
y <- logitxSim(100, ar=0.3)

##estimate and store result:
mymod <- logitx(y, ar=1)

##extract stuff:
coef(mymod)
fitted(mymod)
logLik(mymod)
plot(mymod)
print(mymod)
summary(mymod)
toLatex(mymod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}