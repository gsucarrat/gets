\name{infocrit}
\alias{infocrit}
\alias{info.criterion}

\title{Computes the Average Value of an Information Criterion}

\description{
Given a log-likelihood, the number of observations and the number of estimated parameters, the average value of a chosen information criterion is computed. This facilitates comparison of models that are estimated with a different number of observations, e.g. due to different lags.
}

\usage{
infocrit(x, method=c("sc","aic","aicc","hq"))

info.criterion(logl, n=NULL, k=NULL, method=c("sc","aic","aicc","hq"))
}

\arguments{
  \item{x}{a \code{list} that contains, at least, three items: \code{logl} (a numeric, the log-likelihood), \code{k} (a numeric, usually the number of estimated parameters) and \code{n} (a numeric, the number of observations)}
  \item{method}{character, either "sc" (default), "aic", "aicc" or "hq"}
  \item{logl}{numeric, the value of the log-likelihood}
  \item{n}{integer, number of observations}
  \item{k}{integer, number of parameters}
}

\details{
Contrary to \code{\link{AIC}} and \code{\link{BIC}}, \code{info.criterion} computes the average criterion value (i.e. division by the number of observations). This facilitates comparison of models that are estimated with a different number of observations, e.g. due to different lags.
}

\value{
  \code{infocrit}: a numeric (i.e. the value of the chosen information criterion)

  \code{info.criterion}: a list with elements
    \item{method}{type of information criterion}
    \item{n}{number of observations}
    \item{k}{number of parameters}
    \item{value}{the value on the information criterion}
}

\references{
H. Akaike (1974): 'A new look at the statistical model identification'.
  IEEE Transactions on Automatic Control 19, pp. 716-723

E. Hannan and B. Quinn (1979): 'The determination of the order of an
  autoregression'. Journal of the Royal Statistical Society B 41, pp. 190-195

C.M. Hurvich and C.-L. Tsai (1989): 'Regression and Time Series Model
  Selection in Small Samples'. Biometrika 76, pp. 297-307

Pretis, Felix, Reade, James and Sucarrat, Genaro (2018): 'Automated General-to-Specific (GETS) Regression Modeling and Indicator Saturation for Outliers and Structural Breaks'. Journal of Statistical Software 86, Number 3, pp. 1-44

G. Schwarz (1978): 'Estimating the dimension of a model'. The Annals of
  Statistics 6, pp. 461-464
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
