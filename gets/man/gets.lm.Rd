\name{gets.lm}
\alias{gets.lm}

\title{General-to-Specific (GETS) Modelling 'lm' objects}

\description{
General-to-Specific (GETS) Modelling of objects of class \code{\link{lm}}.
}

\usage{
\method{gets}{lm}(x, keep = NULL, include.1cut = TRUE, print.searchinfo = TRUE, ...) 
}

\arguments{
  \item{x}{an object of class 'lm', see \code{\link{lm}} }
  \item{keep}{\code{NULL} or a vector of integers that determines which regressors to be excluded from removal in the specification search }
  \item{include.1cut}{\code{logical}. If \code{TRUE} (default), then the 1-cut model is added to the list of terminal models. If \code{FALSE}, then the 1-cut is not added, unless it is a terminal model in one of the paths }
  \item{print.searchinfo}{\code{logical}. If \code{TRUE} (default), then selected info is printed during search }
  \item{\dots}{further arguments passed on to \code{\link{getsFun}} }
}


\details{
Internally, \code{gets.lm} invokes \code{\link{getsFun}} for the GETS-modelling, which is also invoked by \code{\link{getsm}}. See their help pages for more information.
}

\value{A list of class \code{\link{lm}}. Note that the 'top' of the list contains information (paths and terminal models) from the GETS modelling, see \code{\link{paths}} and \code{\link{terminals}} }

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{
\code{\link{lm}}, \code{\link{getsFun}}, \code{\link{getsm}}, \code{\link{paths}} and \code{\link{terminals}}
}

\examples{

##generate some data:
set.seed(123) #for reproducibility
y <- rnorm(30) #generate Y
x <- matrix(rnorm(30*10), 30, 10) #matrix of Xs
colnames(x) <- paste0("var", 1:NCOL(x))

##estimate model:
mymod <- lm(y ~ x)

##do gets modelling:
gets(mymod)

##ensure intercept is not removed:
gets(mymod, keep=1)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
