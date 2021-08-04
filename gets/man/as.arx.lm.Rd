\name{as.arx}
\alias{as.arx}
\alias{as.arx.lm}

\title{Convert an object to class 'arx'}

\description{
The function \code{as.arx} is a generic function and its methods returns an object of class \code{\link{arx}}.
}

\usage{
as.arx(object, ...)

##S3 method for objects of class 'lm':
\method{as.arx}{lm}(object, ...) 
}

\arguments{
  \item{object}{object of class \code{\link{lm}} }
  \item{...}{arguments passed on to and from other methods }
}

\value{Object of class \code{\link{arx}} }

\author{
Genaro Sucarrat \url{http://www.sucarrat.net/}\cr\cr
}

\seealso{
\code{\link{lm}}, \code{\link{arx}}
}

\examples{

##generate some data:
set.seed(123) #for reproducibility
y <- rnorm(30) #generate Y
x <- matrix(rnorm(30*10), 30, 10) #create matrix of Xs

##typical situation:
mymodel <- lm(y ~ x)
as.arx(mymodel)
                                 
##use hetero-robust vcov:
as.arx(mymodel, vcov.type="white")

##add ar-dynamics:
as.arx(mymodel, ar=1:2)

##add log-variance specification:
as.arx(mymodel, arch=1:2)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
