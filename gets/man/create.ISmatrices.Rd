\name{create.ISmatrices}
\alias{create.ISmatrices}

\title{
Internal function to construct Indicator Saturation Matrices
}
\description{
Internal function that constructs a list of indicator matrices used in indicator 
saturation (IIS, SIS, TIS, UIS) to be passed back towards \code{\link{isat}}. 
}
\usage{
create.ISmatrices(iis, sis, tis, uis, y.n, y.index.as.char)
}

\arguments{
  \item{iis}{
    Logical values indicating whether to generate Impulse indicators.
}
  \item{sis}{
  Logical values indicating whether to generate Step indicators.
}
  \item{tis}{
  Logical values indicating whether to generate Trend indicators.
}
  \item{uis}{
  Taking the \code{uis} argument from \code{\link{isat}} for to process further.
}
  \item{y.n}{
  Sample size.
}
  \item{y.index.as.char}{
  The index of y as a character vector to retain any implied time series
  structure.
}
}

\value{A list object containing the required indicator matrices. 
}

\seealso{
 See Also \code{\link{isat}}
}


