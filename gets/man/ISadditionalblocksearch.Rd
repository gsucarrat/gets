\name{ISadditionalblocksearch}
\alias{ISadditionalblocksearch}

\title{Internal function to carry out additional Block Search in 
isat: Leave-One-Out Block Search for Excessive Indicators}

\description{
  Internal function that handles excessive indicator retention by re-running 
  GETS on smaller blocks of indicators using a leave-one-out method. This 
  can be necessary when following a block search, e.g. for IIS, more indicators
  than observations have been identified. Then another block search is carried 
  out. This can be avoided by choosing tighter significance levels in 
  \code{t.pval} or improving the specification.
}
\usage{
ISadditionalblocksearch(mXis, isNames, y, y.n, mX, mXnames, mxkeep, 
  ISmatrixname, indicator.set, estimations.counter, getsFun.counter, 
  print.searchinfo, tol, LAPACK, userEstArg, t.pval, gof.method, 
  gofFunArg, max.regs, max.paths, arLjungB, archLjungB, normality.JarqueB, 
  user.diagnostics, include.gum, include.1cut, include.empty, turbo, do.pet, 
  ratio.threshold, max.block.size, wald.pval, additional.block.search)
}

\arguments{
  \item{mXis}{
    Initial regressor matrix including previously selected indicators.
  }
  \item{isNames}{
    Character vector with names of the candidate indicator variables.
  }
  \item{y}{
    Response variable, typically a numeric vector or time-series.
  }
  \item{y.n}{
    Integer indicating the number of observations in the sample.
  }
  \item{mX}{
    Matrix of original regressors, excluding indicators.
  }
  \item{mXnames}{
    Character vector of column names corresponding to \code{mX}.
  }
  \item{mxkeep}{
    Indices of regressors to always retain during selection.
  }
  \item{ISmatrixname}{
    Character string naming the current IS matrix (e.g., "SIS", "TIS").
  }
  \item{indicator.set}{
    Matrix of candidate indicators from which subsets are selected.
  }
  \item{estimations.counter}{
    Running counter of the number of model estimations performed.
  }
  \item{getsFun.counter}{
    Running counter of how many times \code{getsFun} was called.
  }
  \item{print.searchinfo}{
    Logical; if \code{TRUE}, outputs progress messages during execution.
  }
  \item{tol}{
    Numeric tolerance used in collinearity checks (see \code{qr}).
  }
  \item{LAPACK}{
    Logical; if \code{TRUE}, LAPACK is used for QR decomposition.
  }
  \item{userEstArg}{
    List specifying the user-defined estimation function and its arguments.
  }
  \item{t.pval}{
    Significance level for coefficient t-tests in the GETS algorithm.
  }
  \item{gof.method}{
    Character indicating whether to minimize or maximize the 
    goodness-of-fit criterion.
  }
  \item{gofFunArg}{
    List containing the name and arguments of the GOF function used.
  }
  \item{max.regs}{
    Maximum number of regressions allowed per deletion path.
  }
  \item{max.paths}{
    Maximum number of GETS paths to be explored.
  }
  \item{arLjungB}{
    List with order and significance level for the Ljung-Box test on 
    residual autocorrelation.
  }
  \item{archLjungB}{
    List with order and significance level for the Ljung-Box test on 
    squared residuals (ARCH).
  }
  \item{normality.JarqueB}{
    Significance level for the Jarque-Bera test for normality; 
    if \code{NULL}, test is skipped.
  }
  \item{user.diagnostics}{
    List containing user-specified diagnostic function and p-value thresholds.
  }
  \item{include.gum}{
    Logical; if \code{TRUE}, the General Unrestricted Model is included 
    in terminal models.
  }
  \item{include.1cut}{
    Logical; if \code{TRUE}, the 1-cut model is included in terminal models.
  }
  \item{include.empty}{
    Logical; if \code{TRUE}, the empty model is included in terminal models.
  }
  \item{turbo}{
    Logical; avoids redundant GETS path evaluations if \code{TRUE}.
  }
  \item{do.pet}{
    Logical; if \code{TRUE}, performs Parsimonious Encompassing Tests (PETs).
  }
  \item{ratio.threshold}{
    Numeric ratio used to determine block sizes in relation to sample size.
  }
  \item{max.block.size}{
    Maximum number of indicators allowed per block during the 
    leave-one-out search.
  }
  \item{wald.pval}{
    Significance threshold for Wald tests in the GETS algorithm.
  }
  \item{additional.block.search}{Default is \code{FALSE}. If \code{isat()} identifies (retains) too many indicators in the block search, a sensible model might not be estimable (e.g. due to a lack of available degrees of freedom). In this case, setting \code{additional.block.search = TRUE} carries out an additional block search using a Leave-one-Out method (e.g. for three blocks, indicators are assigned to block \code{1,2,3,1,2,3,1,2,3}). This is an experimental feature.}
}

\value{
A list containing:
  \item{ISfinalmodel}{The final selected model.}
  \item{estimations.counter}{Updated estimation counter.}
  \item{getsFun.counter}{Updated function call counter.}
  \item{ISblocks}{Final block structure.}
}

\author{
Moritz Schwarz, \url{https://www.moritzschwarz.org/}\cr\cr
}


\seealso{\code{\link{ISMatricesLoop}}, \code{\link{isat}}}

