\name{Imputation}
\alias{meanModeImputer}
\alias{samplingImputer}

\title{Imputation}
\description{
  Functions to impute missing data in omics data sets.
}
\usage{
meanModeImputer(X)
samplingImputer(X)
}
\arguments{
  \item{X}{A numeric matrix, where the columns represent independent
    observations (patients or samples) and the columns represent measured
    features (genes, proteins, clinical variables, etc).}
}
\details{
  We recommend imputing small amounts of missing data in the input data
  sets when using the \code{plasma} package. The underlying issue is
  that the PLS models we use for individual omics data sets will not be
  able to make predictions on a sample if even one data point is
  missing. As a result, if a sample is missing at least one data point in
  every omics data set, then it will be impossible to use that sample at
  all.

  For a range of available imputation methods and R packages, consult
  the \href{https://CRAN.R-project.org/view=MissingData}{CRAN Task
  View on Missing Data}.  We also recommend the
  \href{https://rmisstastic.netlify.app/}{R-miss-tastic web site on
  missing data}.  Their simulations suggest that, for purposes of
  producing predictive models from omics data, the imputation method is
  not particularly important. Because of the latter finding, we have
  only implemented two simple imputation methods in the \code{plasma}
  package:
  \enumerate{
    \item  The \code{meanModeImputer} function will replace any missing data by the
      mean value of the observed data if there are more than five
      distinct values; otherwise, it will replace missing data by the
      mode.  This approach works relatively well for both continuous
      data and for binary or small categorical data.
    \item The \code{samplingImpute} function replaces missing values by sampling
      randomly from the observed data distribution.
  }
}
\value{
  Both functions return a numeric matrix of the same size and with the
  same row and column names as the input variable
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\examples{
loadESCAdata()
imputed <- with(plasmaEnv, lapply(assemble, samplingImputer) )
imputed <- with(plasmaEnv, lapply(assemble, meanModeImputer))
}
\keyword{array}
