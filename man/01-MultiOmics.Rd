\name{MultiOmics-class}
\docType{class}
\alias{MultiOmics}
\alias{MultiOmics-class}
\alias{summary,MultiOmics-method}
\alias{[,MultiOmics,ANY,ANY,ANY-method}%]
\alias{plot,MultiOmics,missing-method}
\alias{validMultiOmics}
\alias{prepareMultiOmics}

\title{Class \code{"MultiOmics"}}
\description{
  The \code{prepareMultiOmics} function returns a new object of \code{MultiOmics} class for use in \code{fitCoxModel}.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{prepareMultiOmics} constructor. In
  the simplest case, you enter two objects: a list of dataframes and a dataframe of clinical outcomes.
}
\usage{
prepareMultiOmics(datalist, outcome)
\S4method{summary}{MultiOmics}(object, \dots)
\S4method{plot}{MultiOmics,missing}(x, y, \dots)
}
\arguments{
  \item{datalist}{a list of dataframes formatted to have variables as rows (dimension D) and samples as columns (dimension N).}
  \item{outcome}{a dataframe of clinical outcomes formatted to have
    sample names as row indexes and variable names as column indexes}
  \item{object}{An object of the \code{MultiOmics} class.}
  \item{x}{An object of the \code{MultiOmics} class.}
  \item{y}{Nothing; ignored.}
  \item{\dots}{Extra graphical or other parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{data}:}{A list of dataframes with variables as rows or varying length and samples as columns of uniform length N, where N is the maximum value of non-missing samples in any given dataset. Note that \code{NA}s have been added to \dQuote{pad} to make the column length uniform across data types.}
      \item{\code{outcome}:}{A dataframe of clinical outcomes with variables as columns and samples as rows.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link[base]{plot}}:}{Produces a visual representation of the  dimensionalities of each dataframe in datalist. D corresponds to the number of variables in each omics dataframe, and N corresponds to samples (or members) whose variable is not entirely missing. Gray areas correspond to missing samples.}
      \item{\code{\link[base]{summary}}:}{Produces summary tables corresponding to datasets and outcomes.}
  }
}
\value{
The \code{prepareMultiOmics} function returns a new object of the \code{MultiOmics} class.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\examples{
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
MO <- with(plasmaEnv,
           prepareMultiOmics(datalist = assemble, outcome = Outcome))
plot(MO)
summary(MO)
}
\keyword{classes}
