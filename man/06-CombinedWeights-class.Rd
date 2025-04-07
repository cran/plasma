\name{CombinedWeights-class}
\docType{class}
\alias{CombinedWeights}
\alias{CombinedWeights-class}
\alias{summary,CombinedWeights-method}
\alias{image,CombinedWeights-method}
\alias{combineAllWeights}
\alias{stdize}
\alias{interpret}

\title{Class \code{"CombinedWeights"}}
\description{
  The \code{CombinedWeights} object class merges the weight matrices for
  all data sets in a plasma object.
}
\section{Objects from the Class}{
  Objects are defined using the \code{combineAllWeights} functions.
  Simply supply an object of class \code{plasma}. 
}
\usage{
combineAllWeights(pl)
\S4method{summary}{CombinedWeights}(object, \dots)
\S4method{image}{CombinedWeights}(x, \dots)
stdize(object, type = c("standard", "robust"))
interpret(object, component, alpha = 0.05)
}
\arguments{
  \item{pl}{An object of the \code{plasma} class.}
  \item{object}{An object of the \code{CombinedWeightss} class.}
  \item{x}{An object of the \code{CombinedWeightss} class.}
  \item{type}{A single character string indicating how to standardize
    the object. Legal value are "standard" or "robust".}
  \item{component}{A single chaaracter string; which componen should be
    interpreted.}
  \item{alpha}{A single numerical value between 0 and 1; what
    signfiicance value should be used to select important features.}
  \item{\dots}{Ignored; potentially, extra arguments to the summary or
    image methods.}
}
\section{Slots}{
  \describe{
      \item{\code{combined}:}{a matrix of the original variables in
	dataset \code{N} as rows and the PLS components \code{M} as
	columns.}
      \item{\code{featureSize}:}{a numeric (usually integer) vector that
	stores the number of features in each omics data set.}
      \item{\code{dataSource}:}{a factor indicating which omics data set
	each feature came from.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link[base]{summary}}:}{outputs summary statistics for the contributions of dataset \code{N} to components from all datasets in the case of \code{getAllWeights} or dataset \code{M} in the case of \code{getCompositeWeights}.}
  }
}
\value{
  The \code{combineAllWeights} function returns a newly constructed object of the
  \code{CombinedWeights} class. The summary method returna list
  containing four matrices. Each matrix has one row for each omics data
  set and one column for each model component. Each amtric contains
  different summary statistics, including the Mean, SD, Median, and MAD.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\examples{
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
# restrict data set size
MO <- with(plasmaEnv, prepareMultiOmics(
    assemble[c("ClinicalBin", "ClinicalCont", "RPPA")], Outcome))

splitVec <- with(plasmaEnv, rbinom(nrow(Outcome), 1, 0.6))
trainD <- MO[, splitVec == 1]
testD <- MO[, splitVec == 0]

firstPass <- fitCoxModels(trainD, "Days", "vital_status", "dead")
pl <- plasma(object = trainD, multi = firstPass)

getCompositeWeights(object = pl, N = "ClinicalBin", M = "RPPA")

cbin <- getAllWeights(object = pl, N = "ClinicalBin")
summary(cbin)
image(cbin)
heat(cbin, cexCol = 0.5)

cbin01 <- pickSignificant(object = cbin, alpha = 0.01)
image(cbin01)
heat(cbin01, cexCol = 0.5)

getTop(object = cbin01, N = 3)
}
\keyword{classes}
