\name{Contribution-class}
\docType{class}
\alias{Contribution}
\alias{Contribution-class}
\alias{[,Contribution,ANY,ANY,ANY-method}%]
\alias{summary,Contribution-method}
\alias{image,Contribution-method}
\alias{heat}
\alias{heat,Contribution-method}
\alias{getCompositeWeights}
\alias{getAllWeights}
\alias{getFinalWeights}
\alias{getTop}
\alias{pickSignificant}

\title{Class \code{"Contribution"}}
\description{
  The \code{Contribution} object class contains the weight matrix between variables and the PLS components. The values in the weight matrix are a numeric representation of how much a variable from the omics datasets contributed to defining the final PLS components.
}
\section{Objects from the Class}{
  Objects are defined using the \code{getAllWeights}, \code{getCompositeWeights}, \code{getTop}, or \code{pickSignificant} functions. In the simplest scenario, one would enter an object of class \code{plasma} and any specific parameters associated with the function (see arguments section for more info).
}
\usage{
getCompositeWeights(object, N, M)
getAllWeights(object, N)
getFinalWeights(object)
getTop(object, N = 1)
pickSignificant(object, alpha)
\S4method{summary}{Contribution}(object, \dots)
\S4method{image}{Contribution}(x, col = viridis(64), mai = c(1.82, 1.52, 0.32, 0.32), \dots)
\S4method{heat}{Contribution}(object, main = "Contributions", col = viridis(64),
                                           mai = c(1.52, 0.32, 0.82, 1.82), \dots)
}
\arguments{
  \item{object}{In the first four functions, an object of the
    \code{plasma} class. In the methods described here, an object
    of the \code{Contributions} class.}
  \item{N}{in the function \code{getCompositeWeights}, the name of the dataset being modeled. \cr in the function \code{getTop}, the number of significant components you want to print.}
  \item{M}{name of the dataset being modeled pairwise with dataset \code{N} in the \code{getCompositeWeights} function.}
  \item{alpha}{level of significance used in the \code{pickSignificant} function.}
  \item{\dots}{other graphical parameters.}
  \item{x}{an object of the \code{Contributions} class.}
  \item{main}{A character vector of length one; the main plot title.}
  \item{col}{A vector of color descriptors.}
  \item{mai}{A vector of four nonnegative numbers.}
}
\section{Slots}{
  \describe{
      \item{\code{contrib}:}{a matrix of the original variables in dataset \code{N} as rows and the PLS components \code{M} as columns.}
      \item{\code{datasets}:}{a character vector that stores the names of the datasets that were specified for the function.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{summary}}:}{outputs summary statistics for the contributions of dataset \code{N} to components from all datasets in the case of \code{getAllWeights} or dataset \code{M} in the case of \code{getCompositeWeights}.}
      \item{\code{\link{image}}:}{outputs a heatmap of the transposed \code{contrib} matrix.}
      \item{\code{heat}:}{outputs a clustered heatmap of the \code{contrib} matrix.}
  }
}
\value{
  The \code{plasma} function returns a newly constructed object of the
  \code{plasma} class. 
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
