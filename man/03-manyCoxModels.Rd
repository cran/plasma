\name{MultiplePLSCoxModels-class}
\docType{class}
\alias{MultiplePLSCoxModels}
\alias{MultiplePLSCoxModels-class}
\alias{summary,MultiplePLSCoxModels-method}
\alias{[,MultiplePLSCoxModels,ANY,ANY,ANY-method}%]
\alias{plot,MultiplePLSCoxModels,missing-method}
\alias{predict,MultiplePLSCoxModels-method}
\alias{validMultipleCoxModels}
\alias{fitCoxModels}
\alias{getSizes}

\title{Class \code{"MultiplePLSCoxModels"}} 
\description{
  The \code{MultiplePLSCoxModels} object class ...
  The \code{validMultipleCoxModels} function checks if each data set contains the same set of samples.
  The \code{fitCoxModels} function fits many plsRcoxmodels and returns an S4 object of class \code{MultiplePLSCoxModels}.
  The \code{getSizes} function returns a matrix with the list of dataframes of the \code{MultiOmics} object as rownames and columns with NT, cNT, and p-values.
}

\usage{
fitCoxModels(multi, timevar, eventvar, eventvalue, verbose)
\S4method{summary}{MultiplePLSCoxModels}(object, \dots)
\S4method{plot}{MultiplePLSCoxModels,missing}(x, y, col = c("blue", "red"),
      lwd = 2, xlab = "", ylab = "Fraction Surviving",
       mark.time = TRUE, legloc = "topright", \dots)
\S4method{predict}{MultiplePLSCoxModels}(object, newdata, type = c("components", "risk",
                 "split", "survfit"), \dots)
}
\arguments{
  \item{multi}{an object of class \code{MultiOmics} for fitting the model.}
  \item{timevar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the time-to-event.}
  \item{eventvar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the event.}
  \item{eventvalue}{a character string specifying the value of the event in \code{eventvar}.}
  \item{verbose}{logical; should the function report progress.}
  \item{object}{an object of class \code{MultiplePLSCoxModels} for outputting the summary.}
  \item{x}{an object of class \code{MultiplePLSCoxModels} for plotting the Kaplan-Meier curves.}
  \item{y}{An ignored argrument for the plot method.}
  \item{col}{A vector of color specifications. Default is c(\dQuote{blue}, \dQuote{red}).}
  \item{lwd}{A vector specifying the line width. Default is \dQuote{2}.}
  \item{xlab}{A character string to label the x-axis. Default is \dQuote{}.}
  \item{ylab}{A character string to label the y-axis. Default is \dQuote{Fraction Surviving}.}
  \item{mark.time}{A logical value; should tickmarks indicate censored
    data? Default is TRUE.}
  \item{legloc}{A character string indicating where to put the legend. Default is \dQuote{topright}.}
  \item{\dots}{Other graphical parameters.}
  \item{newdata}{A \code{MultiOmics} object with the same structure as
    the training data.}
  \item{type}{An enumerated character value.}

}
\section{Slots}{
  \describe{
      \item{\code{models}:}{A list of \code{SingleModel} objects, one for each assay.}
      \item{\code{timevar}:}{A character matching the name of the column containing the time-to-event.}
      \item{\code{eventvar}:}{A character matching the name of the column containing the event.}
      \item{\code{eventvalue}:}{A character specifying the event in \code{eventvar.} }
  }
}

\section{Methods}{
  \describe{
      \item{\code{\link[base]{plot}}:}{Plots Kaplan-Meier curves for each
	omics dataset split into Low Risk and High Risk groups.}
      \item{\code{\link[base]{summary}}:}{Returns a description of the
	\code{MultiplePLSCoxModels} object and the names of the omics
	datasets used to build the model.}
      \item{\code{\link[stats]{predict}}:}{usually returns a list of numeric
  vectors of predicted risk per data type. When \code{type = "survfit"},
  retuns a list of \code{survfit} objects.}
  }
}
\value{
  The \code{fitCoxModels} function retuns a newly constructed object of
  the \code{MultiplePLSCoxModels} class. The \code{plot} method
  invisibly returns the object on which it was invoked. The
  \code{summary} method returns no value. The predict method returns a
  list of prediction results, each of which comes from the
  \code{predict} method for the \code{\link{SingleModel-class}.}
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\seealso{fitSingleModel}
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

summary(firstPass)
plot(firstPass)
getSizes(firstPass)
pre1 <- predict(firstPass, testD)
}
\keyword{classes}
