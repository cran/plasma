\name{plasma-class}
\docType{class}
\alias{plasma-class}
\alias{plot,plasma,missing-method}
\alias{barplot,plasma-method}
\alias{summary,plasma-method}
\alias{predict,plasma-method}
\alias{plasma}

\title{Class \code{"plasma"}}
\description{
  The \code{plasma} object class is returned after running the \code{plasma} function.
  The \code{plasma} function uses the \code{PLSRCox} components from one
  dataset as the predictor variables and the \code{PLSRCox} components
  of another dataset as the response variables to fit a partial least
  squares regression (\code{plsr}) model. Then, we take the mean of the
  predictions to create a final matrix of samples versus components.

  The matrix of components described earlier is then used to fit a Cox
  Proportional Hazards (\code{coxph}) model with AIC stepwise variable
  selection to return a final object of class \code{plasma} which
  includes a \code{coxph} model with a reduced number of predictors.
}
\section{Objects from the Class}{
Objects should be defined using the \code{plasma} function.
}
\usage{
plasma(object, multi)
\S4method{plot}{plasma,missing}(x, y, \dots)
\S4method{barplot}{plasma}(height, source, n, direction = c("both", "up","down"),
                lhcol = c("cyan", "red"), \dots)
\S4method{predict}{plasma}(object, newdata = NULL, type = c("components", "risk",
                 "split"), \dots)
}
\arguments{
  \item{multi}{an object of the \code{MultiplePLSCoxModels} class.}
  \item{object}{an object of the \code{plasma} class.}
  \item{height}{an object of the \code{plasma} class for the \code{barplot} method.}
  \item{x}{an object of class \code{plasma} for plotting the Kaplan-Meier curves.}
  \item{y}{An ignored argrument for the plot method.}
  \item{source}{A length-one character vector; the name of a data set in
    a \code{plasma} object.}
  \item{n}{A length-one integer vector; the number of high-weight
    features to display.}
  \item{direction}{A length-one character vector; show features with
    positive weights (up), negative (down), or both.}
  \item{lhcol}{A chaacter vector of length 2, indicating the preferred
  colors for low (negative) or high (positive) weights.}
  \item{newdata}{A \code{MultiOmics} object with the same structure as
    the training data.}
  \item{type}{An enumerated character value.}
  \item{\dots}{Additional graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{traindata}:}{An object of class \code{MultiOmics} used for training the model.}
      \item{\code{compModels}:}{A list containing objects in the form of \code{plsr}.}
      \item{\code{fullModel}:}{A coxph object with variables (components) selected via AIC stepwise selection.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{plot}}:}{Plots a Kaplan-Meier curve of the final \code{coxph} model that has been categorized into \dQuote{low risk} and \dQuote{high risk} based whether it is higher or lower, respectively, than the median value of risk.}
      \item{\code{\link{predict}}:}{creates an object of class \code{plasmaPredictions}.}
      \item{\code{\link{barplot}}:}{Produces a barplot of the \code{n} largest
	weights assigned to features from the appropriate data \code{source}.}
  }
}
\value{
  The \code{plasma} function returns a newly constructed object of the \code{plasma} class. The \code{plot} method invisibly returns the object on which it was invoked. The \code{predict} method returns an object of the \code{\link{plasmaPredictions}} class.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\seealso{plasmaPredictions, plsr}
\examples{
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
# restrict data set size
MO <- with(plasmaEnv, prepareMultiOmics(
    assemble[c("ClinicalBin", "ClinicalCont", "RPPA")], Outcome) )

splitVec <- with(plasmaEnv, rbinom(nrow(Outcome), 1, 0.6))
trainD <- MO[, splitVec == 1]
testD <- MO[, splitVec == 0]

firstPass <- fitCoxModels(trainD, "Days", "vital_status", "dead")
pl <- plasma(object = trainD, multi = firstPass)

plot(pl, legloc = "topright", main = "Training Data")
barplot(pl, "RPPA", 6)
barplot(pl, "RPPA", 10, "up")
}
\keyword{classes}
