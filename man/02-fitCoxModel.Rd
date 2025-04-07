\name{SingleModel-class}
\docType{class}
\alias{fitSingleModel}
\alias{SingleModel-class}
\alias{predict,SingleModel-method}
\alias{summary,SingleModel-method}
\alias{plot,SingleModel,missing-method}
\alias{fitSingleModel}

\title{Class \code{"SingleModel"}}
\description{
  The \code{fitSingleModel} function takes in an object of
  \code{MultiOmics} class and returns a new object of
  \code{SingleModel} class. 
}

\usage{
fitSingleModel(multi, N, timevar, eventvar, eventvalue)
\S4method{summary}{SingleModel}(object, \dots)
\S4method{plot}{SingleModel,missing}(x, y, col = c("blue", "red"),
      lwd = 2, xlab = "", ylab = "Fraction Surviving",
       mark.time = TRUE, legloc = "topright", \dots)
\S4method{predict}{SingleModel}(object, newdata, type = c("components", "risk",
                 "split", "survfit"), \dots)
}
\arguments{
  \item{multi}{an object of class \code{MultiOmics} for fitting the model.}
  \item{N}{A character string identifying the data set being modeled.}
  \item{timevar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the time-to-event.}
  \item{eventvar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the event.}
  \item{eventvalue}{a character string specifying the value of the event.}
  
  \item{x}{an object of class \code{plsRcoxmodel} for plotting the Kaplan-Meier curves.}
  \item{y}{An ignored argrument for the plot method.}
  \item{col}{A vector of color specifications.}
  \item{lwd}{A vactor specifying the line width.}
  \item{xlab}{A character string to label the x-axis.}
  \item{ylab}{A character string to label the y-axis.}
  \item{mark.time}{A logical value; should tickmarks indicate censored
    data?}
  \item{legloc}{A character string indicating where to put the legend.}

  \item{object}{an object of class \code{SingleModel}.}
  \item{newdata}{A \code{MultiOmics} object with the same structure as
    the training data.}
  \item{type}{An enumerated character value.}
  \item{\dots}{other parameters used in graphing or prediction.}
  
}
\section{Slots}{
  \describe{
      \item{\code{plsmod}:}{Object of class \code{plsRcoxmodel}
	containing the fitted model.}
      \item{\code{Xout}:}{Object of type \code{data.frame} containing
	the original \code{outcome} dataframe and additional columns for
	"Risk", and "Split", corresponding to the risk of the event
	calculated by the model, and patient assignment to low versus
	high-risk groups, respectively.}
      \item{\code{dsname}:}{A character vector of length one; the name
	of the data set being modeled from a \code{MultiOmics} object.}
      \item{\code{SF}:}{Object of type \code{survfit} which is used by the \code{plot} method to plot Kaplan-Meier curves grouped by predicted Split. See documentation for \code{link{survfit}}.}
      \item{\code{riskModel}:}{Object of type \code{coxph} that uses predicted Risk (continuous) as the predictor variable and survival as the response variable. See documentation for \code{link{coxph}}.}
      \item{\code{splitModel}:}{Object of type \code{coxph} that uses predicted Split (predicted Risk categorized into \dQuote{high} and \dQuote{low} risk by the median predicted Risk) as the predictor variable and survival as the response variable. See documentation for \code{link{coxph}}.}

  }
}


\section{Methods}{
  \describe{
      \item{\code{\link[base]{plot}}:}{Plots Kaplan-Meier curves for each
	omics dataset split into Low Risk and High Risk groups.}
      \item{\code{\link[base]{summary}}:}{Returns a description of the
	\code{MultiplePLSCoxModels} object and the names of the omics
	datasets used to build the model.}
      \item{\code{\link[stats]{predict}}:}{Usually, a numeric vector containing
        the predicted risk values. However, when using \code{type =
	"survfit"}, tghe return value is a \code{survfit} object from
	the\code{survival} package.}
  }
}
\value{
  The \code{fitSingleModel} function returns a newly constructed object
  of the \code{SingleModel} class. The \code{plot} method invisibly
  returns the value on which it was invoked. The \code{summary} method
  returns an object summarizing the final model produced by PLS R cox
  regression. The \code{predict} method returns either a vector or
  matrix depending on the type of predictions requested.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\seealso{
  \code{\link{getSizes}}
}
\examples{
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
MO <- with(plasmaEnv, prepareMultiOmics(assemble, Outcome) )
MO <- MO[c("ClinicalBin", "ClinicalCont", "RPPA"),]
set.seed(98765)
splitVec <- with(plasmaEnv, rbinom(nrow(Outcome), 1, 0.6))
trainD <- MO[, splitVec == 1]
testD <- MO[, splitVec == 0]

zerothPass <- fitSingleModel(trainD, N = "RPPA",
                           timevar = "Days", eventvar = "vital_status",
                           eventvalue = "dead")
summary(zerothPass)
plot(zerothPass)
pre0 <- predict(zerothPass, testD)
}
\keyword{classes}
