\name{plasmaPredictions-class}
\docType{class}
\alias{plasmaPredictions}
\alias{plasmaPredictions-class}
\alias{plot,plasmaPredictions,missing-method}

\title{Class \code{"plasmaPredictions"}}
\description{
  The \code{plasmaPredictions} object class is returned when running the
  \code{predict} method on an object of class \code{plasma}.
}
\section{Objects from the Class}{
  Users shold not create objects of this class directly. They will be
  automatically created when you apply the \code{predict} method to a
  fully worked out \code{plasma} model.
}
\usage{
\S4method{plot}{plasmaPredictions,missing}(x, y, col = c("blue", "red"),
                 lwd = 2, xlab = "", ylab = "Fraction Surviving",
                 mark.time = TRUE, legloc = "topright", \dots)
}
\arguments{
  \item{x}{An object of the \code{plasmaPredictions} class for plotting the Kaplan-Meier curves.}
  \item{y}{An ignored argument for the plot method.}
  \item{col}{A vector of color specifications. Default is c(\dQuote{blue}, \dQuote{red}).}
  \item{lwd}{A vactor specifying the line width. Default is \dQuote{2}.}
  \item{xlab}{A character string to label the x-axis. Default is \dQuote{}.}
  \item{ylab}{A character string to label the y-axis. Default is \dQuote{Fraction Surviving}.}
  \item{mark.time}{A logical value; should tickmarks indicate censored
    data? Default is TRUE.}
  \item{legloc}{A character string indicating where to put the legend. Default is \dQuote{topright}.}
  \item{\dots}{Other graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{meanPredictions}:}{A matrix with samples as rows and
	factors as columns that is a result of taking the mean of the
	PLS component predictions from each dataset.}
      \item{\code{riskDF}:}{Object of type \code{data.frame} containing
	the original \code{outcome} dataframe and additional columns for
	"Risk", and "Split", corresponding to the risk of the event
	calculated by the model, and patient assignment to low versus
	high-risk groups, respectively.}
      \item{\code{riskModel}:}{Object of type \code{coxph} that uses
	predicted Risk (continuous) as the predictor variable and
	survival as the response variable. See documentation for
	\code{link{coxph}}.} 
      \item{\code{splitModel}:}{Object of type \code{coxph} that uses
	predicted Split (predicted Risk categorized into \dQuote{high}
	and \dQuote{low} risk by the median predicted Risk) as the
	predictor variable and survival as the response variable. See
	documentation for \code{link{coxph}}.} 
      \item{\code{SF}:}{Object of type \code{survfit} which is used by
	the \code{plot} method to plot Kaplan-Meier curves grouped by
	predicted Split. See documentation for \code{link{survfit}}.} 
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link[base]{plot}}:}{Produces Kaplan-Meier curves for the
	low risk and high risk groups.}
  }
}
\value{
  The \code{predict} method on an object of the \code{\link{plasma}}
  class returns an object of the \code{\link{plasmaPredictions}}
  class. The \code{plot} method invisibly returns the value on which it
  was invoked.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\note{
  An object of \code{plasmaPredictions} class contains many models that
  are similar to an object of \code{\link{MultiplePLSCoxModels}} class.
}
\seealso{plasma}
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

testpred <- predict(pl, testD)
plot(testpred, main = "Testing", xlab = "Time (Days)")
}
\keyword{classes}
