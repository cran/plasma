library(plasma)
## check data sets
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
ls()
## make sure we can assemble MultiOmics objects
MO <- with(plasmaEnv, prepareMultiOmics(assemble, Outcome))
## Test summary and plot methods
summary(MO)
opar <- par(mai = c(1.02, 1.82, 0.82, 0.42))
plot(MO)
par(opar)
## Make sure we can use the subset operator
train <- rep(c(TRUE, FALSE), times = c(112, 185-112))
MO2 <- MO[, train]
summary(MO2)
## For future reference, show there are no survival differences
if (require("survival")) {
  plot(survfit(Surv(Days, vital_status == "dead") ~ train, data = MO@outcome),
               col=c("cyan4", "magenta"), lwd=3)
}
## And check that the other argument also works
MO3 <- MO[c("ClinicalBin", "RPPA"),]
summary(MO3)

## Things expected to fail
## disorient one of the component datasets.
badinput <- with(plasmaEnv, assemble)
badinput$mRNASeq <- t(badinput$mRNASeq)
pmo <- try( with(plasmaEnv, prepareMultiOmics(badinput, Outcome) ))
