library(plasma)
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
ls()
## prepare MultiOmics
MO <- with(plasmaEnv, prepareMultiOmics(assemble, Outcome))
MO <- MO[c("ClinicalBin", "ClinicalCont", "RPPA"),]
summary(MO)
## test complete cox models
set.seed(12345)
train <- rep(FALSE, 185)
train[sample(185, 113)] <- TRUE
MO2 <- MO[, train]
summary(MO2)
bigfit <- fitCoxModels(MO2, timevar = "Days",
                       eventvar = "vital_status",
                       eventvalue = "dead")
## test core methods
class(bigfit)
summary(bigfit)
getSizes(bigfit)
plot(bigfit)
preds <- predict(bigfit)
## make sure we can predict on new data
testobj <- MO[, !train]

predsB <- predict(bigfit, newdata = testobj)
predsR <- predict(bigfit, newdata = testobj, type = "risk")
predsS <- predict(bigfit, newdata = testobj, type = "split")
