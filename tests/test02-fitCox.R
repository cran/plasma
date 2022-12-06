library(plasma)
# Repeat basic stuff from first test
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
ls()
MO <- with(plasmaEnv, prepareMultiOmics(assemble, Outcome))
train <- rep(c(TRUE, FALSE), times = c(112, 185-112))
MO2 <- MO[, train]
## Fit a survival model on a single data set on the MultiOmics object
fitted <- fitSingleModel(MO2, "ClinicalBin", "Days", "vital_status", "dead")
summary(fitted)
p <- predict(fitted)
summary(p)
plot(fitted, xlab = "Time (Days)", legloc = "topright", main = "Training Data")
## Make sure we can predict on new data
testobj <- MO[, !train]
summary(testobj)
## How do we get rid of the "Missing value" information?
pre <- predict(fitted, newdata = testobj)
prer <- predict(fitted, newdata = testobj, type = "risk")
pres <- predict(fitted, newdata = testobj, type = "split")
pairs(cbind(pre, prer, pres))

## Things that should fail:
p <- try( predict(fitted, "risk") )
p <- try( predict(fitted, type = "riak") )
