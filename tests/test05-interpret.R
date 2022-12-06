library(plasma)
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
ls()
## prepare MultiOmics
MO <- with(plasmaEnv, prepareMultiOmics(assemble, Outcome) )
MO <- MO[c("ClinicalBin", "ClinicalCont", "RPPA"),]
summary(MO)
## split into train and test
set.seed(12345)
train <- rep(FALSE, 185)
train[sample(185, 113)] <- TRUE
MO2 <- MO[, train]
summary(MO2)
## test complete cox models
bigfit <- fitCoxModels(MO2, "Days", "vital_status", "dead")
## extend across dataset pairs
mfm <- plasma(MO2, bigfit)
plot(mfm)

cross <- getCompositeWeights(mfm, "ClinicalBin", "RPPA")
summary(cross)
image(cross)

heat(cross)

cb <- getAllWeights(mfm, "ClinicalBin")
summary(cb)
image(cb)
heat(cb)
getTop(cb)

rppaContrib <- getAllWeights(mfm, "RPPA")
getTop(rppaContrib, 5)
image(rppaContrib)
heat(rppaContrib)

Q <- max(abs(quantile(rppaContrib@contrib, c(0.05, 0.95))))
sig <- apply(abs(rppaContrib@contrib) > Q, 1, any)
smaller <- rppaContrib[sig,]
image(smaller)
heat(smaller)
