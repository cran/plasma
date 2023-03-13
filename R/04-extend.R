extendCoxModels <- function(object, firstPass, verbose = TRUE) {
  ## get the components from each data set
  Components <- lapply(firstPass@models, function(DS) {
    DS@plsmod$tt
  })
  tempComps <- t(sapply(Components, dim))
  dimnames(tempComps) <- list(names(Components), c("nPats", "nComps"))
  sizing <- apply(tempComps, 2, sum)
  ## build a storage area to hold the combined components
  stores <- matrix(NA, nrow(object@outcome), sizing[2])
  rownames(stores) <- rownames(object@outcome)
  compLevel <- unlist(lapply(tempComps[,2], function(K) 1:K))
  NN <- names(compLevel)
  ## deal with special case of only one component
  nodigit <- !grepl("\\d", NN) # doesn't end with a digit
  if (any(nodigit)) {
    NN[nodigit] <- paste(NN[nodigit], "1", sep = "")
  }
  colnames(stores) <- NN
  ## fill that storage area
  for (N in names(Components)) {
    if(verbose) cat(N, "\n", file = stderr())
    X <- Components[[N]]
    colnames(X) <- paste(N, 1:ncol(X), sep = "")
    stores[rownames(X), colnames(X)] <- X
  }
  ## record the data set that each component comes from
  featgroups <- substring(colnames(stores), 1, -1 + nchar(colnames(stores))) #kevin look at later
  ## learn how to cross-predict component values between each pair of data sets
  componentModels <- lapply(names(object@data), function(N) {
    if(verbose) cat(N, "\n", file = stderr())
    X <- object@data[[N]]
    allNA <- apply(X, 2, function(xcol) all(is.na(xcol)))
    X <- X[, !allNA]
    ident <- apply(X, 1, function(x) length(unique(x)))
    X <- t(X[ident > 1, ])
    Xout <- object@outcome[rownames(X),]
    Xstores <- stores[rownames(X),]
    plsRegression <- lapply(names(object@data), function(M) {
      if(verbose) cat("\t", M, "\n", file = stderr())
      Y = Xstores[, featgroups == M]
      learn <- try(plsr(Y ~ X, 2))
      if (inherits(learn, "try-error")) {
        learn <- NA
        extend <- array(0, dim = c(nrow(Y), ncol(Y), 2))
        dimnames(extend) <- list(rownames(Y),
                                 paste(M, 1:ncol(Y), sep = ""),
                                 1:2)
      } else {
        extend <- predict(learn, X)
      }
      list(learn = learn, extend = extend)
    })
    names(plsRegression) <- names(object@data)
    slurp <- lapply(plsRegression, function(x) x$extend[,,2]) # subset the portion of the 3D array that corresponds to the final item of the output of the plsr function, which fits a model iteratively
    allPred <- do.call(cbind, slurp)
    colnames(allPred) <- NN
    list(plsRegression = plsRegression, allPred = allPred)
  })
  names(componentModels) <- names(object@data)
  ## Extract predicted values from all component datasets
  componentPredictions <- sapply(componentModels, function(x) x$allPred)
  ## Combine all the predictions into a 3D array
  myArray <- array(0, dim = c(nrow(stores), ncol(stores), length(componentPredictions)))
  dimnames(myArray) <- list(rownames(stores), colnames(stores),
                            names(object@data))
  for (I in 1:length(componentPredictions)) {
    B <- componentPredictions[[I]]
    myArray[rownames(B), colnames(B), I] <- B
  }
  ## Average the predictions to reduce to samples-x-components matrix
  meanPreds <- apply(myArray, 1:2, mean, na.rm = TRUE)
  list(meanPreds = meanPreds, compModels = componentModels)
}

setClass("plasmaPredictions",
         slots = c(meanPredictions = "matrix",
                   riskDF = "data.frame",
                   riskModel = "coxph",
                   splitModel = "coxph",
                   SF = "survfit"))

## plot method for plasmaPredictions object
setMethod("plot", c("plasmaPredictions", "missing"), function(x, y,  col = c("blue", "red"),
    lwd = 2, xlab = "", ylab = "Fraction Surviving", mark.time = TRUE, legloc = "topright", ...) {
  S <- summary(x@splitModel)
  PT <- S$sctest[3]
  plot(x@SF, col = col, lwd = lwd, xlab = xlab, ylab = ylab, mark.time = mark.time, ...)
  legend(legloc, paste(c("low", "high"), "risk"), col = col, lwd = lwd)
  title(sub = paste("p =", formatC(PT, format = "e", digits = 2)))
})


setClass("plasma",
         contains = "plasmaPredictions",
         slots = c(traindata = "MultiOmics",
                   compModels = "list",
                   fullModel = "coxph"))

## object = MultiOmics
## multi = MultiplePLSCoxModels
plasma <- function(object, multi) {
  temp <- extendCoxModels(object, multi)
  mp <- temp$meanPreds
  colms <- c(multi@timevar, multi@eventvar)
  riskDF <- data.frame(object@outcome[, colms], mp)
  riskDF <- na.omit(riskDF)         # see below
  form <- formula(paste("Surv(", multi@timevar, ",", multi@eventvar, "== \"",
                        multi@eventvalue, "\") ~ .", sep = ""))
  model <- coxph(form, data = riskDF)
  model1 <- step(model, trace = 0)  # could throw an error with NA's in data
  riskDF$Risk <- predict(model1)
  lh <- c("Low", "High")
  riskDF$Split <- factor(lh[1 + 1*(riskDF$Risk >= median(riskDF$Risk))],
                          levels = lh)
  form <- formula(paste("Surv(", multi@timevar, ",", multi@eventvar, "== \"",
                        multi@eventvalue, "\") ~ Risk", sep = ""))
  rmodel <- coxph(form, riskDF)
  form <- formula(paste("Surv(", multi@timevar, ",", multi@eventvar, "== \"",
                        multi@eventvalue, "\") ~ Split", sep = ""))
  smodel <- coxph(form, riskDF)
  SF <- survfit(form, riskDF)
  new("plasma",
      traindata = object,
      compModels = temp$compModels,
      fullModel = model1,
      meanPredictions = temp$meanPreds,
      riskDF = riskDF,
      riskModel = rmodel,
      splitModel = smodel,
      SF = SF)
}

## summary method for plasma objects
setMethod("summary", "plasma", function(object, ...) {
})

## predict method for plasma objects
setMethod("predict", "plasma", function(object, newdata = NULL,
                                        type = c("components", "risk", "split"),
                                        ...) {
  if (is.null(newdata)) {
    newdata <- object@traindata
  }
  componentModels <- object@compModels
  testData <- newdata@data
  tempor <- lapply(names(newdata@data), function(N) {
    wb <- componentModels[[N]]
    goForIt <- lapply(names(newdata@data), function(M) {
      wbp <- wb$plsRegression[[M]]
      localModel <- wbp$learn
      Z <- wbp$extend
      if (inherits(localModel, "mvr")) {
        localTest <- t(testData[[N]])[, dimnames(localModel$coefficients)[[1]]]
        predictions <- predict(localModel, localTest)
      } else {
        Y <- testData[[N]]
        predictions <- array(NA, dim = c(ncol(Y), ncol(Z), 2))
        dimnames(predictions) <- list(colnames(Y),
                                      paste(M, 1:ncol(Z), sep = ""),
                                      1:2)
      }
      predictions
    })
    names(goForIt) <- names(newdata@data)
    ## reintegrate
    lap <- lapply(names(goForIt), function(N) {
      G <- goForIt[[N]]
      top <- dim(G)[3]
      X <- G[,,top]
      if (!inherits(X, "array")) {
        X <- matrix(X, ncol = 1)
        colnames(X) <- paste(N, "1", sep = "")
      }
      X
    })
    do.call(cbind, lap)
  })
  names(tempor) <- names(newdata@data)

  barf <- predict(object@fullModel)
  testOut <- newdata@outcome
  meanPreds <- object@meanPredictions
  tstArray <- array(NA, dim = c(nrow(testOut), ncol(meanPreds), length(newdata@data)))
  dimnames(tstArray) <- list(rownames(testOut), colnames(meanPreds), names(newdata@data))
  for (I in 1:length(tempor)) {
    B <- tempor[[I]]
    if(any(rownames(B) != rownames(tstArray)) ) {
      stop ("bad rownames, loop:", I, "\n")
    }
    if(any(colnames(B) != colnames(tstArray)) ) {
      stop("bad colnames, loop:", I, "\n")
    }
    if (I > dim(tstArray)[3] ) {
      stop("bad thrid dimesnkon", I, "\n")
    }
    tstArray[rownames(B), colnames(B), I] <- B
  }
  meanTestPreds <- as.data.frame(apply(tstArray, 1:2, mean, na.rm = TRUE))
  all(!is.na(meanTestPreds))

  testcox <- predict(object@fullModel, meanTestPreds) # was mod1
  testOut$Risk <- testcox
  testOut$Split <- factor(c("low", "high")[1+1*(testcox > median(barf))],
                          levels = c("low", "high"))

  riskForm <- formula(object@riskModel)
  riskMod <- coxph(riskForm, data = testOut)

  splitForm <- formula(object@splitModel)
  splitMod <- coxph(splitForm, data = testOut)
  SF <- survfit(splitForm, testOut)

  new("plasmaPredictions",
      meanPredictions = as.matrix(meanTestPreds),
      riskDF = testOut,
      riskModel = riskMod,
      splitModel = splitMod,
      SF = SF)
})
