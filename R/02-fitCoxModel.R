setOldClass("plsRcoxmodel")
setOldClass("coxph")
setOldClass("survfit")

setClass("SingleModel",
         slots=c(plsmod = "plsRcoxmodel",
                 dsname = "character",
                 Xout = "data.frame",
                 SF = "survfit",
                 riskModel = "coxph",
                 splitModel = "coxph"))

fitSingleModel <- function(multi, N, timevar, eventvar, eventvalue) {
  ## get the training data
  X <- multi@data[[N]]
  allNA <- apply(X, 2, function(xcol) all(is.na(xcol)))
  X <- X[, !allNA]
  ident <- apply(X, 1, function(x) length(unique(x[!is.na(x)])))
  X <- X[ident > 1, ]
  nz <- apply(X, 1, function(x) sum(x > min(x, na.rm=TRUE), na.rm=TRUE))
  chooser <- nz > max(3, 0.01*nrow(X))
  if (sum(chooser) == 0) {
    chooser <-nz == max(nz, na.rm=TRUE)
  }
  X <- X[chooser,]
  out <-  multi@outcome
  Xout <-out[colnames(X),]
  mynt <- round(1 + log10(nrow(X)))
  ## Fit the PLS model
  tm <- Xout[, timevar]
  ev <- Xout[eventvar] == eventvalue
  plsmod <- plsRcoxmodel(t(X), time = tm, event = ev, nt = mynt)
  Xout$Risk = predict(plsmod, verbose = FALSE)
  Xout$Split <- 1*(Xout$Risk > median(Xout$Risk))
  riskModel <- coxph(formula(paste("Surv(", timevar, ",", eventvar, "== \"",
                                   eventvalue, "\") ~ Risk", sep = "")),
                     data = Xout)
  splitModel <- coxph(formula(paste("Surv(", timevar, ",", eventvar, "==\"",
                                    eventvalue, "\") ~ Split", sep = "")),
                      data = Xout)
  SF <- survfit(formula(paste("Surv(", timevar, ",", eventvar, "==\"",
                                    eventvalue, "\") ~ Split", sep = "")),
                      data = Xout)
  new("SingleModel",
      plsmod = plsmod,
      dsname = N,
      Xout = Xout,
      SF = SF,
      riskModel = riskModel,
      splitModel = splitModel)
}


## predict method for SingleModel objects
setMethod("predict", "SingleModel", function(object, newdata,
                                             type = c("components", "risk", "split", "survfit"),
                                             ...) {
  type <- match.arg(type)
  model <- switch(type,
                  components = object@plsmod,
                  risk = object@riskModel,
                  split = object@splitModel)
  if(missing(newdata)) {
    result <- predict(model, verbose = FALSE)
  } else {
    base <- predict(object@plsmod, verbose = FALSE)
    newd <- as.data.frame(t(newdata@data[[object@dsname]]))
    blankRow <- apply(is.na(newd), 1, all)
    newd <- newd[!blankRow, colnames(object@plsmod$dataX)]
    if (type != "components") {
      rsk <- predict(object@plsmod, newdata = newd, verbose = FALSE)
      newd$Risk <- rsk
    }
    if (type %in% c("split", "survfit")) {
      newd$Split <-  1*(newd$Risk > median(base))
    }
    if (type == "survfit") {
      simpleData <- data.frame(newdata@outcome[!blankRow,], Split = newd$Split)
      result <- survfit(formula(object@splitModel), data = simpleData)
    } else {
      result <- predict(model, newdata = newd, verbose = FALSE)
    }
  }
  result
})

## summary method for SingleModel objects
setMethod("summary", "SingleModel", function(object, ...) {
  cat("Risk Model:\n", file = stdout())
  print(summary(object@riskModel))
  cat("SplitModel:\n", file = stdout())
  print(summary(object@splitModel))
  cat("PLS Model:\n", file = stdout())
  S <- summary(object@plsmod$FinalModel)
  S$call <- NULL
  S
})

## plot method for SingleModel objects
setMethod("plot", c("SingleModel", "missing"), function(x, y,  col = c("blue", "red"), lwd = 2, xlab = "", ylab = "Fraction Surviving", mark.time = TRUE, legloc = "topright", ...) {
  plot(x@SF, col = col, lwd = lwd, xlab = xlab, ylab = ylab,  mark.time = mark.time, ...)
  legend(legloc, paste(c("Low", "High"), "Risk"), col = col, lwd = lwd)
  invisible(x)
})
