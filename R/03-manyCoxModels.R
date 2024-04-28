setClass("MultiplePLSCoxModels",
         slots = c(models = "list",
                   timevar = "character",
                   eventvar = "character",
                   eventvalue = "character"))

## As in MOFA, each data set must contain the same set of samples
validMultipleCoxModels <- function(object) {
  types <- all(sapply(object@models, function(MO) (
    inherits(MO, "SingleModel"))
  ))
  types
}
setValidity("MultiplePLSCoxModels", validMultipleCoxModels)

## Need to define a class for the return value of this function
fitCoxModels <- function(multi, timevar, eventvar, eventvalue, verbose = TRUE) {
  firstPass <- lapply(names(multi@data), function(N) {
    if (verbose) cat("Fitting model with ", N, "\n", file = stderr())
    fitSingleModel(multi, N, timevar, eventvar, eventvalue)
  })
  names(firstPass) <- names(multi@data)
  new("MultiplePLSCoxModels",
      models = firstPass,
      timevar = timevar,
      eventvar = eventvar,
      eventvalue = eventvalue)
}

## summary method for MultiplePLSCoxModels objects
setMethod("summary", "MultiplePLSCoxModels", function(object, ...) {
  cat("An object containing MultiplePLSCoxModels based on:\n",
      file = stdout())
  print(names(object@models))
})

## plot method for MultiplePLSCoxModels objects
setMethod("plot", c("MultiplePLSCoxModels", "missing"), function(x, y,  col = c("blue", "red"), lwd = 2, xlab = "", ylab = "Fraction Surviving", mark.time = TRUE, legloc = "topright", ...) {
  L <- length(x@models)
  W <- ifelse(L > 12, 3, 2)
  H <- (L %/% W) + 1*(L %% W > 0)
  opar <- par(mfrow = c(H, W))
  on.exit(par(opar))
  for (N in names(x@models)) {
    plot(x@models[[N]], main = N,  col = col, lwd = lwd, xlab = xlab, ylab = ylab,
         mark.time = mark.time, legloc = legloc, ...)
    S <- summary(x@models[[N]]@splitModel)
    PT <- S$sctest[3]
    title(sub = paste("p =", formatC(PT, format = "e", digits = 2)))
  }
  invisible(x)
})

## predict method for MultiplePLSCoxModels objects
setMethod("predict", "MultiplePLSCoxModels", function(object, newdata,
                                             type = c("components", "risk", "split", "survfit"),
                                             ...) {
  type <- match.arg(type)
  if (missing(newdata)) {
    result <- lapply(object@models, predict, type = type, ...)
  } else {
    result <- lapply(object@models, predict, newdata = newdata, type = type, ...)
  }
  result
})

## getSizes function cutpasted from 02-fitCoxModel.R
getSizes <- function(object) {
  NT <- sapply(object@models, function(result) {
  S <- summary(result@plsmod$FinalModel)
  PT <- S$sctest[3]
  names(PT)<-"p.value"
  c(NT = result@plsmod$nt, cNT = result@plsmod$computed_nt, PT)
  })
  colnames(NT) <- names(object@models)
  t(NT)
 }

