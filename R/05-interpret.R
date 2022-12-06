setClass("Contribution",
         slots = c("contrib" = "matrix",
                   "datasets" = "character"))

## summary method for MultiOmics objects
setMethod("summary", "Contribution", function(object, ...) {
  cat("Contributions of dataset", object@datasets[1], "to components from ",
      file = stdout())
  if(length(object@datasets) > 1) {
    cat(object@datasets[2], ".\n", file = stdout())
  } else {
    cat("all datasets.\n", file = stdout())
  }
  summary(object@contrib)
})

setMethod("[", "Contribution", function(x, i, j,  ..., drop = FALSE) {
  contrib <- x@contrib
  if (!missing(i)) {
    contrib <- x@contrib[i,]
  }
  if (!missing(j)) {
    contrib <- contrib[, j]
  }
  new("Contribution", datasets = x@datasets, contrib = contrib)
})


setMethod("image", c("Contribution"), function(x, col = viridis(64), mai = c(1.82, 1.52, 0.32, 0.32), ...) {
  NR <- nrow(x@contrib)
  NC <- ncol(x@contrib)
  M <- max(abs(x@contrib)) + 0.001
  opar <- par(mai = mai)
  on.exit(par(opar))
  image(1:NR, 1:NC, x@contrib, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
        zlim = c(-M, M),
        col = col, ...)
  mtext(rownames(x@contrib), side = 1, line = 1, at = 1:NR, las=2)
  mtext(colnames(x@contrib), side = 2, line = 1, at = 1:NC, las=2)
})

if (!isGeneric("heat"))
  setGeneric("heat",
             function(object, ...) { standardGeneric("heat") }
             )

# might want to set 'margins=c(5,15)' in ... for long feature names
setMethod("heat", "Contribution", function(object, main = "Contributions",
                                           col = viridis(64),
                                           mai = c(1.52, 0.32, 0.82, 1.82),
                                           ...) {
  M <- max(abs(object@contrib)) + 0.001
  if (main == "Contributions") main = paste(object@datasets, collapse = " => ")
  opar <- par(mai = mai)
  on.exit(par(opar))
  heatmap(object@contrib, scale = "none",
            main = main, col = col, zlim = c(-M, M), ...)
})

## object is a thing of the "plasma" class
## M and N are names of data sets being modeled
getCompositeWeights <- function(object, N, M) {
  cm <- object@compModels # returns a list of lists for all pairs of data sets
  wb <- cm[[N]]$plsRegression # ge the model for data set M
  inside <- wb[[M]]
  learn <- inside$learn # rest will FAIL if we were unable to construct a PLS model
  V <- as.vector(L <- learn$loadings)
  Y <- learn$Yloadings
  cross <- L %*% t(Y)
  new("Contribution",
      contrib = cross,
      datasets = c(N, M))
}

getAllWeights <- function(object, N) {
  whatever <- lapply(names(object@compModels), function(D) {
    getCompositeWeights(object, N, D)@contrib
  })
  W <- which(sapply(whatever, ncol) == 1)
  names(whatever) <- names(object@compModels)
  if (length(W) == 1) {
    for (I in W) {
      colnames(whatever[[I]]) <- paste(names(whatever)[I], 1, sep = "")
    }
  }
  cont <- do.call(cbind, whatever)
  new("Contribution",
      contrib = cont,
      datasets = N)
}

getTop <- function(object, N = 1) {
  topFeatures <- apply(object@contrib, 2, function(x, N) {
    mx <- max(abs(x))
    sorted <- rank(mx - abs(x))
    S <- sorted %in% 1:N
    rownames(object@contrib)[which(S)][order(sorted[S])]
  }, N = N)
  topFeatures
}


pickSignificant <- function(object, alpha) {
  Q <- max(abs(quantile(object@contrib, c(alpha, 1 - alpha))))
  sig <- apply(abs(object@contrib) > Q, 1, any)
  object[sig,]
}

influencer <- function(object) {
  toto <- sapply(names(object@compModels), function(N) getAllWeights(object, N))
  wolf <- lapply(toto, function(X) X@contrib)
  do.call(rbind, wolf)
}
