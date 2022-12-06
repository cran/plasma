## Need a class to hold the training and test data

## A multiomics data set. The data is a list, with each entry a data set
## from a different assay. Like MOFA, should have the same samples
setClass("MultiOmics",
         slots = c(data = "list",
                   outcome = "data.frame"))
## As in MOFA, each data set must contain the same set of samples
validMultiOmics <- function(object) {
  sampleSizes <- sapply(object@data, ncol)
  okSS <- all(sampleSizes == nrow(object@outcome))
  if (okSS) {
    namesOK <- sapply(object@data, function(DS) {
      all(colnames(DS) == rownames(object@outcome))
    })
    valid <- ifelse(all(namesOK), TRUE, "Sample names in all datasets must agree.")
  } else {
    valid <- "All datasets must have the same number of samples."
  }
  valid
}
setValidity("MultiOmics", validMultiOmics)

## summary method for MultiOmics objects
setMethod("summary", "MultiOmics", function(object, ...) {
  cat("Datasets:\n", file = stdout())
  print(sapply(object@data, dim))
  cat("Outcomes:\n", file = stdout())
  summary(object@outcome)
})

## plot method for MultiOmics objects.
setMethod("plot", c("MultiOmics", "missing"), function(x, y, ...) {
  binmat <- 1*sapply(x@data, function(DS) {
    useless <- apply(DS, 2, function(samp) {
      all(is.na(samp))
    })
    !useless
  })
  NF <- sapply(x@data, nrow)
  memberPlot(t(binmat), features = NF, ylab = "",  ...)
})

setMethod("[", "MultiOmics", function(x, i, j,  ..., drop = FALSE) {
  if (!missing(i)) {
    dataslice <- x@data[i]
  } else {
    dataslice <- x@data
  }
  if (!missing(j)) {
    dataslice <- lapply(dataslice, function(X) X[,j])
  }
  new("MultiOmics", data = dataslice, outcome = x@outcome[j, ])
})


## Here we assume that we have a list of datasets (on patient subsets)
## and a data.frame of outcomes thatn includes all patients.
## We assume that outcome is organized as patients x variables, but
## each dataset is organized as features x patients.
##
## IMPORTANT: Calling plsRcox will convert column names to syntactically
## valid ones. That means some of the RPPA names (referring to proteins
## like 14-3-3 wo't match correctly unless you convert them yourself
## beforehand. The safest way (and easiewst on users) is to make sure
## when creating this object that you have converted them.
prepareMultiOmics <- function(datalist, outcome) {
  ready <- lapply(as.list(names(datalist)), function(N) {
    DS <- datalist[[N]]
    ## remove odd punctuation from feature names.
    rownames(DS) <- make.names(rownames(DS))
    ## merge each dataset with the first two columns.
    filled <- merge(outcome[, 1:2], t(DS), by = "row.names", all.x = TRUE)
    ## first column should be the row names (patients)
    rownames(filled) <- filled[, 1]
    ## match the order
    filled <- filled[rownames(outcome),]
    ## columns two and three are copies of the outcome data
    filled <- filled[, -(1:3)]
    if (all(is.na(filled))) {
      stop("No matching sample names in dataset '", N, "'.\n", sep = "")
      }
    t(filled)
  })
  names(ready) <- names(datalist)
  new("MultiOmics", data = ready, outcome = outcome)
}


