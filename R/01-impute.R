## Routines for imputation
##
## We assume X is a numeric matrix, with rows = features and columns = samples.
## So, imputation happens one row at a time. 

## Here we want to replace missing data by either the mean or the mode of
## the observed values. We decide which one to use based on how many distinct
## values there are. We arbitrasriyl say that anything with five or fewer
## distinct values is likely to be a caregorical variable, and so use the mode.
## Otherwise, we use the mean.
meanModeImputer <- function(X) {
  ## row-by-row
  for (i in 1L:nrow(X)) {
    if (sum(is.na(X[i,])) > 0){ # Only update a row if there is missing data
      ux <- unique(na.omit(X[i,]))
      if (length(ux) < 6) {
        mymode <- ux[which.max(tabulate(match(X[i,], ux)))]
        X[i, is.na(X[i,])] <- mymode
      } else {
        mymean <- mean(X[i,], na.rm = TRUE) 
        X[i, is.na(X[i,])] <- mymean
      }
    }
  }
  return(X)
}

## We wnat to reaplce missing data by actual observed values, sampled according
## to their empirical distribution.
samplingImputer <- function(X) {
  randomSampler <- function(V) { # V is a vector
    V[is.na(V)] <- sample(na.omit(V), sum(is.na(V)), replace = TRUE)
    V
  }
  t(apply(X, 1, randomSampler))
}
