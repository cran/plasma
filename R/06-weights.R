##

setClass("CombinedWeights",
         slots = c(combined = "matrix",
                   featureSize = "numeric",
                   dataSource = "factor"))

combineAllWeights <- function(pl) {
  combined <- lapply(names(pl@compModels), function(N) getAllWeights(pl, N)@contrib)
  names(combined) <- names(pl@compModels)
  featsize <- sapply(combined, nrow)
  datasrc <- factor(rep(names(combined), times = featsize))
  combined <- do.call(rbind, combined)
  new("CombinedWeights",
      combined = combined,
      featureSize = featsize,
      dataSource = datasrc)
}

## summary method for CombinedWeights objects
setMethod("summary", "CombinedWeights", function(object, ...) {
  contra <- object@combined
  temp <- list(object@dataSource)
  aggie <- function(X, L, F) {
    mu <- aggregate(X, L, F)
    rownames(mu) <- mu[, 1]
    mu$Group.1 <- NULL
    mu <- as.matrix(mu)
  }

  mean <- aggie(contra, temp, mean)
  sd <- aggie(contra, temp, sd)
  med <- aggie(contra, temp, median)
  mad <- aggie(contra, temp, mad)

  list(Mean = mean, SD = sd, Median = med, MAD = mad)
})

setMethod("image", "CombinedWeights", function(x, ...) {
  summer <- summary(x)
  mx <- max(abs(summer$Mean))
  sx <- max(summer$SD)
  dx <- max(abs(summer$Median))
  ax <- max(summer$MAD)
  N <- ncol(summer$Mean)
  M <- nrow(summer$Mean)
  opar <- par(mfrow = c(2,2))
  on.exit(par(opar))
  image(1:N, 1:M, t(summer$Mean), xlab = "Component", ylab = "Dataset",
        main = "Mean", zlim = c(-mx, mx), col = redgreen(64))
  image(1:N, 1:M, t(summer$SD), xlab = "Component", ylab = "Dataset",
        main = "SD", zlim = c(0, sx), col = bluescale(64))
  image(1:N, 1:M, t(summer$Median), xlab = "Component", ylab = "Dataset",
        main = "Median", zlim =c(-dx, dx), col = redgreen(64))
  image(1:N, 1:M, t(summer$MAD), xlab = "Component", ylab = "Dataset",
        main = "Mad", zlim = c(0, ax),col = bluescale(64))
  invisible(x)
})

stdize <- function(object, type = c("standard", "robust")) {
  type <- match.arg(type)
  summer <- summary(object)
  mu <- switch(type,
               standard = summer$Mean,
               robust = summer$Median)
  sigma <- switch(type,
                  standard = summer$SD,
                  robust = summer$MAD)
  M <- apply(mu, 2, function(X) rep(X, times = object@featureSize))
  S <- apply(sigma, 2, function(X) rep(X, times = object@featureSize))
  (object@combined - M)/S
}

interpret <- function(object, component, alpha = 0.05) {
  brute <- stdize(object)
  Q <- qnorm(1 - alpha/2) # two-sided 5% cutoff
  topA <- which(abs(brute[, component]) > Q)
  data.frame(Feature = rownames(brute)[topA],
             Source = object@dataSource[topA],
             Weight = brute[topA, component])
}

# merge feature contributiosn across all components
getFinalWeights <- function(object) {
  fm <- object@fullModel
  coy <- matrix(fm$coefficients, ncol = 1)
  mainterms <- attr(terms(fm), "term.labels")
  combined <- lapply(names(object@compModels), function(N) {
    getAllWeights(object, N)@contrib[, mainterms, drop = FALSE]
  })
  names(combined) <- names(object@compModels)
  runthrough <- lapply(names(combined), function(N) {
    X <- combined[[N]] %*% coy
    data.frame(Weight = X, Source = N, Feature = rownames(X))
  })
  FW <- do.call(rbind, runthrough)
  mu <- aggregate(FW$Weight, list(FW$Source), mean)
  mu2 <- rep(mu$x, times = as.vector(table(FW$Source)))
  sigma <- aggregate(FW$Weight, list(FW$Source), sd)
  sigma2 <- rep(sigma$x, times = as.vector(table(FW$Source)))
  FW$Standard <- (FW$Weight - mu2)/sigma2
  FW
}

setMethod("barplot", c("plasma"), function(height, source, n,
                                           direction = c("both", "up","down"),
                                           lhcol = c("cyan", "red"),
                                           ...) {
  direction <- match.arg(direction)
  wws <- getFinalWeights(height) # default name from function, but a plasma object

  ## Get positive and negative weights. Start with negative.
  wmutDn <- wws[wws$Source  == source, ]
  wmutDn <- wmutDn[order(wmutDn$Weight),]
  wmutDn <- wmutDn[wmutDn$Weight < 0,]
  m <- min(nrow(wmutDn), n)
  wmutDn <- wmutDn[1:m,]

  ## Then add positive weights.
  wmutUp <- wws[wws$Source  == source,]
  wmutUp <- wmutUp[order(wmutUp$Weight, decreasing = TRUE),]
  wmutUp <- wmutUp[wmutUp$Weight > 0,]
  m <- min(nrow(wmutUp), n)
  wmutUp <- wmutUp[1:m,]

  ##
  wmut <- switch(direction,
                 up = wmutUp,
                 down = wmutDn,
                 both = rbind(wmutUp, wmutDn))
  wmut <- wmut[order(abs(wmut$Weight), decreasing = FALSE),]
  wmut$Feature <- factor(as.character(wmut$Feature),
                         levels = unique(wmut$Feature))
##  print(summary(wmut))
  MX <- max(abs(wmut$Weight), na.rm = TRUE)
  ptr <- switch(direction,
                both = painter(c(-MX, MX), c(lhcol[1], "white", lhcol[2])),
                up = painter(c(0, MX), c("white", lhcol[2])),
                down = painter(c(-MX, 0), c(lhcol[1], "white")))
  mycol <- ptr(wmut$Weight)
  spar <- par(mai = c(0.82, 2.02, 0.42, 0.42))
  on.exit(par(spar))
  pts <- barplot(wmut$Weight, horiz = TRUE, col = mycol, xlim = c(-MX, 2*MX))
  mtext(rownames(wmut), side = 2, at = pts, las = 2, line = 1)
  par(new = TRUE)
  pin <- par("pin")
  mai <- par("mai")
  plt <- par("plt")
  weird <- c(mai[1] + 1/8*pin[2],   # bottom
             mai[2] + 13/16*pin[1], # left
             mai[3] + 1/2*pin[2],   # top
             mai[4] + 2/16*pin[1])  # right
  opar <- par(mai = weird)
  on.exit(par(opar))
  S <- switch(direction,
              both = seq(-MX, MX, length = 64),
              up = seq(0, MX, length = 64),
              down = seq(-MX, 0, length = 64))
  mat <- matrix(S, nrow = 1)
  image(1, S, mat, col = ptr(mat),
        xaxt = "n", xlab = "", ylab = "Weight", main = source)
})

painter <- function(range, colors, N = 64) {
##  cat("painter called with", range, "\n", file = stderr())
  crp <- colorRampPalette(colors)
  pal <- crp(N)
  bks <- seq(min(range), max(range), length = N)
  function(X) {
    idx <- cut(X, breaks = bks, labels = FALSE,
               include.lowest = TRUE)
    pal[idx]
  }
}
