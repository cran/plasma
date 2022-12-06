plasmaEnv <- new.env()

loadESCAdata <- function(env = plasmaEnv) {
  U <- url("http://silicovore.com/data/TCGA-ESCA.RData")
  on.exit(close(U))
  if (!inherits(U, "connection")) {
    stop("Could not create connection.")
  }
  created <- load(U, env)
  expect <- c("assemble", "m450info", "Outcome")
  if (!(all(expect %in% created))) {
    stop("Could not load all objects.")
  }
  if (!(all(created %in% expect))) {
    stop("Somehow loaded unexpected objects.")
  }
  invisible(created)
}

loadLUSCdata <- function(env = plasmaEnv) {
  U <- url("http://silicovore.com/data/TCGA-LUSC1.RData")
  on.exit(close(U))
  if (!inherits(U, "connection")) {
    stop("Could not create connection.")
  }
  created <- load(U, env)
  expect <- c("assemble", "m450info", "Outcome")
  if (!(all(expect %in% created))) {
    stop("Could not load all objects.")
  }
  if (!(all(created %in% expect))) {
    stop("Somehow loaded unexpected objects.")
  }
  invisible(created)
}
