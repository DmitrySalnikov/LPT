source("L.normalizer.R")
source("exact.permutations.R")
source("L.test.R")

L1.stat <- function(x, y, normalizer) {
  sum(sapply(x, function(x.i) { log( 1 + abs(x.i - y) / normalizer ) } ))
}

L1.test <- function(x, y, normalization = FALSE, n.permutations = 10000, exact = NULL, permutations = NULL, check_permutations = TRUE, ...) {
  L.test(L1.stat, "L1", x, y, normalization, n.permutations, exact, permutations, check_permutations)
}