L1.test <- function(x, y, normalization = FALSE, n.permutations = 10000, exact = NULL, permutations = NULL, check_permutations = TRUE, ...) {
  L.test(x, y, 1, normalization, n.permutations, exact, permutations, check_permutations)
}

L2.test <- function(x, y, normalization = FALSE, n.permutations = 10000, exact = NULL, permutations = NULL, check_permutations = TRUE, ...) {
  L.test(x, y, 2, normalization, n.permutations, exact, permutations, check_permutations)
}
