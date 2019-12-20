L.test <- function(test.stat, method_prefix, x, y, normalization = FALSE, n.permutations = 10000, exact = NULL, permutations = NULL, check_permutations = TRUE) {
  n.x <- length(x)
  if(n.x < 2) stop("not enough 'x' observations")
  n.y <- length(y)
  if(n.y < 2) stop("not enough 'y' observations")

  z0 <- c(x, y)
  n.z <- n.x + n.y

  if (is.null(permutations)) {
    if (is.null(exact)) {
      exact = if (n.exact.perms(n.x, n.y, n.z) <= n.permutations) TRUE else FALSE
    }

    if (exact) {
      permutations <- exact.permutations(x, y, n.x, n.y)
      n.permutations = ncol(permutations)
    } else {
      if (n.permutations < 1) stop("permutation set can not be generated")
      permutations <- replicate(n.permutations, sample(z0))
    }
  } else {
    perm.dim <- dim(permutations)
    if (is.null(perm.dim) || length(perm.dim) != 2 || perm.dim[1] != n.z) stop("permutations size is not equal to pulled sample size")
    if (check_permutations) {
      if (sum(apply(permutations, 2, function(z) { sort(z) != sort(z0) } )) != 0) stop("permutations contain values not coincided with the given samples")
    }
    n.permutations = ncol(permutations)
    exact = NA
  }

  if (normalization) {
    normalizer <- L.normalizer(z0, n.z)
    if (normalizer == 0) {
      normalizer <- 1
      normalization <- FALSE
      warning("normalizer has calculated as 0, normalization is canceled")
    }
  } else {
    normalizer <- 1
  }

  stat <- test.stat(x, y, normalizer)
  perm.stat <- apply(permutations, 2, function(z) { test.stat(z[1:n.x], z[(n.x+1):n.z], normalizer) } )
  p.value <- mean(perm.stat > stat)

  data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative <- "two-sided"
  method <- paste(method_prefix, "two-sample permutation test")

  rval <- list(
    data.name = data.name,
    alternative = alternative,
    method = method,
    normalization = normalization,
    normalizer = normalizer,
    n.x = n.x,
    n.y = n.y,
    n.permutations = n.permutations,
    exact = exact,
    p.value = p.value
  )
  class(rval) <- "htest"
  return(rval)
}