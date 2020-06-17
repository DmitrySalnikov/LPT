L.stat <- function(x, y, gamma, normalizer) {
  sum(sapply(x, function(x.i) { log( 1 + (abs(x.i - y) / normalizer)**gamma ) } ))
}

L.test <- function(x, y, gamma, normalization = FALSE, n.permutations = 10000, exact = NULL, permutations = NULL, check_permutations = TRUE) {
  n.x <- length(x)
  if(n.x < 2) stop("not enough 'x' observations")
  n.y <- length(y)
  if(n.y < 2) stop("not enough 'y' observations")
  if (gamma <= 0) stop("gamma parameter must be strictly positive (gamma > 0)")
  
  z0 <- c(x, y)
  n.z <- n.x + n.y

  if (is.null(permutations)) {
    if (is.null(exact)) {
      n.exact <- n.exact.perms(n.x, n.y)
      exact = if (is.nan(n.exact) || n.exact > n.permutations) FALSE else TRUE
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
    exact = NULL
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

  stat <- L.stat(x, y, gamma, normalizer)
  perm.stat <- apply(permutations, 2, function(z) { L.stat(z[1:n.x], z[(n.x+1):n.z], gamma, normalizer) } )
  p.value <- mean(perm.stat > stat)

  data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative <- "two-sided"
  method <- paste0("L", gamma, " two-sample permutation test")
  permutations <- cbind(z0, permutations)
  stat.values <- c(z0 = stat, perm.stat)

  rval <- list(
    data.name = data.name,
    alternative = alternative,
    method = method,
    normalization = normalization,
    normalizer = normalizer,
    n.x = n.x,
    n.y = n.y,
    n.permutations = n.permutations,
    permutations = permutations,
    stat.values = stat.values,
    exact = exact,
    p.value = p.value
  )
  class(rval) <- "htest"
  return(rval)
}
