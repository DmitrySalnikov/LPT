L1.stat <- function(x, y, normalizer) {
  sum(sapply(x, function(x.i) { log( 1 + abs(x.i - y) / normalizer ) } ))
}

L1.normalizer <- function(z, n.z) {
  normalizer <- 0
  for (i in 1:(n.z-1))
    for (j in (i+1):(n.z))
      normalizer <- normalizer + abs(z[i] - z[j])
  2 * normalizer / (n.z * (n.z - 1))
}

next.perm.idx <- function(perm.idx, k, n) {
  perm.idx
  for (i in k:1) {
    if (perm.idx[i] < n+i-k) {
      perm.idx[i] = perm.idx[i] + 1
      if (i != k) {
        for (j in (i+1):k) {
          perm.idx[j] = perm.idx[j-1] + 1
        }
      }
      return (perm.idx)
    }
  }
  NULL
}

generate.perm.idx <- function(k, n) {
  idx <- list(1:k)
  i <- 1
  while (TRUE) {
    next.idx <- next.perm.idx(idx[[i]], k, n)
    if (!is.null(next.idx)) {
      i <- i + 1
      idx[[i]] <- next.idx
    } else {
      break
    }
  }
  idx
}

exact.permutations <- function(x, y, n.x, n.y) {
  n <- if (n.x == n.y) n.x%/%2 else min(n.x, n.y)
  perms <- c(x, y)
  for (k in 1:n) {
    x.idx <- generate.perm.idx(k, n.x)
    y.idx <- generate.perm.idx(k, n.y)
    for (x.i in x.idx) {
      for (y.i in y.idx) {
        perms <- cbind(perms, c(x[-x.i], y[y.i], x[x.i], y[-y.i]))
      }
    }
  }
  perms
}

L1.test <- function(x, y, normalization = FALSE, n.permutations = 800, exact = FALSE, permutations = NULL, check_permutations = TRUE, ...) {
  n.x <- length(x)
  if(n.x < 2) stop("not enough 'x' observations")
  n.y <- length(y)
  if(n.y < 2) stop("not enough 'y' observations")

  z0 <- c(x, y)
  n.z <- n.x + n.y

  if (is.null(permutations)) {
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
    normalizer <- L1.normalizer(z0, n.z)
    if (normalizer == 0) {
      normalizer <- 1
      normalization <- FALSE
      warning("normalizer has calculated as 0, normalization is canceled")
    }
  } else {
    normalizer <- 1
  }


  stat <- L1.stat(x, y, normalizer)
  perm.stat <- apply(permutations, 2, function(z) { L1.stat(z[1:n.x], z[(n.x+1):n.z], normalizer) } )
  p.value <- mean(perm.stat > stat)

  data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative <- "two-sided"
  method <- "L1 two-sample permutation test"

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
