L1.stat <- function(x, y, normalizer) {
  sum(sapply(x, function(x.i) { log( 1 + abs(x.i - y) / normalizer ) } ))
}

L1.test <- function(x, y, normalization = FALSE, n.permutations = 800, exact = FALSE) {
  n.x <- length(x)
  n.y <- length(y)
  
  normalizer <- 1
  # for (i in 1:(2*n-1))
  #   for (j in (i+1):(2*n))
  #     A <- A + abs(Z[i]-Z[j])
  # A <- A/(n*(2*n-1))
  
  stat <- L1.stat(x, y, normalizer)
  
  z0 <- c(x, y)
  # z.perm <- if (exact) exact.permutations(x, y) else replicate(n.permutations, sample(c(x, y)))
  z.perm <- replicate(n.permutations, sample(z0))
  perm.stat <- apply(z.perm, 2, function(z) { L1.stat(z[1:n.x], z[(n.x+1):(n.x+n.y)], normalizer) } )
  p.value <- mean(perm.stat > stat)
  
  data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative <- "two-sided"
  method <- "L1 two-sample permutation test"
  normalizer <- if(normalization) normalizer else normalization
  n.permutations <- if(exact) "exact" else n.permutations
  
  rval <- list(
    data.name = data.name,
    alternative = alternative,
    method = method, 
    normalizer = normalizer,
    n.x = n.x,
    n.y = n.y,
    n.permutations = n.permutations,
    p.value = p.value
  )
  class(rval) <- "htest"
  return(rval)
}
