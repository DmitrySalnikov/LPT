L.normalizer <- function(z, n.z) {
  normalizer <- 0
  for (i in 1:(n.z-1))
    for (j in (i+1):(n.z))
      normalizer <- normalizer + abs(z[i] - z[j])
  2 * normalizer / (n.z * (n.z - 1))
}
