scaleZ <- function(Z, freq1 = NULL) {
  #Scale Z according to van Raden 2008
  Z0 <- Z
  if ( is.null(freq1) ) freq1 <- colMeans(Z0)/2
  sum2pq <- sum(2*freq1*(1-freq1))
  for (i in 1:nrow(Z)) {
    Z[i,] <- (Z0[i,] - 2 * freq1)/sqrt(sum2pq)
  }
  return(Z)
}