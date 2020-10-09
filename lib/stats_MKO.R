#! R
# Filename: stats_MKO.R

stat.MKO.W <- function(Z_2k){
  # Z_2k is a p*2k matrix, whose first column is the variable Z based on X, and other columns are the variables based on knockoffs \tilde{X}_i.
  p <- nrow(Z_2k)
  k <- ceiling(ncol(Z_2k)/2)
  W <- matrix(0, p, k)
  for(j in 1 : p){
    for (i in 1:k) {
      W[j, i] <- Z_2k[j, i] - sum(Z_2k[j, (k+2) : (2*k)])/(k-1)
    }
  }
  return(W)
}
  