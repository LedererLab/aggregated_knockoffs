# Code written by: Lun Li and Johannes Lederer
library(mclust)
centerByGMM <- function(x, g){
  # Perform a gmm clustering using function 'Mclust' from 'mclust' package
  # on the data passed in, then re-center the data from each cluster and 
  # return the result
  #
  # Args:
  #   x: a matrix of data
  #   g: number of clusters. if > 1: apply the gmm clustering with g clusters.
  #      Default = 1.
  #
  # Returns: A matrix of data undergoes the recentering process
  if (is.null(dim(x)) | class(x[1]) != "numeric" | class(g) != "numeric") {
    stop("argument is not valid")
  }
  if(g == 1) { # for g = 1, do not apply mclust functions
    # for(i in 1:dim(x)[2]) {
    #   foo <- x[which(x[,i]!=0),i]
    #   x[which(x[,i]!=0),i] <- foo - mean(foo)
    # }
    # return(x)
    return(apply(x, 2, function(y) y - mean(y)))
  }
  rslt <- x
  for (i in 1:ncol(x)) {
    if(sum(x[,i] != 0) < 10) next # skip columns with non-zero entries < 10
    mc.temp <- Mclust(data.frame(x[,i]), G = g)
    for(j in 1:g) {
      foo <-  x[mc.temp$classification == j, i]
      rslt[mc.temp$classification == j, i] <- foo - mean(foo)
    }
  }
  colnames(rslt) <- colnames(x)
  return(rslt)
}