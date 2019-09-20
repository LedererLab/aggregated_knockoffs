# Code written by: Lun Li and Johannes Lederer
TransformData <- function(rawdata, method = 'LiLederer', g = 1) {
  # This function accepts a set of data and perform request transformation.
  # There are two options: our method of transformation (Lilederer), or Kurtz's
  # method (Kurtzetal).
  #
  # Args:
  #   rawdata: A matrix of data.
  #   method : Transformation with 2 options: 'LiLederer' and 'Kurtzetal'.
  #   g      : Integer indicates the number of clusters used for the GMM model.
  #            Only used for 'LiLederer'.
  #             
  # Returns  : A matrix of transformed data using the requested method.
  if (class(rawdata) != "matrix") {
    stop("invalid data type")
  }
  if (method == 'LiLederer') {
    return(centerByGMM(getClr(rawdata), g)) # our method
  } else if(method == 'Kurtzetal') {
    return(getClr(rawdata, FALSE))
  } else{
    stop("unknown method")
  }
}


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

# Code written by: Lun Li and Johannes Lederer
getClr <- function(comp.mtx, nonzero.only = TRUE, pc = 1e-7){
  # Perform a clr transformation to the given matrix with compositional data.
  # If the argument nonzero.only is set to TRUE, this function will only 
  # consider non-zeros when computing the geometric means.
  #
  # If the argument nonzero.only is set to FALSE, then this becomes a 
  # modification of Kurtz's way to generalize Aitchison's method. Kurtz, in his
  # paper, added a psedo count of one unit to the orininal cound data to avoid 
  # problems with zero counts. We here add a small pseudo count (pc) to the rows 
  # that contain zeros to avoid problems with zeros. The same value (pc) will be
  # subtracted from the computed means for those rows.
  #
  # Args:
  #   comp.mtx    : matrix with compositional data
  #   nonzero.only: whether only consider non-zero data are considere or all 
  #                  data default TRUE
  #   pc          : pseudo count used to handle zeros in Muller's way. only
  #                used when nonzero.only is set to FALSE. default = 1e-7.
  #
  # Returns: A matrix with clr-transformed data
  #
  # Reference: 
  # Kurtz, Zachary D., et al. "Sparse and Compositionally Robust Inference of 
  # Microbial Ecological Networks". Equition (1) on pg 5.
  if (is.null(dim(comp.mtx)) | class(comp.mtx[1]) != "numeric"|
      !is.logical(nonzero.only) | !is.numeric(pc) | pc <=0) {
    stop("invalid argument")
  }
  rslt <- matrix(nrow=nrow(comp.mtx), ncol=ncol(comp.mtx))
  for (i in 1:nrow(comp.mtx)) {
    x <- comp.mtx[i,] # for each row of the compositional data matrix
    if (nonzero.only) {
      gx <- exp(mean(log(x[x > 0])))
    } else {
      gx <- ifelse(!any(x == 0), exp(mean(log(x))), exp(mean(log(x + pc))) - pc)
    }
    rslt.temp <- log(x/gx)
    # set all the data that are originally 0 (-Inf after log) back to 0
    rslt.temp[rslt.temp == -Inf] <- 0
    rslt[i, ] <- rslt.temp
  }
  colnames(rslt) <- colnames(comp.mtx)
  return(rslt)
}


