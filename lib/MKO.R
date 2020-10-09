#! R
# Filename: MultipleKnockoffs.R

#---MKO: multiple knockoffs (Holden and Helton 2018)---#
# Args:
#    W: a matrix whose length of row is the number of variables (p) and length of column is the number of knockoffs
#    fdr: a value between 0 and 1
#    offset：0 or 1. 0：knockoff; 1: knockoff+
#
# Return: p-dimensional vector whose elements take value 0 or 1

MKO <- function(W, fdr = 0.2, offset = 0){
  # Step 1: using ksteps times standard knockoff with corresponding FDRs
  p <- nrow(W)
  k <- ncol(W)
  if(k <= 1){
    stop("The number of knockoffs (ncol(W)) should be larger than 1 for MKO")
  }
  t <- MKO.threshold(W, fdr, offset)
  hatS <- as.numeric(W[, 1] >= t)
  return(hatS)
}

#---MKO.threshold: threshold of multiple knockoffs (Holden and Helton 2018)---#
# Args:
#    W: a matrix whose length of row is the number of variables (p) and length of column is the number of knockoffs
#    fdr: a value between 0 and 1
#    offset：0 or 1. 0：knockoff; 1: knockoff+
#
# Return: the threshold of MKO

MKO.threshold <- function (W, fdr = 0.1, offset = 1) 
{
  if (offset != 1 && offset != 0) {
    stop("Input offset must be either 0 or 1")
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t) (offset + sum(W[, -1] >= t))/((ncol(W) - 1) * max(1, 
                                                                                   sum(W[, 1] >= t))))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}


