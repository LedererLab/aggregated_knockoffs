# Code written by: Fang Xie and Johannes Lederer
# select only independently linear columns of a matrix
selindlinearcols <- function(mat, tol = 1e-4){ #1e-4
  ColNames <- colnames(mat)
  tt <- qr(mat, tol = tol)
  full.new <- mat[, tt$pivot[seq_len(tt$rank)]]
  names(full.new) <- ColNames[tt$pivot[seq_len(tt$rank)]]
  
  rem.mat <- mat[, -tt$pivot[seq_len(tt$rank)]]
  names(rem.mat) <- ColNames[-tt$pivot[seq_len(tt$rank)]]
  return(list(full.new = full.new, rem.mat = rem.mat))
}
