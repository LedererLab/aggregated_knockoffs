# Centered log-ratio transformation for compostional data
# Codes written by Fang and Johannes

# Input: x is a matrix with compostional inputs
#
# Output: A centered log-ratio transformed matrix

clr.FJ <- function(X, pc = 1e-7){
  for (i in 1 : nrow(X)) {
    if(sum(X[i, ] == 0) > 0){
      gx = exp(mean(log(X[i, ] + pc))) - pc
      X[i, ] = log(X[i, ] / gx)
    }else{
      gx = exp(mean(log(X[i, ])))
      X[i, ] = log(X[i, ] / gx)
    }
  }
  return(X)
}

