# A function to generate Sigma with a special structure 
# Input:
#    p: dimension of output matrix
#    rho: parameter to generate elements of Sigma (rho^(abs(i-j)))
#    correlation: "equal" means all the elements are equal except diagonal ones
#                 "decreasing" means all the elements equal to rho^|i-j|
#                  default is "decreasing"
# Output:
#    Sigma: A positive definite matrix with diagonal elements 1 and others less than 1

Sig <- function(p, rho, correlation = NULL){
  if(is.null(correlation)){
    correlation <- "decreasing"
  }
  if(correlation == "equal"){
    Sigma <- matrix(1, p, p)
    for (i in 1 : p) {
      for (j in 1 : p) {
        if(j != i){
          Sigma[i, j] <- rho
        }
      }
    }
  }
  if(correlation == "decreasing"){
    Sigma <- matrix(0, p, p)
    for (i in 1 : p) {
      for (j in 1 : p) {
        Sigma[i, j] <- rho^(abs(i - j))
      }
    }
  }
  
  return(Sigma)
}

