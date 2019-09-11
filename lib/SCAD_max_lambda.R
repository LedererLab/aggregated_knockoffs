# SCAD_max_lambda

SCAD_max_lambda <- function(X, y, nlambda = 500){
  if (!requireNamespace('ncvreg', quietly=T))
    stop('ncvreg is not installed', call.=F)
  out <- ncvreg(X, y, family = "gaussian", penalty = "SCAD", nlambda = nlambda)
  betaM <- out$beta[-1, ]
  lambda <- out$lambda    # a decreasing sequence

  max_lambda <- rep(0, ncol(X))
  for (i in 1 : ncol(X)) {
    nonzero <- which(betaM[i, ] != 0)
    if(length(nonzero) == 0){
      max_lambda[i] <- 0
    }else{
      max_lambda[i] <- lambda[nonzero[1]]
    }
  }
  return(max_lambda)
}

#--------------------------------------------------------------
# stat.SCAD_lambdasmax
stat.SCAD_lambdasmax <- function (X, Xk, y) {
  swap = rbinom(ncol(X), 1, 0.5)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap),
                  byrow = TRUE)
  X.swap = X * (1 - swap.M) + Xk * swap.M
  Xk.swap = X * swap.M + Xk * (1 - swap.M)
  Z = SCAD_max_lambda(cbind(X.swap, Xk.swap), y)
  p = ncol(X)
  orig = 1 : p
  W = pmax(Z[orig], Z[orig + p]) * sign(Z[orig] - Z[orig + p])
  W = W * (1 - 2 * swap)
  
  tem <- Z[which(swap == 1)]
  Z[which(swap == 1)] <- Z[which(swap == 1) + p]
  Z[which(swap == 1) + p] <- tem
  
  return(list(W = W, Z = Z))
}



