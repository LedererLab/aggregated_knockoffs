#---AKO: aggregation knockoffs---#
# Args:
#    W: a matrix
#    ksteps: a number larger than 1
#    fdr: a value between 0 and 1
#    offset：0 or 1. 0：knockoff; 1: knockoff+
#    r：a value in [0, 1). r = 0 if taking union; r in (0,1) if taking fraction
#
# Return: p-dimensional vector whose elements take value 0 or 1

AKO <- function(W, ksteps = 5, fdr = 0.2, offset = 0, r = 1){
  # Step 1: using ksteps times standard knockoff with corresponding FDRs
  p <- nrow(W)
  if(ksteps <= 1){
    stop("ksteps should be larger than 1 for AKO")
  }
  hatS.all <- matrix(0, p, ksteps)
  for(i in 1 : ksteps){
    # Threshold
    t <- knockoff.threshold(W[ , i], fdr = fdr / (2^i), offset = offset)
    # active set 
    as <- which(W[ , i] >= t)
    hatS.all[as, i] <- 1
  }
  S <- hatS.all
  
  # Step 2: aggregate the k estimated active sets by taking the union
  if(length(which(S != 0 & S != 1)) >0){
    stop("S is a matrix taking value only 0 and 1 in which 1 indicates the position of nonzero component.")
  }
  ksteps <- ncol(S)
  hatS <- rep(0, nrow(S))
  if(ksteps == 1){
    N <- S
  }else{
    N <- apply(S, 1, sum)
  }
  if(r != 0){
    hatS[which(N >= ksteps * r)] <- 1 # criteria----partial Union sets
  }else{
    hatS[which(N > 0)] <- 1 # criteria----full Union sets
  }
  return(hatS)
}

#--------------------------------------------------------------
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


#--------------------------------------------------------------
# modify the function stat.glmnet_lambdasmax a little (return Z)
# Lasso_max_lambda
Lasso_max_lambda <- function(X, y, nlambda=500, intercept=T, standardize=T, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  # Standardize the variables
  if( standardize ){
    X = scale(X)
  }
  
  n = nrow(X); p = ncol(X)
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  
  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept,
                        standardize=F, standardize.response=F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices])
}

# Lasso_max_lambda <- function(X, y, nlambda = 500, intercept=T, standardize=T, ...){
#   if (!requireNamespace('glmnet', quietly=T))
#     stop('glmnet is not installed', call.=F)
#   # Standardize the variables
#   if( standardize ){
#     X = scale(X)
#   }
# 
#   lambda_max = max(abs(t(X) %*% y)) / n
#   lambda_min = lambda_max / 2e3
#   k = (0:(nlambda-1)) / nlambda
#   lambda = lambda_max * (lambda_min/lambda_max)^k
# 
#   out <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept,
#                         standardize=F, standardize.response=F, ...)
#   betaM <- out$beta
# 
#   max_lambda <- rep(0, ncol(X))
#   for (i in 1 : ncol(X)) {
#     nonzero <- which(betaM[i, ] != 0)
#     if(length(nonzero) == 0){
#       max_lambda[i] <- 0
#     }else{
#       max_lambda[i] <- lambda[nonzero[1]]
#     }
#   }
#   return(max_lambda)
# }


# stat.Lasso_lambdasmax
stat.Lasso_lambdasmax <- function (X, Xk, y) {
  swap = rbinom(ncol(X), 1, 0.1)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap),
                  byrow = TRUE)
  X.swap = X * (1 - swap.M) + Xk * swap.M
  Xk.swap = X * swap.M + Xk * (1 - swap.M)
  Z = Lasso_max_lambda(cbind(X.swap, Xk.swap), y)
  p = ncol(X)
  
  # W <- rep(0, p)
  # for(orig in 1 : p){
  #   if(sign(Z[orig] - Z[orig + p]) == 0){
  #     W[orig] <- Z[orig]
  #   }else{
  #     W[orig] <- pmax(Z[orig], Z[orig + p]) * sign(Z[orig] - Z[orig + p])
  #   }
  # }
  
  orig <- 1 : p
  W = pmax(Z[orig], Z[orig + p]) * sign(Z[orig] - Z[orig + p])
  W = W * (1 - 2 * swap)
  
  tem <- Z[which(swap == 1)]
  Z[which(swap == 1)] <- Z[which(swap == 1) + p]
  Z[which(swap == 1) + p] <- tem
  
  return(list(W = W, Z = Z))
}








