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

#--------------------------------------------------------------
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



