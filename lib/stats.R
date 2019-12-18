#! R
# Filename: stats.R

# required libraries
library(glmnet)
library(ncvreg)

#--lasso_max_lambda--#
Lasso_max_lambda <- function(X, y, nlambda = 500, intercept = T, standardize = T, ...) {
  if (!requireNamespace('glmnet', quietly = T))
    stop('glmnet is not installed', call. = F)
  
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
  
  fit <- glmnet::glmnet(X, y, lambda = lambda, intercept = intercept, standardize = F, standardize.response = F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices])
}

#--SCAD_max_lambda--#
SCAD_max_lambda <- function(X, y, nlambda = 500){
  if (!requireNamespace('ncvreg', quietly = T))
    stop('ncvreg is not installed', call. = F)
  
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

#--logistic regression estimator--#
logReg_est <- function(X, y){
  if (!requireNamespace('glmnet', quietly = T))
    stop('glmnet is not installed', call. =  F)
  
  cvfit <- cv.glmnet(X, y, family = "binomial")
  lamb_cv <- cvfit$lambda.min
  fit <- glmnet(X, y, family = "binomial", lambda = lamb_cv)
  beta <- fit$beta
}

logReg_max_lambda <- function(X, y, nlambda = 500) {
  if (!requireNamespace('glmnet', quietly = T))
    stop('glmnet is not installed', call. = F)
  
  n = nrow(X); p = ncol(X)
  fit <- glmnet(X, y, family = "binomial", nlambda = nlambda)
  betaM <- fit$beta
  lambda <- fit$lambda    # a decreasing sequence
  
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

#--statistic after pluging in knockoffs--#
stat.ToDoFDRCon <- function (X, Xk, y, stat.type = "max_lambda_lasso") {
  p = ncol(X)
  orig = 1 : p
  swap = rbinom(p, 1, 0.5)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap),
                  byrow = TRUE)
  X.swap = X * (1 - swap.M) + Xk * swap.M
  Xk.swap = X * swap.M + Xk * (1 - swap.M)
  
  if(stat.type == "max_lambda_scad"){
    Z = SCAD_max_lambda(cbind(X.swap, Xk.swap), y)
    W = pmax(Z[orig], Z[orig + p]) * sign(Z[orig] - Z[orig + p])
  }else if(stat.type == "max_lambda_lasso"){
    Z = Lasso_max_lambda(cbind(X.swap, Xk.swap), y)
    W = pmax(Z[orig], Z[orig + p]) * sign(Z[orig] - Z[orig + p])
  }else if(stat.type == "estimator_logistic"){
    Z = logReg_est(cbind(X.swap, Xk.swap), y)
    W = pmax(abs(Z[orig]), abs(Z[orig + p])) * sign(abs(Z[orig]) - abs(Z[orig + p]))
  }else if(stat.type == "max_lambda_logistic"){
    Z = logReg_max_lambda(cbind(X.swap, Xk.swap), y)
    W = pmax(Z[orig], Z[orig + p]) * sign(Z[orig] - Z[orig + p])
  }else{
    stop("stat.type should be one of 'max_lambda_lasso', 'max_lambda_scad', 
         'estimator_logistic' and 'max_lambda_logistic'")
  }
  
  W = W * (1 - 2 * swap)
  tem <- Z[which(swap == 1)]
  Z[which(swap == 1)] <- Z[which(swap == 1) + p]
  Z[which(swap == 1) + p] <- tem
  
  return(list(W = W, Z = Z))
}