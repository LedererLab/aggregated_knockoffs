#! R
# Filename: dataGen.R

# load libraries
library(MASS)

#---data generated function for linear Gaussian model---#
linGauData <- function(n, p, s, sigma, rho, snr = 5){
  CovMat = matrix(0, p, p)
  for (i in 1 : p) {
    for (j in 1 : p) {
      CovMat[i, j] = rho^abs(i-j)
    }
  }
  mean = rep(0, p)
  ## 1: X ##
  X = scale(mvrnorm(n, mu = mean, Sigma = CovMat))
  ## 2: noise ##
  noise = rnorm(n, 0, sigma)
  ## 3: beta0 ##
  beta0 = rep(0, p)
  beta0[sample(1 : p, s)] = 1
  SNR = sum((X %*% beta0)^2) / (n * var(noise))
  beta0 = beta0 * sqrt(snr / SNR)
  ## 4: y ##
  y = X %*% beta0 + noise
  # return X, y, beta0
  return(list(X = X, y = y, beta0 = beta0))
}

#---data generated function for logistic regression---#
logRegData <- function(n, p, s, sigma, rho, snr = 5){
  CovMat = matrix(0, p, p)
  for (i in 1 : p) {
    for (j in 1 : p) {
      CovMat[i, j] = rho^abs(i-j)
    }
  }
  mean = rep(0, p)
  ## 1: X ##
  X = scale(mvrnorm(n, mu = mean, Sigma = CovMat))
  ## 2: noise ##
  noise = rnorm(n, 0, sigma)
  ## 3: beta0 ##
  beta0 = rep(0, p)
  beta0[sample(1 : p, s)] = 1
  SNR = sum((X %*% beta0)^2) / (n * var(noise))
  beta0 = beta0 * sqrt(snr / SNR)
  ## 4: y ##
  y = rep(0, n)
  for (i in 1 : n) {
    prob = exp(X[i, ] %*% beta0 + noise[i]) / (1 + exp(X[i, ] %*% beta0 + noise[i]))
    y[i] = rbinom(1, size = 1, p = prob)
  }
  # return X, y, beta0
  return(list(X = X, y = y, beta0 = beta0))
}


#---Generate the compositional data from the abundance data---#
# 1. generete abundance data W from a log-normal dististrion
# 2. transform abundance data W into compositional data Z through Z_ij = log(W_ij/sum(W_i))
# 3. generate noise
# 4. generate parameter beta
# 5. generate outcome Y by Y_i = prod(Z_i, beta) + noise
logComposData <- function(n, p, s, sigma, rho, snr = 5, data.type = "linear"){
  CovMat = matrix(0, p, p)
  for (i in 1 : p) {
    for (j in 1 : p) {
      CovMat[i, j] = rho^abs(i-j)
    }
  }
  mean = rep(0, p)
  ## 1 ##
  abundanceMat = exp(mvrnorm(n, mu = mean, Sigma = CovMat))
  ## 2 ##
  logComposData = scale(t(apply(abundanceMat, 1, function(x) log(x/sum(x)))))
  
  ## 3 ##
  noise = rnorm(n, 0, sigma)
  if(data.type == "linear"){
    ## 4 ##
    beta0 = rep(0, p)
    beta0[sample(1 : p, s)] = 1
    SNR = sum((logComposData %*% beta0)^2) / (n * var(noise))
    beta0 = beta0 * sqrt(snr / SNR)
    ## 5 ##
    outcome = logComposData %*% beta0 + noise
  }else if(data.type == "logistic"){
    ## 4 ##
    beta0 = rep(0, p)
    beta0[sample(1 : p, s)] = rnorm(s, 0, 1)
    SNR = sum((logComposData %*% beta0)^2) / (n * var(noise))
    beta0 = beta0 * sqrt(snr / SNR)
    ## 5 ##
    outcome = rep(0, n)
    for (i in 1 : n) {
      prob = exp(logComposData[i, ] %*% beta0 + noise[i]) / (1 + exp(logComposData[i, ] %*% beta0 + noise[i]))
      outcome[i] = rbinom(1, size = 1, p = prob)
    }
  }else{
    stop("data.type should be 'linear' or 'logistic'")
  }
  # return beta, logComposData, outcome
  return(list(logComposData = logComposData, outcome = outcome, beta0 = beta0))
}




