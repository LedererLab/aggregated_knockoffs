#! R
# Filename: Simulation1.R
# File discription: This file is to display the results of aggregation knockoffs 
#                   in the linear regression case. 

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load sources
source("lib/AggregatedKnockoffs.R")
source("lib/dataGen.R")
source("lib/stats.R")

# loading libraries
library(pracma)
library(knockoff)

# generating data
# initial parameters
n <- 300 # sample size
p <- 100 # number of variables
s <- 30 # number of nonzero variables
sigma <- 1 # standard deviation of noise
rho <- 0.5 # correlation parameter in Sigma
snr <- 5 # signal to noise ratio
 
# number for Aggregation Knockoffs
ksteps = 5
# repeating number to get the average of FDR
numRep = 50
# fdr sequence
fdr = seq(0.001, 1, length.out = 100)


#---generate linear Gaussian data---#
Data = linGauData(n, p, s, sigma, rho, snr)
X = Data$X
y = Data$y
beta0 = Data$beta0


# calculate mu and Sigma based on X
mu <- colMeans(X)
Sigma <- cov(X)


# matries to collect the results
a <- matrix(0, length(fdr), numRep)
pwr <- FDP <- a # KO
pwr.AKO.m <- FDP.AKO.m <- a # modified AKO


#--repeation start--#
tic()
for(i in 1 : numRep){
  W <- matrix(0, p, ksteps)
  for (w in 1 : ksteps) {
    Xk <- create.gaussian(X = X, mu = mu, Sigma = Sigma, method = "sdp")
    stat <- stat.ToDoFDRCon(X, Xk, y, stat.type = "max_lambda_lasso")
    W[ , w] <- stat$W
  }           
  v <- 1 # counter for length(fdr)
  for(j in fdr){
    #--KO--#
    t <- knockoff.threshold(W = W[ , 1], fdr = j, offset = 0)
    as <- which(W[ , 1] >= t)
    FDP[v, i] <- (sum(beta0[as] == 0))/max(length(as), 1)
    pwr[v, i] <- (sum(beta0[as] != 0))/s
    #--modified AKO--#
    as.AKO.m <- which(AKO.m(W = W, ksteps = ksteps, fdr = j, offset = 0, r = 0) == 1)
    FDP.AKO.m[v, i] <- (sum(beta0[as.AKO.m] == 0)) / max(length(as.AKO.m), 1)
    pwr.AKO.m[v, i] <- (sum(beta0[as.AKO.m] != 0)) / s
    
    v <- v + 1
  }
}
toc()
#--repeation end--#


# results:
#-----FDR, Power-----#
FDR <- apply(FDP, 1, mean)
Pwr <- apply(pwr, 1, mean)

FDR.AKO.m <- apply(FDP.AKO.m, 1, mean)
Pwr.AKO.m <- apply(pwr.AKO.m, 1, mean)


# results display
source("lib/myPlot.R")
# source("lib/plotlinear.R")





