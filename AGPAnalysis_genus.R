#------------------------------------------------#
#---A real application on American Gut Project---#
#------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load sources
source("lib/ImportData_genus.R")
source("lib/SelectIndLinCols.R")
source("lib/AggregatedKnockoffs.R")
source("lib/stats.R")
# source("lib/clr.R")

# libraries
library(pracma)
library(knockoff)
library(matrixcalc)  # to use the function is.positive.definite

# Import Genus level data before January 2018
source("lib/DataUS012018_genus.R")

# X and Y
X.raw <- as.matrix(X.raw.US)
Y <- Y.BMI.US

# log-transformation
X.new <- selindlinearcols(X.raw, tol = 5e-2)$full.new
X.new <- t(apply(X.new, 1, function(x) x/sum(x)))
X <- log(X.new)

# # CLR
# X <- clr.FJ(X.raw)
# X <- selindlinearcols(X, tol = 5e-2)$full.new

# Transform Y into 0, 1 for logistic regression
Y[which(Y < 30)] <- 0
Y[which(Y >= 30)] <- 1


# n, p, mu, Sigma
n <- nrow(X)
p <- ncol(X)
mu <- colMeans(X)
Sigma <- cov(X)
if(!is.positive.definite(Sigma)){
  stop("Sigma should be positive definite")
}

# ksteps, fdr
ksteps <- 5
fdr <- 0.1

# number of repeations
numRep <- 20

# results collectors
as.KO <- matrix(0, p, numRep)   # KO: Barber & Candes (knockoff filter); la: Lasso
as.AKO.m <- matrix(0, p, numRep) # modified AKO: modified aggregated FDR control

# methods
tic()
for(rep in 1 : numRep){
  #----------Lasso-----------
  W <- matrix(0, p, ksteps)
  for (w in 1 : ksteps) {
    Xk <- create.gaussian(X = X, mu = mu, Sigma = Sigma, method = "sdp")
    stat <- stat.ToDoFDRCon(X, Xk, Y, stat.type = "max_lambda_logistic")
    W[ , w] <- stat$W
  }
  # standard knockoff
  t <- knockoff.threshold(W = W[ , 1], fdr = fdr, offset = 0)
  as.KO[ , rep] <- as.numeric(W[ , 1] >= t)
  # modified AKO
  as.AKO.m[ , rep] <- AKO.m(W = W, ksteps = ksteps, fdr = fdr, offset = 0, r = 0)
}
toc()

#---display the results---#
cat("n = ", n, "\n")
cat("p = ", p, "\n")
cat("fdr = ", fdr, "\n")
colnames(X)[which(apply(as.KO, 1, sum) / numRep >= 0.5)]
colnames(X)[which(apply(as.AKO.m, 1, sum) / numRep >= 0.5)]



