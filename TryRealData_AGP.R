#------------------------------------------------#
#---A real application on American Gut Project---#
#------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# load sources
source("lib/ImportData.R")
source("lib/TransformData.R")
source("lib/CLRTransformation.R")
source("lib/CenteringByGMM.R")
source("lib/SelectIndLinCols.R")
source("lib/AKO.R")
source("lib/SCAD_max_lambda.R")
source("lib/Lasso_max_lambda.R")

# libraries
library(pracma)
library(knockoff)
library(MASS)
library(LaplacesDemon)
library(ncvreg)
library(matrixcalc)  # to use the function is.positive.definite


# Import Phylum level data before January 2018
source("lib/DataUS012018.R")

# X and Y
X.raw <- as.matrix(X.raw.US)
X.trans <- TransformData(X.raw, method='LiLederer')
X <- selindlinearcols(X.trans)$full.new
Y <- Y.BMI.US

# n, p, mu, Sigma
n <- nrow(X)
p <- ncol(X)
mu <- colMeans(X)
Sigma <- (t(X) %*% X) / n
if(!is.positive.definite(Sigma)){
  stop("Sigma should be positive definite")
}

# ksteps, kc, fdr
ksteps <- 20
kc <- 20
fdr <- 0.2

# results collectors
as.KO.la <- matrix(0, p, length(fdr))   # KO: Barber & Candes (knockoff filter); la: Lasso
as.AKO.la <- matrix(0, p, length(fdr)) # AKO: aggregated FDR control; la: Lasso

as.KO.SCAD <- matrix(0, p, length(fdr))   # KO: Barber & Candes (knockoff filter)
as.AKO.SCAD <- matrix(0, p, length(fdr)) # AKO: aggregated FDR control

# methods
tic()
#----------Lasso-----------
W <- matrix(0, p, ksteps)
set.seed(999)
for (w in 1 : ksteps) {
  Xk <- create.gaussian(X = X, mu = mu, Sigma = Sigma, method = "sdp")
  stat <- stat.Lasso_lambdasmax(X, Xk, Y)
  W[ , w] <- stat$W
}
v <- 1                    # counter for length(fdr)
for(j in fdr){
  t <- knockoff.threshold(W = W[ , 1], fdr = j * (1 - 1 / 2^kc), offset = 0)
  as.KO.la[ , v] <- as.numeric(W[ , 1] >= t)
  
  as.AKO.la[ , v] <- AKO(W = W, ksteps = ksteps, fdr = j, offset = 0, r = 0)
  v <- v + 1
}
#----------SCAD-----------
W <- matrix(0, p, ksteps)
set.seed(999)
for (w in 1 : ksteps) {
  Xk <- create.gaussian(X = X, mu = mu, Sigma = Sigma, method = "sdp")
  stat <- stat.SCAD_lambdasmax(X, Xk, Y)
  W[ , w] <- stat$W
}
v <- 1                    # counter for length(fdr)
for(j in fdr){
  t <- knockoff.threshold(W = W[ , 1], fdr = j * (1 - 1 / 2^kc), offset = 0)
  as.KO.SCAD[ , v] <- as.numeric(W[ , 1] >= t)

  as.AKO.SCAD[ , v] <- AKO(W = W, ksteps = ksteps, fdr = j, offset = 0, r = 0)
  v <- v + 1
}
toc()

#---display the results---#
colnames(X)[which(as.KO.la[,1] != 0)]
colnames(X)[which(as.AKO.la[,1] != 0)]
colnames(X)[which(as.KO.SCAD[,1] != 0)]
colnames(X)[which(as.AKO.SCAD[,1] != 0)]









