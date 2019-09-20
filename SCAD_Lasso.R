# simulation examples for displaying the relationship between FDR and TargetFDR, selection accuracy and TargetFDR
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("lib/AggregatedKnockoffs.R")
source("lib/myPlot.R")
source("lib/Sigma.R")
source("lib/SCAD_max_lambda.R")
source("lib/Lasso_max_lambda.R")

# loading libraries
library(pracma)
library(knockoff)
library(MASS)
library(LaplacesDemon)
library(ncvreg)

# generating data
# initial parameters
n <- 200                                 # sample size
p <- 100                                 # number of variables
s0 <- 20                                 # number of nonzero variables
rho <- 0                                 # correlation parameter in Sigma
sigma <- 1                               # standard deviation of noise
trueS <- sample(1 : p, s0)
beta0 <- rep(0, p)                       # true parameter
beta0[trueS] <- 1
# Sigma
Sigma = Sig(p, rho, correlation = "decreasing")
# noise
noise <- rnorm(n, 0, sigma)
SNR <- as.numeric((t(beta0) %*% Sigma %*% beta0) / sigma^2)
snr <- 5
# scaled beta0
beta0 <- beta0 * sqrt(snr / SNR)
# X, y
X <- scale(mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
y <- as.vector(X %*% beta0 + noise)
# calculate mu and Sigma based on X
mu <- colMeans(X)
Sigma <- (t(X) %*% X) / n

ksteps = 5
kc = 5
numFDR = 100
fdr = seq(0.001, 1, length.out = 100) / (1 - 1 / 2^kc)


tic()
lf <- length(fdr)
a <- matrix(0, lf, numFDR)
pwr.SCAD <- FDP.SCAD <- a # "SCAD": SCAD-KO
Pwr.SCAD <- FDR.SCAD <- rep(0, lf)
pwr.SCAD.p <- FDP.SCAD.p <- a # "SCAD.p": SCAD-KO+
Pwr.SCAD.p <- FDR.SCAD.p <- rep(0, lf)
pwr.SCAD.2s <- FDP.SCAD.2s <- a # "SCAD.2s": SCAD-AKO
Pwr.SCAD.2s <- FDR.SCAD.2s <- rep(0, lf)
pwr.SCAD.2s.p <- FDP.SCAD.2s.p <- a # "SCAD.2s.p": SCAD-AKO+
Pwr.SCAD.2s.p <- FDR.SCAD.2s.p <- rep(0, lf)

pwr.Lasso <- FDP.Lasso <- a # "Lasso": Lasso-KO
Pwr.Lasso <- FDR.Lasso <- rep(0, lf)
pwr.Lasso.p <- FDP.Lasso.p <- a # "Lasso.p": Lasso-KO+
Pwr.Lasso.p <- FDR.Lasso.p <- rep(0, lf)
pwr.Lasso.2s <- FDP.Lasso.2s <- a # "Lasso.2s": Lasso-AKO
Pwr.Lasso.2s <- FDR.Lasso.2s <- rep(0, lf)
pwr.Lasso.2s.p <- FDP.Lasso.2s.p <- a # "Lasso.2s.p": Lasso-AKO+
Pwr.Lasso.2s.p <- FDR.Lasso.2s.p <- rep(0, lf)

trueS <- as.numeric(beta0 > 0)

for(i in 1 : numFDR){
  #----------Lasso-----------#
  W <- matrix(0, p, ksteps)
  for (w in 1 : ksteps) {
    Xk <- create.gaussian(X = X, mu = mu, Sigma = Sigma, method = "sdp")
    stat <- stat.Lasso_lambdasmax(X, Xk, y)
    W[ , w] <- stat$W
  }
  offset <- 0               # knockoff
  v <- 1                    # counter for length(fdr)
  for(j in fdr){
    t <- knockoff.threshold(W = W[ , 1], fdr = j * (1 - 1 / 2^kc), offset = offset)
    as <- which(W[ , 1] >= t)
    FDP.Lasso[v, i] <- (sum(beta0[as] == 0))/max(length(as), 1)
    pwr.Lasso[v, i] <- (sum(beta0[as] != 0))/s0
    
    as.Lasso.2s <- which(AKO(W = W, ksteps = ksteps, fdr = j, offset = offset, r = 0) == 1)
    FDP.Lasso.2s[v, i] <- (sum(beta0[as.Lasso.2s] == 0))/max(length(as.Lasso.2s), 1)
    pwr.Lasso.2s[v, i] <- (sum(beta0[as.Lasso.2s] != 0))/s0
    v <- v + 1
  }
  offset <- 1               # knockoff+
  v <- 1                    # counter for length(fdr)
  for(j in fdr){
    t <- knockoff.threshold(W = W[ , 1], fdr = j * (1 - 1 / 2^kc), offset = offset)
    as <- which(W[ , 1] >= t)
    FDP.Lasso.p[v, i] <- (sum(beta0[as] == 0))/max(length(as), 1)
    pwr.Lasso.p[v, i] <- (sum(beta0[as] != 0))/s0
    
    as.2s <- which(AKO(W = W, ksteps = ksteps, fdr = j, offset = offset, r = 0) == 1)
    FDP.Lasso.2s.p[v, i] <- (sum(beta0[as.2s] == 0))/max(length(as.2s), 1)
    pwr.Lasso.2s.p[v, i] <- (sum(beta0[as.2s] != 0))/s0
    v <- v + 1
  }
  
  #----------SCAD-----------#
  W <- matrix(0, p, ksteps)
  for (w in 1 : ksteps) {
    Xk <- create.gaussian(X = X, mu = mu, Sigma = Sigma, method = "sdp")
    stat <- stat.SCAD_lambdasmax(X, Xk, y)
    W[ , w] <- stat$W
  }
  offset <- 0               # knockoff
  v <- 1                    # counter for length(fdr)
  for(j in fdr){
    t <- knockoff.threshold(W = W[ , 1], fdr = j * (1 - 1 / 2^kc), offset = offset)
    as <- which(W[ , 1] >= t)
    FDP.SCAD[v, i] <- (sum(beta0[as] == 0))/max(length(as), 1)
    pwr.SCAD[v, i] <- (sum(beta0[as] != 0))/s0
    
    as.2s <- which(AKO(W = W, ksteps = ksteps, fdr = j, offset = offset, r = 0) == 1)
    FDP.SCAD.2s[v, i] <- (sum(beta0[as.2s] == 0))/max(length(as.2s), 1)
    pwr.SCAD.2s[v, i] <- (sum(beta0[as.2s] != 0))/s0
    v <- v + 1
  }
  offset <- 1               # knockoff+
  v <- 1                    # counter for length(fdr)
  for(j in fdr){
    t <- knockoff.threshold(W = W[ , 1], fdr = j * (1 - 1 / 2^kc), offset = offset)
    as <- which(W[ , 1] >= t)
    FDP.SCAD.p[v, i] <- (sum(beta0[as] == 0))/max(length(as), 1)
    pwr.SCAD.p[v, i] <- (sum(beta0[as] != 0))/s0
    
    as.2s <- which(AKO(W = W, ksteps = ksteps, fdr = j, offset = offset, r = 0) == 1)
    FDP.SCAD.2s.p[v, i] <- (sum(beta0[as.2s] == 0))/max(length(as.2s), 1)
    pwr.SCAD.2s.p[v, i] <- (sum(beta0[as.2s] != 0))/s0
    v <- v + 1
  }
}
toc()

#-------------Lasso------------#
FDR.Lasso <- apply(FDP.Lasso, 1, mean)
Pwr.Lasso <- apply(pwr.Lasso, 1, mean)
SA.Lasso <- - FDR.Lasso + Pwr.Lasso

FDR.Lasso.2s <- apply(FDP.Lasso.2s, 1, mean)
Pwr.Lasso.2s <- apply(pwr.Lasso.2s, 1, mean)
SA.Lasso.2s <- - FDR.Lasso.2s + Pwr.Lasso.2s

FDR.Lasso.p <- apply(FDP.Lasso.p, 1, mean)
Pwr.Lasso.p <- apply(pwr.Lasso.p, 1, mean)
SA.Lasso.p <- - FDR.Lasso.p + Pwr.Lasso.p

FDR.Lasso.2s.p <- apply(FDP.Lasso.2s.p, 1, mean)
Pwr.Lasso.2s.p <- apply(pwr.Lasso.2s.p, 1, mean)
SA.Lasso.2s.p <- - FDR.Lasso.2s.p + Pwr.Lasso.2s.p

#------------SCAD-----------#
FDR.SCAD <- apply(FDP.SCAD, 1, mean)
Pwr.SCAD <- apply(pwr.SCAD, 1, mean)
SA.SCAD <- - FDR.SCAD + Pwr.SCAD       # SA: selection accuracy

FDR.SCAD.p <- apply(FDP.SCAD.p, 1, mean)
Pwr.SCAD.p <- apply(pwr.SCAD.p, 1, mean)
SA.SCAD.p <- - FDR.SCAD.p + Pwr.SCAD.p

FDR.SCAD.2s <- apply(FDP.SCAD.2s, 1, mean)
Pwr.SCAD.2s <- apply(pwr.SCAD.2s, 1, mean)
SA.SCAD.2s <- - FDR.SCAD.2s + Pwr.SCAD.2s

FDR.SCAD.2s.p <- apply(FDP.SCAD.2s.p, 1, mean)
Pwr.SCAD.2s.p <- apply(pwr.SCAD.2s.p, 1, mean)
SA.SCAD.2s.p <- - FDR.SCAD.2s.p + Pwr.SCAD.2s.p


#--------------------------------------------------------------------#
#---plot the curves of FDR and selection accuracy w.r.t Target FDR---#
#--------------------------------------------------------------------#
#----Lasso----#
#-------------#
# define file names and legend names
filename1 = "Results/Lasso-SA-n=200-p=100-snr=5.pdf"  
filename2 = "Results/Lasso-FDR-n=200-p=100-snr=5.pdf"     
filename3 = "Results/Lasso+-SA-n=200-p=100-snr=5.pdf"  
filename4 = "Results/Lasso+-FDR-n=200-p=100-snr=5.pdf"      
legendnames = c("lasso-KO", "lasso-AKO")
legendnames.p = c("lasso-KO+", "lasso-AKO+")

#-------------Lasso-knockoff-------------#
FDR <- cbind(FDR.Lasso, FDR.Lasso.2s)
Pwr <- cbind(Pwr.Lasso, Pwr.Lasso.2s)
SA <- cbind(SA.Lasso, SA.Lasso.2s)

# display the relationship between selection accuracy and TargetFDR 
pdf(file = filename1, width = 5, height = 5)
myPlot.SA(x = fdr * (1 - 1 / 2^kc), y = SA, main = "lasso", legendnames = legendnames)
dev.off()


# display the relationship between FDR and TargetFDR
pdf(file = filename2, width = 5, height = 5)
myPlot.FDR(x = fdr * (1 - 1 / 2^kc), y = FDR, main = "lasso",legendnames = legendnames)
dev.off()

#----------Lasso-knockoff+-----------#
FDR <- cbind(FDR.Lasso.p, FDR.Lasso.2s.p)
Pwr <- cbind(Pwr.Lasso.p, Pwr.Lasso.2s.p)
SA <- cbind(SA.Lasso.p, SA.Lasso.2s.p)

# display the relationship between selection accuracy and TargetFDR 
pdf(file = filename3, width = 5, height = 5)
myPlot.SA(x = fdr * (1 - 1 / 2^kc), y = SA, main = "lasso+",legendnames = legendnames.p)
dev.off()

# display the relationship between FDR and TargetFDR
pdf(file = filename4, width = 5, height = 5)
myPlot.FDR(x = fdr * (1 - 1 / 2^kc), y = FDR, main = "lasso+", legendnames = legendnames.p)
dev.off()

#------------#
#----SCAD----#
#------------#
# define file names and legend names
filename1 = "Results/SCAD-SA-n=200-p=100-snr=5.pdf"  
filename2 = "Results/SCAD-FDR-n=200-p=100-snr=5.pdf"       
filename3 = "Results/SCAD+-SA-n=200-p=100-snr=5.pdf"  
filename4 = "Results/SCAD+-FDR-n=200-p=100-snr=5.pdf"     
legendnames = c("scad-KO", "scad-AKO")
legendnames.p = c("scad-KO+", "scad-AKO+")

#-----------SCAD-knockoff-----------#
FDR <- cbind(FDR.SCAD, FDR.SCAD.2s)
Pwr <- cbind(Pwr.SCAD, Pwr.SCAD.2s)
SA <- cbind(SA.SCAD, SA.SCAD.2s)     # SA: selection accuracy

# display the relationship between selection accuracy and TargetFDR 
pdf(file = filename1, width = 5, height = 5)
myPlot.SA(x = fdr * (1 - 1 / 2^kc), y = SA, main = "scad", legendnames = legendnames)
dev.off()

# display the relationship between FDR and TargetFDR
pdf(file = filename2, width = 5, height = 5)
myPlot.FDR(x = fdr * (1 - 1 / 2^kc), y = FDR, main = "scad", legendnames = legendnames)
dev.off()

#-------------SCAD-knockoff+-------------#
FDR <- cbind(FDR.SCAD.p, FDR.SCAD.2s.p)
Pwr <- cbind(Pwr.SCAD.p, Pwr.SCAD.2s.p)
SA <- cbind(SA.SCAD.p, SA.SCAD.2s.p)

# display the relationship between selection accuracy and TargetFDR 
pdf(file = filename3, width = 5, height = 5)
myPlot.SA(x = fdr * (1 - 1 / 2^kc), y = SA, main = "scad+", legendnames = legendnames.p)
dev.off()

# display the relationship between FDR and TargetFDR
pdf(file = filename4, width = 5, height = 5)
myPlot.FDR(x = fdr * (1 - 1 / 2^kc), y = FDR, main = "scad+", legendnames = legendnames.p)
dev.off()





