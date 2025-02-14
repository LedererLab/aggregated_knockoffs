#! R
# Filename: AggregatedKnockoffs.R

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
  # calculate fdr
  fdr = fdr / (1 - 1 / 2^ksteps)
  for(i in 1 : ksteps){
    # Threshold
    t <- knockoff.threshold(W[ , i], fdr = fdr / (2^i), offset = offset)
    # active set 
    as <- which(W[ , i] >= t)
    hatS.all[as, i] <- 1
  }
  S <- hatS.all
  
  # Step 2: aggregate the k estimated active sets by taking the union
  hatS <- rep(0, p)
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


#---A modified AKO function: in the first step replace fdr/2^i by fdr/2^(i-1)---#
# Args:
#    W: a matrix
#    ksteps: a number larger than 1
#    fdr: a value between 0 and 1
#    offset：0 or 1. 0：knockoff; 1: knockoff+
#    r：a value in [0, 1). r = 0 if taking union; r in (0,1) if taking fraction
#
# Return: p-dimensional vector whose elements take value 0 or 1

AKO.m <- function(W, ksteps = 5, fdr = 0.2, offset = 0, r = 1){
  # Step 1: using ksteps times standard knockoff with corresponding FDRs
  p <- nrow(W)
  if(ksteps <= 1){
    stop("ksteps should be larger than 1 for AKO")
  }
  hatS.all <- matrix(0, p, ksteps)
  for(i in 1 : ksteps){
    # Threshold
    t <- knockoff.threshold(W[ , i], fdr = fdr / (2^(i-1)), offset = offset)
    # active set 
    as <- which(W[ , i] >= t)
    hatS.all[as, i] <- 1
  }
  S <- hatS.all
  
  # Step 2: aggregate the k estimated active sets by taking the union
  hatS <- rep(0, p)
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


# q_i=q/k
AKO.aver <- function(W, ksteps = 5, fdr = 0.2, offset = 0, r = 1){
  # Step 1: using ksteps times standard knockoff with corresponding FDRs
  p <- nrow(W)
  if(ksteps <= 1){
    stop("ksteps should be larger than 1 for AKO")
  }
  hatS.all <- matrix(0, p, ksteps)
  
  for(i in 1 : ksteps){
    # Threshold
    t <- knockoff.threshold(W[ , i], fdr = fdr / ksteps, offset = offset)
    # active set 
    as <- which(W[ , i] >= t)
    hatS.all[as, i] <- 1
  }
  S <- hatS.all
  
  # Step 2: aggregate the k estimated active sets by taking the union
  hatS <- rep(0, p)
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

# q_i=2*q/k
AKO.aver.m <- function(W, ksteps = 5, fdr = 0.2, offset = 0, r = 1){
  # Step 1: using ksteps times standard knockoff with corresponding FDRs
  p <- nrow(W)
  if(ksteps <= 1){
    stop("ksteps should be larger than 1 for AKO")
  }
  hatS.all <- matrix(0, p, ksteps)
  
  for(i in 1 : ksteps){
    # Threshold
    t <- knockoff.threshold(W[ , i], fdr = 2 * fdr / ksteps, offset = offset)
    # active set 
    as <- which(W[ , i] >= t)
    hatS.all[as, i] <- 1
  }
  S <- hatS.all
  
  # Step 2: aggregate the k estimated active sets by taking the union
  hatS <- rep(0, p)
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
