# DistributedPCA: computes two-step covariance estimators for
# coupled data matrices
#
# inputs:
#   dat = K-list of n x pk CENTERED data matrices
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#   Delts = K-list of pk x pk column covariance matrix estimates
#
library(parallel)

DistributedPCA <- function(dat) {

  Xs <- dat
  
  K <- length(Xs)
  n_cores <- 1 # number of cores
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  
  # center data
  for (k in 1:K) {
    Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
  }
  
  # Step 1: Compute individual Sigmas
  Sig_indiv <- mclapply(X = Xs, FUN = function(X) {return((X %*% t(X)))}, mc.cores = n_cores)
  for (k in 1:K) {
    Sig_indiv[[k]] <- 1/pks[k] * (Sig_indiv[[k]])
  }
  V_indiv <- mclapply(Sig_indiv, FUN = function(Sig) {return(eigen(Sig)$vectors)}, mc.cores = n_cores)
  
  # Step 2: Combine eigenvectors
  Sig <- 1/K * Reduce("+", mclapply(X = V_indiv, FUN = function(X) {return(X %*% t(X))}, mc.cores = n_cores))
  Sig <- 1/2 * (Sig + t(Sig)) # to ensure symmetric
  
  return(list(Sig = Sig))
}
