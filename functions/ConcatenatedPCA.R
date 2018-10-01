# ConcatenatedPCA: computes concatenated PCA 
#     (i.e. naive covariance estimators) for coupled data matrices
#
# inputs:
#   dat = K-list of n x pk CENTERED data matrices
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#   Delts = K-list of pk x pk column covariance matrix estimates
#
library(parallel)

ConcatenatedPCA <- function(dat) {
  
  Xs <- dat
  
  K <- length(Xs)
  
  # check and record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  ns <- sapply(Xs, FUN = nrow)
  if (all.equal(ns, rep(n,K)) == F) {
    stop('Sample size n must be constant among X_k.')
  }
  
  # center data
  for (k in 1:K) {
    Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
  }
  
  # Compute Sig = 1/p * sum_k Xk %*% t(Xk)
  Sig <- 1/p * Reduce("+", lapply(X = Xs, FUN = function(X) {return(X %*% t(X))}))
  Sig <- 1/2 * (Sig + t(Sig)) # to ensure symmetric
    
  # Compute Delts = 1/n * t(Xk) %*% Xk
  Delts <- lapply(X = Xs, FUN = function(X) {return(1/n * t(X) %*% X)})
  Delts <- lapply(X = Delts, FUN = function(X) {return(1/2 * (X + t(X)))}) # to ensure symmetric
  
  return(list(Sig = Sig, Delts = Delts))

}
                  
