# IndividualPCA: computes ordinary PCA on a single data matrix
#
# inputs:
#   dat = K-list of n x pk CENTERED data matrices
#   k = index of data matrix to perform PCA on
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#

IndividualPCA <- function(dat, k) {
  
  Xk <- dat[[k]]
  Xk <- scale(Xk, center = T, scale = F)
  
  Sig <- 1/ncol(Xk) * Xk %*% t(Xk)
  Delts <- 1/nrow(Xk) * t(Xk) %*% Xk
  
  return(list(Sig = Sig, Delts = Delts))
}
                  
