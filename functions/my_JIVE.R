# my_JIVE: computes JIVE decomposition
#
# inputs:
#   dat = K-list of n x pk CENTERED data matrices
#   rankJ, rankA = see jive() documentation
#   method = (default = "perm") see jive() documentation
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#   Delts = K-list of pk x pk column covariance matrix estimates
#   JIVE_results = results from JIVE function jive()

library(r.jive)

my_JIVE <- function(dat, rankJ, rankA, method = "perm", maxiter = 1000, nperms = 100) {
  
  Xs <- dat
  
  Xs_t <- lapply(Xs, FUN = t)
  
  if (!(missing(rankJ))) {
    results <- my_r_jive(data = Xs_t, rankJ = rankJ, rankA = rankA,
                       method = method, maxiter = maxiter, nperms = nperms)
  }else {
    results <- my_r_jive(Xs_t, method = method, maxiter = maxiter, nperms = nperms)
  }
  # summary(results)
  
  J <- t(do.call(rbind, results$joint))
  Sig <- J %*% t(J)
  Delts <- lapply(results$individual, FUN = function(A) {return(A %*% t(A))})

  return(list(Sig = Sig, Delts = Delts, JIVE_results = results))
}
