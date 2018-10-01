# FFmleMultFrob: computes iPCA covariance estimates via Flip-Flop MLE
# with Multiplicative Frobenius norm penalty
#
# inputs:
#   dat = K-list of n x pk CENTERED data matrices
#   lamDs = K-length vector of penalty constants for Delti_k's
#   maxit = maximum number of iterations
#   thr = threshold for stopping criterion
#   init = initialization; K+1 cell of symmetric positive definite matrices [[Sigi, Delti1, ..., DeltiK]]
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#   Delts = K-list of pk x pk column covariance matrix estimates
#

library(parallel)

FFmleMultFrob <- function(dat, lamDs,
                         maxit = 1e2, thr = 1e-6, init) {
  
  Xs <- dat
  
  if (missing(init)) {
    init <- c(list(diag(nrow(Xs[[1]]))), lapply(X = Xs, FUN = function(X) {return(diag(ncol(X)))}))
  }
  
  K <- length(Xs)
  n_cores <- 1 # number of cores to use
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  
  # center data
  for (k in 1:K) {
    Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
  }
  
  # initialize
  Sigi <- init[[1]]
  Deltis <- init[2:(K+1)]
  
  ind <- 1; iter <- 1; inds <- rep(0, each = maxit)
  while (ind > thr & iter < maxit) {
    oldS <- Sigi
    
    # update Sig
    inS <- Reduce("+", mcmapply(X = Xs, D = Deltis,
                                FUN = function(X, D) {return(X %*% D %*% t(X))},
                                SIMPLIFY = FALSE,
                                mc.cores = n_cores))
    inS <- 1/2 * (inS + t(inS)) # to ensure symmetry
    
    Sig_eigs <- eigen(inS)
    Sig_V <- Sig_eigs$vectors; Sig_d <- Sig_eigs$values
    
    sumDs <- sum(lamDs * sapply(X = Deltis, FUN = function(X) {return(norm(X, "F")^2)}))
    gams <- 1/(2*p) * (Sig_d + sqrt(Sig_d^2 + 8*p*sumDs))
    
    Sigi <- Sig_V %*% diag(1/gams) %*% t(Sig_V)
    Sig <- Sig_V %*% diag(gams) %*% t(Sig_V)
    
    # update Deltks
    inD <- mclapply(X = Xs, 
                    FUN = function(X) {return(t(X) %*% Sigi %*% X)},
                    mc.cores = n_cores)
    inD <- mclapply(X = inD,
                    FUN = function(X) {return(1/2 * (X + t(X)))},
                    mc.cores = n_cores) # to ensure symmetry
    
    Delt_eigs <- mclapply(X = inD, FUN = eigen, mc.cores = n_cores)
    Delt_Vs <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$vectors)}, mc.cores = n_cores)
    Delt_ds <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$values)}, mc.cores = n_cores)
    
    for (k in 1:K) {
      d <- Delt_ds[[k]]
      gams <- (1/(2*n)) * (d + sqrt(d^2 + 8*n*lamDs[k]*norm(Sigi, "F")^2))
      Deltis[[k]] <- Delt_Vs[[k]] %*% diag(1/gams) %*% t(Delt_Vs[[k]])
    }
    Deltis <- mclapply(X = Deltis, FUN = function(X) {return(1/2 * (X + t(X)))}, mc.cores = n_cores) # to ensure symmetry
    
    iter <- iter + 1
    ind <- sqrt(mean(lamDs))*norm(oldS - Sigi, type = "F") / norm(oldS, type = "F") # to account for differences in penalties
    inds[iter] <- ind
  }
  
  if (ind > thr) {
    message('Warning: Multiplicative Frobenius Sigma estimate did not converge.')
  }else {
    cat(paste0('Multiplicative Frobenius Estimate converged after ', iter-1, ' iterations. \n'))
  }
  
  # compute Delts
  Delts <- mclapply(X = Deltis, FUN = function(X) {chol2inv(chol(X))}, mc.cores = n_cores) # invert Deltis
  
  # to ensure symmetry
  Sig <- 1/2 * (Sig + t(Sig))
  Delts <- mclapply(X = Delts, FUN = function(X) {return(1/2 * (X + t(X)))}, mc.cores = n_cores)
  
  return(list(Sig = Sig, Delts = Delts, lamDs = lamDs))
}
