# FFmleGlasso: computes iPCA covariance estimates via Flip-Flop MLE
# with additive Glasso penalty
#
# inputs: (must supply either draw or dat)
#   draw = list of 2 elements:
#             1) sim_data: K-list of n x pk CENTERED data matrices
#             2) truth: list of true Sigma and true Deltks
#   dat = K-list of n x pk CENTERED data matrices
#   lamDs = K-length vector of penalty constants for Delti_k's
#   lamS = positive number; penalty constant for Sigi
#   maxit = maximum number of iterations
#   thr = threshold for stopping criterion
#   init = initialization; K+1 cell of symmetric positive definite matrices [[Sigi, Delti1, ..., DeltiK]]
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#   Delts = K-list of pk x pk column covariance matrix estimates
#

library(QUIC)

FFmleGlasso <- function(draw, dat, lamDs, lamS, 
                        maxit = 10, thr = 1e-6, init,
                        maxit.glasso = 50, thr.glasso = 1e-2,
                        pen_diag = T) {
  
  if (!(missing(dat))) {
    Xs <- dat
  }else {
    Xs <- draw$sim_data
  }
  
  if (missing(init)) {
    init <- c(list(diag(nrow(Xs[[1]]))), lapply(X = Xs, FUN = function(X) {return(diag(ncol(X)))}))
  }
  
  K <- length(Xs)
  
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
  Sig <- solve(Sigi)
  Deltis <- init[2:(K+1)]
  Delts <- lapply(X = Deltis, FUN = solve)
  
  ind <- 1; iter <- 0; inds <- rep(0, each = maxit)
  
  if (K/sqrt(n) >= sum(sqrt(pks))/p) {
    # if in the large p scenario, estimate sig first
    inS <- Reduce("+", mapply(X = Xs, D = Deltis,
                                FUN = function(X, D) {return(X %*% D %*% t(X))},
                                SIMPLIFY = FALSE))
    inS <- 1/(2*p) * (inS + t(inS)) # to ensure symmetry
    
    # cat(paste0('Glasso Iteration ', iter, ': Sigma \n'))
    if (pen_diag == T) {
      quic_res <- QUIC(S = inS, rho = lamS/p, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }else {
      rho_off <- lamS/p * (matrix(1, nrow = n, ncol = n) - diag(n))
      quic_res <- QUIC(S = inS, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }
    Sigi <- quic_res$X
    Sig <- quic_res$W
    
    # to ensure symmetry
    Sigi <- 1/2 * (Sigi + t(Sigi))
    Sig <- 1/2 * (Sig + t(Sig))
  }
  
  iter <- iter + 1
  while (ind > thr & iter <= maxit) {
    oldS <- Sigi
    
    # update Deltks
    inD <- lapply(X = Xs,
                    FUN = function(X) {return(t(X) %*% Sigi %*% X)})
    inD <- lapply(X = inD,
                    FUN = function(X) {return(1/(2*n) * (X + t(X)))}) # to ensure symmetry
    
      for (k in 1:K) {
        # cat(paste0('Glasso Iteration ', iter, ': Delta', k, '\n'))
        if (pen_diag == T) {
          quic_res <- QUIC(S = inD[[k]], rho = lamDs[k]/n, tol = thr.glasso, msg = 0, maxIter = maxit.glasso,
                           X.init = Deltis[[k]], W.init = Delts[[k]])
        }else {
          pk <- nrow(inD[[k]])
          rho_off <- lamDs[k]/n * (matrix(1, nrow = pk, ncol = pk) - diag(pk))
          quic_res <- QUIC(S = inD[[k]], rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso,
                           X.init = Deltis[[k]], W.init = Delts[[k]])
        }
        Deltis[[k]] <- quic_res$X
        Delts[[k]] <- quic_res$W
      }

    Deltis <- lapply(X = Deltis, FUN = function(X) {return(1/2 * (X + t(X)))}) # to ensure symmetry
    Delts <- lapply(X = Delts, FUN = function(X) {return(1/2 * (X + t(X)))}) # to ensure symmetry
    
    # update Sig
    inS <- Reduce("+", mapply(X = Xs, D = Deltis,
                                FUN = function(X, D) {return(X %*% D %*% t(X))},
                                SIMPLIFY = FALSE))
    inS <- 1/(2*p) * (inS + t(inS)) # to ensure symmetry
    
    # cat(paste0('Glasso Iteration ', iter, ': Sigma \n'))
    if (pen_diag == T) {
      quic_res <- QUIC(S = inS, rho = lamS/p, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }else {
      rho_off <- lamS/p * (matrix(1, nrow = n, ncol = n) - diag(n))
      quic_res <- QUIC(S = inS, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }
    Sigi <- quic_res$X
    Sig <- quic_res$W
    
    # to ensure symmetry
    Sigi <- 1/2 * (Sigi + t(Sigi))
    Sig <- 1/2 * (Sig + t(Sig))
    
    iter <- iter + 1
    ind <- norm(oldS - Sigi, type = "F") / norm(oldS, type = "F")
    # print(ind)
    inds[iter] <- ind
  }
  
  if (ind > thr) {
    message('Warning: Glasso Sigma estimate did not converge.')
  }else {
    cat(paste0('Glasso Estimate converged after ', iter-1, ' iterations. \n'))
  }
  
  # to ensure symmetry
  Sig <- 1/2 * (Sig + t(Sig))
  Delts <- lapply(X = Delts, FUN = function(X) {return(1/2 * (X + t(X)))})
  
  if (!(missing(draw))) {
    if ("truth" %in% names(draw)) {
      return(list(Sig = Sig, Delts = Delts,
                  truth = draw$truth))
    }
  }
  
  return(list(Sig = Sig, Delts = Delts))
  
}
