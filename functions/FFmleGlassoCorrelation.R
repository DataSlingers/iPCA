# FFmleGlassoCorrelation: computes iPCA covariance estimates via Flip-Flop MLE
# with additive Glasso penalty (using correlations)
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

FFmleGlassoCorrelation <- function(draw, dat, lamDs, lamS, 
                                   maxit = 1, thr = 1e-3, init,
                                   maxit.glasso = 250, thr.glasso = 1e-4,
                                   pen_diag = F) {
  
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
    Shat_Sig <- Reduce("+", mapply(X = Xs, D = Deltis,
                                   FUN = function(X, D) {return(X %*% D %*% t(X))},
                                   SIMPLIFY = FALSE))
    Shat_Sig <- 1/(2*p) * (Shat_Sig + t(Shat_Sig)) # to ensure symmetry
    SD_Sig <- diag(sqrt(diag(Shat_Sig)))
    Sphat_Sig <- solve(SD_Sig) %*% Shat_Sig %*% solve(SD_Sig)
    
    # cat(paste0('Glasso Iteration ', iter, ': Sigma \n'))
    if (pen_diag == T) {
      quic_res <- QUIC(S = Shat_Sig, rho = lamS/p, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }else {
      rho_off <- lamS/p * (matrix(1, nrow = n, ncol = n) - diag(n))
      quic_res <- QUIC(S = Shat_Sig, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }
    Sigi <- solve(SD_Sig) %*% quic_res$X %*% solve(SD_Sig)
    Sig <- SD_Sig %*% quic_res$W %*% SD_Sig
    
    # to ensure symmetry
    Sigi <- 1/2 * (Sigi + t(Sigi))
    Sig <- 1/2 * (Sig + t(Sig))
  }
  
  iter <- iter + 1
  while (ind > thr & iter <= maxit) {
    oldS <- Sigi
    
    # update Deltks
    Shat_ks <- lapply(X = Xs,
                      FUN = function(X) {
                        Shat_k <- 1/n * t(X) %*% Sigi %*% X
                        Shat_k <- 1/2 * (Shat_k + t(Shat_k)) # to ensure symmetry
                        return(Shat_k)
                      })
    SD_ks <- lapply(X = Shat_ks,
                    FUN = function(X) {
                      return(diag(sqrt(diag(X))))
                    })
    Sphat_ks <- mapply(X = Shat_ks, W = SD_ks,
                       FUN = function(X, W) {
                         return(solve(W) %*% X %*% solve(W))
                       },
                       SIMPLIFY = F)
    
    for (k in 1:K) {
      # cat(paste0('Glasso Iteration ', iter, ': Delta', k, '\n'))
      if (pen_diag == T) {
        quic_res <- QUIC(S = Sphat_ks[[k]], rho = lamDs[k]/n, tol = thr.glasso, msg = 0, maxIter = maxit.glasso,
                         X.init = Deltis[[k]], W.init = Delts[[k]])
      }else {
        pk <- pks[k]
        rho_off <- lamDs[k]/n * (matrix(1, nrow = pk, ncol = pk) - diag(pk))
        quic_res <- QUIC(S = Sphat_ks[[k]], rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso,
                         X.init = Deltis[[k]], W.init = Delts[[k]])
      }
      Deltis[[k]] <- solve(SD_ks[[k]]) %*% quic_res$X %*% solve(SD_ks[[k]])
      Delts[[k]] <- SD_ks[[k]] %*% quic_res$W %*% SD_ks[[k]]
    }
    
    Deltis <- lapply(X = Deltis, FUN = function(X) {return(1/2 * (X + t(X)))}) # to ensure symmetry
    Delts <- lapply(X = Delts, FUN = function(X) {return(1/2 * (X + t(X)))}) # to ensure symmetry
    
    # update Sig
    Shat_Sig <- Reduce("+", mapply(X = Xs, D = Deltis,
                                   FUN = function(X, D) {return(X %*% D %*% t(X))},
                                   SIMPLIFY = FALSE))
    Shat_Sig <- 1/(2*p) * (Shat_Sig + t(Shat_Sig)) # to ensure symmetry
    SD_Sig <- diag(sqrt(diag(Shat_Sig)))
    Sphat_Sig <- solve(SD_Sig) %*% Shat_Sig %*% solve(SD_Sig)
    
    # cat(paste0('Glasso Iteration ', iter, ': Sigma \n'))
    if (pen_diag == T) {
      quic_res <- QUIC(S = Shat_Sig, rho = lamS/p, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }else {
      rho_off <- lamS/p * (matrix(1, nrow = n, ncol = n) - diag(n))
      quic_res <- QUIC(S = Shat_Sig, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigi, W.init = Sig)
    }
    Sigi <- solve(SD_Sig) %*% quic_res$X %*% solve(SD_Sig)
    Sig <- SD_Sig %*% quic_res$W %*% SD_Sig
    
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
