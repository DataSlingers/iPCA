# CMF: computes collective matrix factorization estimates
#
# inputs: (must supply either draw or dat)
#   draw = list of 2 elements:
#             1) sim_data: K-list of n x pk CENTERED data matrices
#             2) truth: list of true Sigma and true Deltks
#   dat = K-list of n x pk CENTERED data matrices
#   r = rank of output matrices
#   lamDs = K-length vector of penalty constants for Delti_k's
#   lamS = positive number; penalty constant for Sigi
#   maxit = maximum number of iterations
#   thr = threshold for stopping criterion
#   init = initialization; K+1 cell of symmetric positive definite matrices [[U, V1, ..., VK]]
#
# outputs:
#   U = n x r row covariance matrix estimate
#   Vs = K-list of pk x r column covariance matrix estimates
#
# usage:
# dat <- list()
# dat[[1]] <- matrix(rnorm(50), nrow = 5, ncol = 10)
# dat[[2]] <- matrix(rnorm(75), nrow = 5, ncol = 15)
# r <- 2
# lamDs <- c(2,1)
# lamS <- 1
# ans <- CMF(dat = dat, r = r, lamDs = lamDs, lamS = lamS)

library(parallel)
library(pracma)

CMF_eval <- function(Xs, U, Vs, lamDs, lamS) {
  # evaluate CMF objective function with missing or non-missing data
  K <- length(Xs)
  val <- 0
  for (k in 1:K) {
    val <- val + sum((Xs[[k]] - U %*% t(Vs[[k]]))^2, na.rm = T) + lamDs[k] * norm(Vs[[k]], "F")^2
  }
  val <- val + lamS * norm(U, "F")^2
  return(val)
}


CMF <- function(draw, dat, r, lamDs, lamS,
                maxit = 1000, thr = 1e-6, init) {
  
  if (!(missing(dat))) {
    Xs <- dat
  }else {
    Xs <- draw$sim_data
  }
  
  K <- length(Xs)
  
  if (missing(init)) {
    init <- c(list(diag(nrow(Xs[[1]]))[,1:r]), lapply(X = Xs, FUN = function(X) {return(diag(ncol(X))[,1:r])}))
  }
  if (missing(lamDs)) {
    lamDs <- rep(0, K)
    lamS <- 0
  }
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  
  # initialize
  U <- init[[1]]
  Vs <- init[2:(K+1)]
  
  ind <- 1; iter <- 1; inds <- rep(0, each = maxit)
  new_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs, lamDs = lamDs, lamS = lamS)
  while (ind > thr & iter < maxit) {
    # print(iter)
    # oldU <- U
    old_val <- new_val
    
    # update U
    Uterm1 <- Reduce("+", mapply(X = Xs, V = Vs,
                                 FUN = function(X, V) {return(X %*% V)},
                                 SIMPLIFY = FALSE))
    Uterm2 <- Reduce("+", lapply(X = Vs, FUN = function(X) {return(t(X) %*% X)}))
    if (lamS == 0) { # no regularization
      U <- Uterm1 %*% pinv(Uterm2)
    }else {
      U <- Uterm1 %*% solve(Uterm2 + lamS * diag(r))
    }
    
    if (identical(lamDs, rep(0, K))) { # no regularization
      Vs <- mapply(X = Xs, lamD = as.list(lamDs),
                   FUN = function(X, lamD) {return(t(X) %*% U %*% pinv(t(U) %*% U))},
                   SIMPLIFY = FALSE)
    }else {
      Vs <- mapply(X = Xs, lamD = as.list(lamDs),
                   FUN = function(X, lamD) {return(t(X) %*% U %*% solve(t(U) %*% U + lamD * diag(r)))},
                   SIMPLIFY = FALSE)
    }
    
    new_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs, lamDs = lamDs, lamS = lamS)
    ind <- old_val - new_val
    # ind <- norm(oldU - U, type = "F") / norm(oldU, type = "F")
    inds[iter] <- ind
    iter <- iter + 1
  }
  
  if (ind > thr) {
    message('Warning: CMF estimate did not converge.')
    # plot(log10(inds))
  }else {
    cat(paste0('CMF estimate converged after ', iter-1, ' iterations. \n'))
  }
  
  if (!(missing(draw))) {
    if ("truth" %in% names(draw)) {
      return(list(U = U, Vs = Vs,
                  truth = draw$truth))
    }
  }
  
  return(list(U = U, Vs = Vs))
}

# CMFGradDescent <- function(draw, dat, r, lamDs, lamS,
#                            maxit = 1000, thr = 1e-6, init) {
#   
#   if (!(missing(dat))) {
#     Xs <- dat
#   }else {
#     Xs <- draw$sim_data
#   }
#   
#   K <- length(Xs)
#   
#   if (missing(init)) {
#     init <- c(list(diag(nrow(Xs[[1]]))[,1:r]), lapply(X = Xs, FUN = function(X) {return(diag(ncol(X))[,1:r])}))
#   }
#   if (missing(lamDs)) {
#     lamDs <- rep(0, K)
#     lamS <- 0
#   }
#   
#   # record dimensions
#   n <- nrow(Xs[[1]])
#   pks <- sapply(Xs, FUN = ncol)
#   p <- sum(pks)
#   
#   # initialize
#   U <- init[[1]]
#   Vs <- init[2:(K+1)]
#   
#   ind <- 1; iter <- 1; inds <- rep(0, each = maxit)
#   new_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs, lamDs = lamDs, lamS = lamS)
#   while (ind > thr & iter < maxit) {
#     # print(iter)
#     # oldU <- U
#     old_val <- new_val
#     
#     # compute gradient of U
#     Uterm1 <- Reduce("+", mapply(X = Xs, V = Vs,
#                                  FUN = function(X, V) {
#                                    return((X - U %*% t(V)) %*% V)
#                                  }, SIMPLIFY = FALSE))
#     Ugrad <- -2 * Uterm1 + 2 * lamS * U
#     
#     # perform backtracking line search to choose U's step size
#     Ualpha <- 1
#     search_val <- CMF_eval(Xs = Xs, U = U - Ualpha * Ugrad,
#                            Vs = Vs, lamDs = lamDs, lamS = lamS)
#     orig_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs, lamDs = lamDs, lamS = lamS)
#     while (search_val > orig_val - Ualpha / 2 * sum(Ugrad^2)) {
#       Ualpha <- Ualpha / 2
#       search_val <- CMF_eval(Xs = Xs, U = U - Ualpha * Ugrad,
#                              Vs = Vs, lamDs = lamDs, lamS = lamS)
#     }
#     
#     # update U
#     U <- U - Ualpha * Ugrad
#     
#     # update V
#     for (k in 1:K) {
#       # compute gradient of Vk
#       Vgrad <- -2 * t(Xs[[k]] - U %*% t(Vs[[k]])) %*% U + 2 * lamDs[k] * Vs[[k]]
#       
#       # perform backtracking line search to choose Vk's step size
#       Valpha <- 1
#       Vs_tmp <- Vs
#       Vs_tmp[[k]] <- Vs[[k]] - Valpha * Vgrad
#       search_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs_tmp, lamDs = lamDs, lamS = lamS)
#       orig_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs, lamDs = lamDs, lamS = lamS)
#       while (search_val > orig_val - Valpha / 2 * sum(Vgrad^2)) {
#         Valpha <- Valpha / 2
#         Vs_tmp <- Vs
#         Vs_tmp[[k]] <- Vs[[k]] - Valpha * Vgrad
#         search_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs_tmp, 
#                                lamDs = lamDs, lamS = lamS)
#       }
# 
#       # update Vk
#       Vs[[k]] <- Vs[[k]] - Valpha * Vgrad
#     }
#     
#     new_val <- CMF_eval(Xs = Xs, U = U, Vs = Vs, lamDs = lamDs, lamS = lamS)
#     ind <- old_val - new_val
#     # ind <- norm(oldU - U, type = "F") / norm(oldU, type = "F")
#     inds[iter] <- ind
#     iter <- iter + 1
#   }
#   
#   if (ind > thr) {
#     message('Warning: CMF estimate did not converge.')
#     plot(log10(inds))
#   }else {
#     cat(paste0('CMF estimate converged after ', iter-1, ' iterations. \n'))
#   }
#   
#   if (!(missing(draw))) {
#     if ("truth" %in% names(draw)) {
#       return(list(U = U, Vs = Vs,
#                   truth = draw$truth))
#     }
#   }
#   
#   return(list(U = U, Vs = Vs))
# }
