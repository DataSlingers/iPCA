TRCMAimpute_iPCA <- function(x, lamS = NULL, lamDs, q = "multfrob", rmvn.ed1s,
                             seed=1, maxit=50, trace=FALSE, 
                             thr.em=1e-6, thr.glasso=1e-2, maxit.glasso=1) {
  # TRCMAimpute_iPCA: TRCMAimpute modified for iPCA
  # 
  # Inputs:
  #  -x = list of K data matrices with missing values
  #  -lamS = penalty for row covariance matrix (only needed if q == "multfrob")
  #  -lamDs = penalty for column covariance matrix
  #  -q = type of penalty for (column) covariance matrix, either 1, "addfrob", or "multfrob"
  #  -rmvn.ed1s = list of K saved eigenvalue decomposition for Eig(t(xxc) %*% xxc) to speed up computations (optional)
  #  -seed = seed
  #  -maxit = maximum number of iterations
  #  -trace = T/F (default F) whether or not to show loglikelihood stuff
  #  -thr.em = threshold for EM algorithm
  #  -thr.glasso = threshold for GLASSO
  #  -maxit.glasso = max number of iterations for GLASSO
  #
  # Outputs:
  #  -xhat = list of K imputed data matrices
  #  -xh1 = initial imputation assuming Delta = I (not computed)
  #  -xh.init = initial imputation assuming Sig = I
  #  -Sigma = estimated Sigma row covariance
  #  -Delta = estimated Delta column covariance
  #  -Sigma.init = initial estimate of Sigma
  #  -Delta.init = initial estimates for Deltas_ks
  #  -M = estimated mean matrix
  #  -loglike = loglikelihood values after initial imputation and one-step EM algorithm
  #  -rmvn.ed1s = list of K saved eigenvalue decomposition for Eig(t(xxc) %*% xxc)
  
  # record dimensions
  n <- nrow(x[[1]])
  pks <- sapply(x, FUN = ncol)
  p <- sum(pks)
  K <- length(x)
  
  set.seed(seed)
  
  # initial imputation: impute missing values assuming Sig = I
  cat("Initial Imputation... \n")
  rmvn.ans <- xinit <- list()
  Mhatinit <- lapply(X = x, 
                    FUN = function(X) {
                      mu <- colMeans(X, na.rm = T)
                      return(rep(1,n) %*% t(mu))
                    })
  Siginit <- diag(n)
  if (missing(rmvn.ed1s)) {
#     for (k in 1:K) {
#       cat(paste0("... for k = ",k, "\n"))
#       rmvn.ans[[k]] <- RMVNimpute(x=x[[k]], rho=lamDs[k], q=q, trace=F,
#                                   maxit=maxit, thr.em=thr.em, maxit.glasso=maxit.glasso, thr.glasso=thr.glasso)
#     }
    rmvn.ans <- mapply(X = x, D = lamDs,
                         FUN = function(X,D) {return(RMVNimpute(X,D,q=q,trace=F,maxit=maxit,thr.em=thr.em,maxit.glasso=25,thr.glasso=thr.glasso))},
                         SIMPLIFY = FALSE)
    xinit <- lapply(X = rmvn.ans, FUN = function(X) {return(X$xhat)})
    Delthats <- Deltinits <- lapply(X = rmvn.ans, FUN = function(X) {return(X$Sig)})
    rmvn.ed1s <- lapply(X = rmvn.ans, FUN = function(X) {return(X$rmvn.ed1)})
  }else {
#     for (k in 1:K) {
#       cat(paste0("... for k = ",k, "\n"))
#       xinit[[k]] <- RMVNimpute(x=x[[k]], rho=lamDs[k], q=q, ed1=rmvn.ed1s[[k]], trace=F,
#                                maxit=maxit, thr.em=thr.em, maxit.glasso=maxit.glasso, thr.glasso=thr.glasso)
#     }
    xinit <- mapply(X = x, D = lamDs, E = rmvn.ed1s,
                      FUN = function(X,D,E) {return(RMVNimpute(X,D,q=q,ed1=E,trace=F,maxit=maxit,thr.em=thr.em,maxit.glasso=25,thr.glasso=thr.glasso)$xhat)},
                      SIMPLIFY = FALSE)
    Delthats <- Deltinits <- lapply(X = rmvn.ans, FUN = function(X) {return(X$Sig)})
  }
  
  xhat <- xinit
  like <- NA
  like2 <- NA
  like_idx <- 1
  if (trace) {
    like[like_idx] <- MNloglike_iPCA(Xs = x, Ms = Mhatinit, Sig = diag(n), Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like2[like_idx] <- MNloglike_iPCA(Xs = xhat, Ms = Mhatinit, Sig = diag(n), Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like_idx <- like_idx + 1
  }
  
  # compute MLEs given the initial imputed xhat
  cat("Computing MLEs... \n")
  xc <- lapply(xhat, FUN = function(X) {return(scale(X,center = T, scale = F))}) # centered data
  muhat <- lapply(xc, FUN = function(X) {return(attr(X, "scaled:center"))}) # estimated column means
  Mhat <- lapply(muhat, FUN = function(X) {return(matrix(rep(X, times = n), nrow = n, byrow = T))})
  
  if (q == 1) {
    ans <- FFmleGlasso(dat = xc, lamDs = lamDs, lamS = lamS, maxit = maxit.glasso, thr = thr.glasso)
  }else if (q == "1_off") {
    ans <- FFmleGlasso(dat = xc, lamDs = lamDs, lamS = lamS, maxit = maxit.glasso, thr = thr.glasso, pen_diag = F)
  }else if (q == "corr1_off"){
    ans <- FFmleGlassoCorrelation(dat = xc, lamDs = lamDs, lamS = lamS, maxit = maxit.glasso, thr = thr.glasso, pen_diag = F)
  }else if (q == "addfrob") {
    ans <- FFmleAddFrob(dat = xc, lamDs = lamDs, lamS = lamS)
  }else if (q == "multfrob") {
    ans <- FFmleMultFrob(dat = xc, lamDs = lamDs)
  }else {
    stop("q must be either 1, 1_off, addfrob, or multfrob.")
  }
  
  Sighat <- ans$Sig
  Delthats <- ans$Delts

  # given MLEs for covariances, impute X using ECM algorithm from GA's paper
  cat("Running ECM Algorithm... \n")

#   for (k in 1:K) {
#     X <- x[[k]]
#     pk <- pks[k]
#     xhatk <- xhat[[k]]
#     Delthatk <- Delthats[[k]]
#     Mhatk <- Mhat[[k]]
#     
#     ind <- 1; iter <- 1
#     rmi <- (1:n)[apply(is.na(X),1,sum)>0] # indices of rows with missing values
#     cmj <- (1:pk)[apply(is.na(X),2,sum)>0] # indices of columns with missing values
#     na_idx <- is.na(X)
#     
#     while (ind > thr.em & iter < maxit) { # same algorithm as GA but specialized to the case where nu = 0
#       oldx <- xhatk
#       oldxc <- scale(oldx, center = T, scale = F)
#       
#       #by rows
#       for(i in rmi) {
#         mi <- na_idx[i,]
#         sinvs <- t(solve(Sighat[-i,-i], Sighat[-i,i]))
#         psi <- Mhatk[i,] + sinvs %*% (xhatk[-i,] - Mhatk[-i,])
#         Gam <- (Sighat[i,i] - sinvs %*% Sighat[-i,i]) %x% Delthatk
#         ginvg <- t(solve(Gam[!mi,!mi], Gam[!mi,mi]))
#         xhatk[i,mi] <- psi[mi] + ginvg %*% (xhatk[i,!mi] - psi[!mi])
#       }
#       
#       #by cols
#       for(j in cmj) {
#         mj <- na_idx[,j] # T/F: missing or not
#         invdd <- solve(Delthatk[-j,-j], Delthatk[-j,j])
#         nu <- Mhatk[,j] + (xhatk[,-j] - Mhatk[,-j]) %*% invdd
#         Phi <- Sighat %x% (Delthatk[j,j] - Delthatk[j,-j] %*% invdd)
#         pinvp <- t(solve(Phi[!mj,!mj], Phi[!mj,mj]))
#         xhatk[mj,j] <- nu[mj] + pinvp %*% (xhatk[!mj,j] - nu[!mj])
#       }
#       
#       ind <- sum((oldx[na_idx] - xhatk[na_idx])^2)/sum(oldxc[na_idx]^2)
#       iter <- iter + 1
#     }
#     xhat[[k]] <- xhatk
#   }

  xhat <- mapply(X = x, pk = as.list(pks), xhatk = xhat, Delthatk = Delthats, Mhatk = Mhat, SIMPLIFY = F,
                   FUN = function(X, pk, xhatk, Delthatk, Mhatk) {
                     ind <- 1; iter <- 1
                     rmi <- (1:n)[apply(is.na(X),1,sum)>0] # indices of rows with missing values
                     cmj <- (1:pk)[apply(is.na(X),2,sum)>0] # indices of columns with missing values
                     na_idx <- is.na(X)
                     
                     while (ind > thr.em & iter < maxit) { # same algorithm as GA but specialized to the case where nu = 0
                       oldx <- xhatk
                       oldxc <- scale(oldx, center = T, scale = F)
                       
                       #by rows
                       for(i in rmi) {
                         mi <- na_idx[i,]
                         sinvs <- t(solve(Sighat[-i,-i], Sighat[-i,i]))
                         psi <- Mhatk[i,] + sinvs %*% (xhatk[-i,] - Mhatk[-i,])
                         Gam <- (Sighat[i,i] - sinvs %*% Sighat[-i,i]) %x% Delthatk
                         ginvg <- t(solve(Gam[!mi,!mi], Gam[!mi,mi]))
                         xhatk[i,mi] <- psi[mi] + ginvg %*% (xhatk[i,!mi] - psi[!mi])
                       }
                       
                       #by cols
                       for(j in cmj) {
                         mj <- na_idx[,j] # T/F: missing or not
                         invdd <- solve(Delthatk[-j,-j], Delthatk[-j,j])
                         nu <- Mhatk[,j] + (xhatk[,-j] - Mhatk[,-j]) %*% invdd
                         Phi <- Sighat %x% (Delthatk[j,j] - Delthatk[j,-j] %*% invdd)
                         pinvp <- t(solve(Phi[!mj,!mj], Phi[!mj,mj]))
                         xhatk[mj,j] <- nu[mj] + pinvp %*% (xhatk[!mj,j] - nu[!mj])
                       }
                       
                       ind <- sum((oldx[na_idx] - xhatk[na_idx])^2)/sum(oldxc[na_idx]^2)
                       iter <- iter + 1
                     }
                     return(xhatk)
                   })
  
  if (trace) {
    like[like_idx] <- MNloglike_iPCA(Xs = x, Ms = Mhat, Sig = Sighat, Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like2[like_idx] <- MNloglike_iPCA(Xs = xhat, Ms = Mhat, Sig = Sighat, Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like_idx <- like_idx + 1
  }
  
  return(list(xhat = xhat, xh.init = xinit, M.init = Mhatinit,
              Sigma.init = Siginit, Delta.init = Deltinits,
              Sigma = Sighat, Delta = Delthats,M = Mhat,
              loglike = like, loglike2 = like2, rmvn.ed1s = rmvn.ed1s))
}
