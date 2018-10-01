RMNimpute_iPCA <- function(x, lamS, lamDs, q, cov.init, xhat.init, Mhat.init,
                           seed=1, maxit=50, trace=FALSE,
                           thr.em=1e-4,thr.glasso=1e-4,maxit.glasso=10) {
  
  # RMNimpute_iPCA: RMNimpute (from GA) modified to impute matrix-variate normal model for iPCA
  #
  # Inputs:
  #  -x = list of K data matrices with missing values
  #  -lamS = penalty for row covariance matrix (only needed if q == "multfrob")
  #  -lamDs = penalty for column covariance matrix
  #  -q = type of penalty for (column) covariance matrix, either 1, "addfrob", or "multfrob"
  #  -cov.init = initialization for covariances; K+1 cell of symmetric positive definite matrices [[Sig, Delt1, ..., DeltK]]; optional but default is the identity
  #  -xhat.init = initialization for xhat 
  #  -Mhat.init = initialization for Mhat 
  #  -seed = seed
  #  -maxit = maximum number of iterations
  #  -trace = T/F (default F) whether or not to show loglikelihood stuff
  #  -thr.em = threshold for EM algorithm
  #  -thr.glasso = threshold for GLASSO
  #  -maxit.glasso = max number of iterations for GLASSO
  #
  # Outputs:
  #  -xhat = list of K imputed data matrices
  #  -Sigma = estimated Sigma row covariance
  #  -Delta = list of K estimated Delta_k column covariance
  #  -M = list of K estimated mean matrices
  #  -loglike = loglikelihood values after each iteration
  
  # record dimensions
  n <- nrow(x[[1]])
  pks <- sapply(x, FUN = ncol)
  p <- sum(pks)
  K <- length(x)
  n_cores <- K
  set.seed(seed)
  
  # estimate column means
  if (!missing(Mhat.init) | !missing(xhat.init)) {
    Mhat <- Mhat.init
    xhat <- xhat.init
  }else {
    xcna <- lapply(X = x, FUN = function(X) {return(scale(X, center = T, scale = F))})
    muhat <- lapply(X = xcna, FUN = function(X) {return(attr(X, "scaled:center"))}) # estimated column means
    Mhat <- lapply(muhat, FUN = function(X) {return(matrix(rep(X, times = n), nrow = n, byrow = T))})
    xc <- lapply(X = xcna, 
                 FUN = function(X) {
                   X[is.na(X)] <- 0
                   return(X) })
    xhat <- Map("+", xc, Mhat)
  }
  
  if (missing(cov.init)) {
    cov.init <- c(list(diag(n)), lapply(X = x, FUN = function(X) {return(diag(ncol(X)))}))
  }
  
  # initialize
  Sighat <- cov.init[[1]]
  Sigihat <- chol2inv(chol(Sighat))
  Delthats <- cov.init[2:(K+1)]
  Deltihats <- mclapply(X = Delthats, FUN = function(X) {chol2inv(chol(X))}, mc.cores = n_cores) # invert Deltis
  
  like <- NA
  if (trace) {
    cat("Computing loglike... ")
    like[1] <- MNloglike_iPCA(Xs = x, Ms = Mhat, Sig = Sighat, Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like_idx <- 2
  }
  
  ind <- 1; iter <- 1
  cat("...Starting while loop... \n")
  while(ind>thr.em & iter<=maxit) {
    oldxh <- xhat
    print(iter)
    iter <- iter + 1
    
    #by cols
    ##E step (Sig)
    cat("E step for Sig... \n")
    covcorrs <- list()
    for (k in 1:K) { # compute E(...Xk...) for each k
      # define useful variables
      pk <- pks[k]
      xk <- x[[k]]
      xhatk <- xhat[[k]]
      Delthatk <- Delthats[[k]]
      Mhatk <- Mhat[[k]]
      covcorr <- matrix(0,n,n)
      na_idx <- is.na(xk) # logical matrix of missing/observed values
      cmj <- (1:pk)[apply(is.na(xk),2,sum)>0] # indices of columns with missing values
      for(j in cmj) {
        mj <- na_idx[,j] # T/F: missing or not
        invdd <- solve(Delthatk[-j,-j], Delthatk[-j,j])
        nu <- Mhatk[,j] + (xhatk[,-j] - Mhatk[,j]) %*% invdd
        Phi <- Sighat %x% (Delthatk[j,j] - Delthatk[j,-j] %*% invdd)
        pinvp <- t(solve(Phi[!mj,!mj], Phi[!mj,mj]))
        xhatk[mj,j] <- nu[mj] + pinvp %*% (xhatk[!mj,j] - nu[!mj])
        covcorr[mj,mj] <- covcorr[mj,mj] + Phi[mj,mj] - pinvp %*% Phi[!mj,mj]
      }
      xhat[[k]] <- xhatk
      covcorrs[[k]] <- covcorr
    }

    ##M step (update Sig)
    cat("M step for Sig... \n")
    
    # first, update means
    xhc <- lapply(X = xhat, FUN = function(X) {return(scale(X, center = T, scale = F))})
    muhat <- lapply(X = xhc, FUN = function(X) {return(attr(X, "scaled:center"))}) # estimated column means
    Mhat <- lapply(muhat, FUN = function(X) {return(matrix(rep(X, times = n), nrow = n, byrow = T))})
    
    inS <- Reduce("+", mcmapply(X = xhc, D = Deltihats, C = covcorrs,
                                FUN = function(X,D,C) {return(X %*% D %*% t(X) + C)},
                                SIMPLIFY = FALSE, mc.cores = n_cores))
    
    # update Sig
    if (q == 1) {
      quic_res <- QUIC(S = inS/p, rho = lamS/(2*p), tol = 1e-2, msg = 0, maxIter = 100, 
                       X.init = Sigihat, W.init = Sighat)
      Sigihat <- quic_res$X; Sigihat <- 1/2 * (Sigihat + t(Sigihat))
      Sighat <- quic_res$W; Sighat <- 1/2 * (Sighat + t(Sighat))
    }else if (q == "addfrob") {
      Sig_eigs <- eigen(inS)
      Sig_V <- Sig_eigs$vectors; Sig_d <- Sig_eigs$values
      gams <- 1/(2*p) * (Sig_d + sqrt(Sig_d^2 + 8*p*lamS))
      Sigihat <- Sig_V %*% diag(1/gams) %*% t(Sig_V)
      Sighat <- Sig_V %*% diag(gams) %*% t(Sig_V)
    }else if (q == "multfrob") {
      Sig_eigs <- eigen(inS)
      Sig_V <- Sig_eigs$vectors; Sig_d <- Sig_eigs$values
      sumDs <- sum(lamDs * sapply(X = Deltihats, FUN = function(X) {return(norm(X, "F")^2)}))
      gams <- 1/(2*p) * (Sig_d + sqrt(Sig_d^2 + 8*p*sumDs))
      Sigihat <- Sig_V %*% diag(1/gams) %*% t(Sig_V)
      Sighat <- Sig_V %*% diag(gams) %*% t(Sig_V)
    }else {
      stop("q must be either 1, addfrob, or multfrob.")
    }
    
    covcorcs <- list()
    for (k in 1:K) {
      # define useful variables
      pk <- pks[k]
      xk <- x[[k]]
      xhatk <- xhat[[k]]
      Delthatk <- Delthats[[k]]
      Mhatk <- Mhat[[k]]
      na_idx <- is.na(xk) # logical matrix of missing/observed values
      
      #by rows
      ##Estep (Deltk)
      cat(paste0("E step for Delta_",k, "... \n"))
      covcorc <- matrix(0,pk,pk);
      rmi <- (1:n)[apply(is.na(xk),1,sum)>0] # indices of rows with missing values
      for(i in rmi) {
        mi <- na_idx[i,]
        sinvs <- t(solve(Sighat[-i,-i], Sighat[-i,i]))
        psi <- Mhatk[i,] + sinvs %*% (xhatk[-i,] - Mhatk[-i,])
        Gam <- (Sighat[i,i] - sinvs %*% Sighat[-i,i]) %x% Delthatk
        ginvg <- t(solve(Gam[!mi,!mi], Gam[!mi,mi]))
        xhatk[i,mi] <- psi[mi] + ginvg %*% (xhatk[i,!mi] - psi[!mi])
        covcorc[mi,mi] <- covcorc[mi,mi] + Gam[mi,mi] - ginvg %*% Gam[!mi,mi]
      }
      xhat[[k]] <- xhatk
      covcorcs[[k]] <- covcorc
      
      ##M step (Deltk)...
      # let's continue the M step for Deltk outside of the foor loop and do everything in parallel
    }
    
    cat("M step for All Delta's... \n")
    # put means in lists
    xhc <- lapply(X = xhat, FUN = function(X) {return(scale(X, center = T, scale = F))})
    muhat <- lapply(X = xhc, FUN = function(X) {return(attr(X, "scaled:center"))}) # estimated column means
    Mhat <- lapply(muhat, FUN = function(X) {return(matrix(rep(X, times = n), nrow = n, byrow = T))})
    
    inD <- mcmapply(X = xhc, C = covcorcs, 
                    FUN = function(X,C) {return(t(X) %*% Sigihat %*% X + C)},
                    SIMPLIFY = FALSE, mc.cores = n_cores)
    
    # update all Deltks in parallel
    if (q == 1) {
      quic_res <- mcmapply(D = inD, Rho = as.list(lamDs), XInit = Deltihats, WInit = Delthats,
                           FUN = function(D, Rho, XInit, WInit) {
                             return(QUIC(S = D/n, rho = Rho/(2*n), tol = 1e-2, msg = 0,
                                         X.init = XInit, W.init = WInit))},
                           SIMPLIFY = FALSE, mc.cores = n_cores)
      Deltihats <- mclapply(X = quic_res, FUN = function(X) {return(X$X)}, mc.cores = n_cores)
      Delthats <- mclapply(X = quic_res, FUN = function(X) {return(X$W)}, mc.cores = n_cores)
    }else if (q == "addfrob") {
      Delt_eigs <- mclapply(X = inD, FUN = eigen, mc.cores = n_cores)
      Delt_Vs <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$vectors)}, mc.cores = n_cores)
      Delt_ds <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$values)}, mc.cores = n_cores)
      for (k in 1:K) {
        d <- Delt_ds[[k]]
        gams <- (1/(2*n)) * (d + sqrt(d^2 + 8*n*lamDs[k]))
        Deltihats[[k]] <- Delt_Vs[[k]] %*% diag(1/gams) %*% t(Delt_Vs[[k]])
        Delthats[[k]] <- Delt_Vs[[k]] %*% diag(gams) %*% t(Delt_Vs[[k]])
      }
    }else if (q == "multfrob") {
      Delt_eigs <- mclapply(X = inD, FUN = eigen, mc.cores = n_cores)
      Delt_Vs <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$vectors)}, mc.cores = n_cores)
      Delt_ds <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$values)}, mc.cores = n_cores)
      for (k in 1:K) {
        d <- Delt_ds[[k]]
        gams <- (1/(2*n)) * (d + sqrt(d^2 + 8*n*lamDs[k]*norm(Sigihat, "F")^2))
        Deltihats[[k]] <- Delt_Vs[[k]] %*% diag(1/gams) %*% t(Delt_Vs[[k]])
        Delthats[[k]] <- Delt_Vs[[k]] %*% diag(gams) %*% t(Delt_Vs[[k]])
      }
    }else {
      stop("q must be either 1, addfrob, or multfrob.")
    }
    
    if (trace) {
      cat("Computing loglike... ")
      like[like_idx] <- MNloglike_iPCA(Xs = x, Ms = Mhat, Sig = Sighat, Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
      like_idx <- like_idx + 1
    }
    
    ind <- mean(mapply(X = x, Xh = xhat, Xo = oldxh, SIMPLIFY = T,
                       FUN = function(X, Xh, Xo) {
                         miss_idx <- is.na(X)
                         val <- sum((Xh[miss_idx] - Xo[miss_idx])^2)/sum(Xo[miss_idx]^2)
                         }))
  
  }
  
  return(list(xhat=xhat,Sigma=Sighat,Delta=Delthats,M=Mhat,loglike=like))
}
