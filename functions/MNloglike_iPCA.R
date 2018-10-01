MNloglike_iPCA <- function(Xs,Ms,Sig,Delts,lamS,lamDs,q,Sigi=NULL,Deltis=NULL,print=F)
{
  # MNloglike_iPCA: compute matrix-variate log-likelihood for iPCA model
  #
  # Inputs:
  #  -Xs: K-list of n x pk data matrices
  #  -Ms: K-list of mean matrices
  #  -Sig: n x n row covariance matrix
  #  -Delts: K-list of pk x pk column covariance matrix
  #  -lamS: a scalar; penalty constant for Sigi
  #  -lamDs: K-length vector of penalty constants for Delti_k's
  #  -q: type of penalty; either 1, "addfrob" or "multfrob"
  #  -Sigi: inverse of Sig; optional
  #  -Deltis: K-length listo f Deltk inverses; optional
  #  -print: logical T/F; whether or not to print log-likelihood value
  #
  # Output:
  #  -val: value of the iPCA log-likelihood
  
  K <- length(Xs)
  n_cores <- K
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  
  # fill in missing function parameters
  if (length(Sigi)==0){ Sigi <- solve(Sig) }
  if (length(Deltis)==0){ Deltis <- mclapply(X = Delts, FUN = function(X) {chol2inv(chol(X))}, mc.cores = n_cores) } # invert Delts 
  
  # center matrix
  if (missing(Ms)) {
    Xc <- lapply(X = Xs, FUN = function(X) {return(scale(X, center = T, scale = F))})
  }else {
    Xc <- mapply(X = Xs, M = Ms, SIMPLIFY = F,
                 FUN = function(X,M) {return(X-M)})
  }
  
  # compute penalty terms
  if (q == 1) {
    pen <- lamS*sum(abs(Sigi)) + sum(lamDs * mcmapply(X = Deltis, 
                                                      FUN = function(X) {return(sum(abs(X)))}, 
                                                      SIMPLIFY = TRUE, mc.cores = n_cores))
  }else if (q == "addfrob") {
    pen <- lamS*sum(Sigi^2) + sum(lamDs * mcmapply(X = Deltis, 
                                                   FUN = function(X) {return(sum(X^2))}, 
                                                   SIMPLIFY = TRUE, mc.cores = n_cores))
  }else if (q == "multfrob") {
    pen <- sum(Sigi^2) * sum(lamDs * mcmapply(X = Deltis, 
                                              FUN = function(X) {return(sum(X^2))}, 
                                              SIMPLIFY = TRUE, mc.cores = n_cores))
  }else {
    stop("q must be either 1, addfrob, or multfrob.")
  }
  
  if (sum(mapply(Xs, FUN = function(X) {return(sum(is.na(X)))})) == 0) { # if no missing data
    t1 <- determinant(Sigi, logarithm = T)$modulus
    t2 <- mcmapply(X = Deltis,
                   FUN = function(X) {return(determinant(X, logarithm = T)$modulus)},
                   SIMPLIFY = TRUE, mc.cores = n_cores)
    t3 <- mcmapply(X = Xc, D = Deltis, 
                   FUN = function(X,D) {return(sum(diag(Sigi %*% X %*% D %*% t(X))))},
                   SIMPLIFY = TRUE, mc.cores = n_cores)
    val <- p*t1 + n*sum(t2) - sum(t3) - pen
    
  }else {
    tlds <- 0
    tldd <- trt <- 0
    for (k in 1:K) {
      Xk <- x <- Xs[[k]]
      Deltk <- Delts[[k]]
      for (j in 1:pks[k]) {
        oj <- !is.na(x[,j])
        sigioj <- solve(Sig[oj,oj])
        tlds <- tlds + determinant(sigioj, logarithm = T)$modulus
        Xk[oj, j] <- chol(sigioj) %*% Xc[[k]][oj, j]
      } 
      for (i in 1:n) {
        oi <- !is.na(x[i,])
        deltioi <- solve(Deltk[oi, oi])
        tldd <- tldd + determinant(deltioi, logarithm = T)$modulus
        trt <- trt + sum(diag(Xk[i, oi] %*% t(Xk[i, oi]) %*% deltioi))
      }
    }
    val <- tlds + tldd - trt - pen
  }
  
  if (print) {
    cat(val, fill = T)
  }
  
  return(val)
}
  
  