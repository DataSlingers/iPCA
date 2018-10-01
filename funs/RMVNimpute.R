RMVNimpute <- function(x,rho,cov.corr=TRUE,q,ed1,
                       maxit=10,thr.em=1e-4,trace=FALSE,maxit.glasso=50,thr.glasso=1e-2) {
  # RMVNimpute: function imputes missing values for multivariate normal case
  # 
  # Inputs:
  #  -x = n x p data matrix which arises from multivariate normal
  #  -rho = penalty parameter
  #  -cov.corr = T/F logical; whether or not to compute corrected covariance matrix
  #  -q = penalty type; either 1, "addfrob", "multfrob"
  #  -ed1 = saved eigenvalue decomposition for Eig(t(xc) %*% xc) to speed up computations (optional)
  #  -maxit = max number of iterations
  #  -thr.em = threshold for EM algorithm
  #  -trace = T/F; whether or not to show loglikilihood stuff
  #  -maxit.glasso = max number of iterations for GLASSO
  #  -thr.glasso = threshold for GLASSO
  #
  # Outputs:
  #  -xhat = imputed data matrix
  #  -Sig = estimated Sigma (p x p matrix)
  #  -mu = estimated mean
  #  -loglike = loglikelihood value
  #  -rmvn.ed1 = saved eigenvalue decomposition for Eig(t(xc) %*% xc)
  
  n <- nrow(x)
  p <- ncol(x)
  rmi <- (1:n)[apply(is.na(x),1,sum)>0]  # indices of rows with missing values
  
  if (n > 1) {muhat <- colMeans(x,na.rm=TRUE)} else{muhat <- rep(0,p)}
  
  xc <- t(t(x) - muhat) # center by column means
  xc[is.na(x)] <- 0 # initialize missing data to 0 (equivalent to imputing column means)
  
  # Initialize Sighat by solving regularization optimization problem
  iter <- 1
  Sighat <- t(xc)%*%xc
  if (q==1){ 
    cat(paste0('Glasso Iteration: ', iter, '\n'))
    quic_res <- QUIC(S = Sighat/n, rho = rho/n, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                     X.init = diag(p), W.init = diag(p))
    Sigihat <- quic_res$X
    Sighat <- quic_res$W
    ed1 <- NULL
  }else if (q == "1_off"){
    cat(paste0('Glasso Iteration: ', iter, '\n'))
    tmp_p <- nrow(Sighat)
    rho_off <- rho/n * (matrix(1, nrow = tmp_p, ncol = tmp_p) - diag(tmp_p))
    quic_res <- QUIC(S = Sighat/n, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                     X.init = diag(p), W.init = diag(p))
    Sigihat <- quic_res$X
    Sighat <- quic_res$W
    ed1 <- NULL
  }else if (q == "addfrob"){
    if (missing(ed1)) {ed1 <- eigen(Sighat)}
    theta <- (ed1$values + sqrt(ed1$values^2 + 8*n*rho))/(2*n)
    Sighat <- ed1$vectors%*%diag(theta)%*%t(ed1$vectors)
  }else if (q == "multfrob") {
    if (missing(ed1)) {ed1 <- eigen(Sighat)}
    theta <- (ed1$values + sqrt(ed1$values^2 + 8*n*rho*n))/(2*n) # here, n = ||Sigi||_F^2 since we assume Sig = I
    Sighat <- ed1$vectors%*%diag(theta)%*%t(ed1$vectors)
  }else {
    stop("q must be 1, 1_off addfrob, or multfrob.")
  }
  
  # to ensure symmetry
  Sighat <- 1/2*(Sighat + t(Sighat))
  
  if(is.nan(sum(Sighat)) || is.na(sum(Sighat))) {maxit <- 0} # error checking
  
  like <- NA
  if(trace){
    like[iter] <- MVNloglike(x = x, mu = muhat, Sig = Sighat, rho = rho, q = q, print = T)
  }
  
  # repeat until EM converges
  ind <- 1
  xhat <- t(t(xc) + muhat) # undo centering
  na_idx <- is.na(x)
  while (ind > thr.em & iter < maxit) {
    oldxh <- xhat
    oldxc <- scale(oldxh, center = T, scale = F)
    
    iter <- iter + 1
    if (cov.corr) {covc <- matrix(0,p,p)}
    
    #E step
    for(i in rmi) {
      mi <- na_idx[i,]
      sinvs <- t(solve(Sighat[!mi,!mi], Sighat[!mi,mi]))
      xhat[i,mi] <- muhat[mi] + sinvs %*% (xhat[i,!mi] - muhat[!mi])
      if(cov.corr){ # compute corrected covariance matrix (see GA supplemental materials pg2)
        covc[mi,mi] <- covc[mi,mi] + Sighat[mi,mi] - sinvs %*% Sighat[!mi,mi]
      }          
    }
    
    #M Step
    if (n > 1) {muhat <- colMeans(xhat)} # update column means
    xhc <- t(t(xhat) - muhat) # center data
    if (cov.corr){ Sighat <- t(xhc)%*%xhc + covc}  else { Sighat <- t(xhc)%*%xhc }
    if(q == 1) {
      cat(paste0('Glasso Iteration: ', iter, '\n'))
      quic_res <- QUIC(S = Sighat/n, rho = rho/n, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigihat, W.init = Sighat)
      Sigihat <- quic_res$X
      Sighat <- quic_res$W
    }else if (q == "1_off"){
      cat(paste0('Glasso Iteration: ', iter, '\n'))
      tmp_p <- nrow(Sighat)
      rho_off <- rho/n * (matrix(1, nrow = tmp_p, ncol = tmp_p) - diag(tmp_p))
      quic_res <- QUIC(S = Sighat/n, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigihat, W.init = Sighat)
      Sigihat <- quic_res$X
      Sighat <- quic_res$W
    }else if (q == "addfrob"){
      ed <- eigen(Sighat)
      theta <- (ed$values + sqrt(ed$values^2 + 8*n*rho))/(2*n)
      Sighat <- ed$vectors%*%diag(theta)%*%t(ed$vectors)
    }else if (q == "multfrob") {
      ed <- eigen(Sighat)
      theta <- (ed$values + sqrt(ed$values^2 + 8*n*rho*n))/(2*n) # here, n = ||Sigi||_F^2 since we assume Sig = I
      Sighat <- ed$vectors%*%diag(theta)%*%t(ed$vectors)
    }else {
      stop("q must be 1, addfrob, or multfrob.")
    }
    
    # to ensure symmetry
    Sighat <- 1/2*(Sighat + t(Sighat))
    
    if(trace) {
      like[iter] <- MVNloglike(x = x, mu = muhat, Sig = Sighat, rho = rho, q = q, print = T) 
    }
    
    if (sum((oldxc[na_idx])^2, na.rm = T) == 0) {
      scaled_denom <- 1e-16 # to avoid dividing by 0
    }else {
      scaled_denom <- sum((oldxc[na_idx])^2,na.rm=TRUE)
    }
    
    ind <- sum((xhat[na_idx] - oldxh[na_idx])^2)/scaled_denom
    if(is.nan(sum(Sighat)) || is.na(sum(Sighat))){ind <- 0}
  }
  
  if(!trace){like=NULL}
  return(list(xhat=xhat,Sig=Sighat,mu=muhat,loglike=like,rmvn.ed1=ed1))
}
