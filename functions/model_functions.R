## @knitr models

library(rARPACK) # eigs
library(expm) # sqrtm
library(rmutil) # rlaplace
library(pracma)

# cov_types can be "base", "identity", or a vector of types

ipca_base_sim <- function(K = 3, dim_U = 2, n_clust = 3, 
                          n = 150, pks = c(300, 500, 400), pks_add = 0, # for base sims
                          # snr = 3, Dsig = n/snr, Dks = c(snr*pks[1], 3*snr*pks[2], 2*snr*pks[3]),
                          snr = 3, Dsig = 50, Dks = c(900, 4500, 2400),
                          cov_types = 'base', cov_params = NA,  # cov_types = "random" for vary_K simulation
                          err_type = 'none', err_params = NA, nsim) {
  
  pks <- pks + pks_add
  
  Xsims <- list()
  
  for (sim in 1:nsim) {
    # simulate U (not necessarily orthogonal to begin with), Sigma, and ylabels
    sigu = 1; # influences sd of u
    if (dim_U == 2) {
      U_all <- matrix(c(1,5,5,3,3,1,-1,2,4,7,0,8,10,3), ncol = 2, byrow = T) # cluster centers for 2d case
      U_subspace <- U_all[1:n_clust, 1:dim_U] # get cluster centers
    }else { # for d > 2
      U_subspace <- matrix(sample(x = 1:10, size = n_clust*dim_U, replace = T),
                           nrow = n_clust, ncol = dim_U)
    }
    
    # simulate Sigma
    uu <- do.call(rbind, replicate(ceiling(n/n_clust), U_subspace, simplify=FALSE))
    uu <- uu[1:n,] + matrix(rnorm(n*dim_U, 0, 1), nrow = n, ncol = dim_U)*sigu
    uu <- scale(uu, center = F, scale = sqrt(colSums(uu^2))) # normalize columns to norm 1
    y <- rep(1:n_clust, times = ceiling(n/n_clust))
    y <- y[1:n] # set cluster labels
    tt <- diag(n)
    UU <- cbind(uu, tt[,(dim_U+1):n])
    Dtrue <- c(Dsig, Dsig*.5, Dsig*.4, Dsig*.25, Dsig*.15)
    Duu <- diag(c(Dtrue[1:dim_U], rep(1,n-dim_U)))
    Sig <- UU %*% Duu %*% t(UU); Sig <- 1/2 * (Sig + t(Sig)) # ensure symmetry
    Sig_eigs <- eigs(Sig, dim_U, 'LM') # want true orthogonal components
    d <- Sig_eigs$values; U <- Sig_eigs$vectors
    ordered <- sort(x = d, decreasing = T, index.return = T)
    Du <- diag(ordered$x); U <- U[,ordered$ix] # order eigenvalues and vectors
    
    # simulate Delta_ks
    Deltks <- list()
    
    if (identical(cov_types, 'base')) {
      # toeplitz
      Deltx <- my_toeplitz(.9, pks[1])
      Delt_eigs <- eigen(Deltx)
      aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
      Deltks[[1]] <- aa %*% diag(dd*Dks[1]/max(dd)) %*% t(aa) # incorporate snr
      
      # real data cov
      dat <- log(read.csv("../data/TCGA_ov_mirna.csv", header = F))
      Delty <- my_realcov(dat, 25, pks[2])
      Delt_eigs <- eigen(Delty)
      aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
      Deltks[[2]] <- aa %*% diag(dd*Dks[2]/max(dd)) %*% t(aa) # incorporate snr
      
      # block diagonal
      Deltz <- my_blkdiag(5, c(.6, .4, .6, .2, .8), pks[3])
      Delt_eigs <- eigen(Deltz)
      aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
      Deltks[[3]] <- aa %*% diag(dd*Dks[3]/max(dd)) %*% t(aa) # incorporate snr
      
    }else if (identical(cov_types, "random")) {
      pks <- sample(x=200:500, size=K, replace=T) # randomly choose pks
      rand_max <- 4
      for (k in 1:K) {
        rand_num <- sample(1:rand_max,1)
        if (rand_num == 4) {
          rand_max <- 3 # only allow real data covariance to be used once
        }
        Delt <- sim_rand_cov(num = rand_num, p = pks[k]) # simulate random deltk
        Delt_eig <- eigen(Delt)
        aa <- Delt_eig$vectors; dd <- Delt_eig$values
        Deltks[[k]] <- aa %*% diag(dd*runif(n=1, min=500, max=2500)/max(dd)) %*% t(aa)
      }
    }else { # if specify covariance types directly
      for (k in 1:K) {
        if (identical(cov_types[[k]], 'toeplitz')) {
          Delta <- my_toeplitz(cov_params[[k]][1], pks[k])
        }else if (identical(cov_types[[k]], 'blkdiag')) {
          Delta <- my_blkdiag(cov_params[[k]][1], cov_params[[k]][2], pks[k])
        }else if (identical(cov_types[[k]], 'realcov')) {
          Delta <- my_realcov(cov_params[[k]][1], cov_params[[k]][2], pks[k])
        }else {
          stop('Unknown covariance type.')
        }
        
        Delt_eigs <- eigen(Delta)
        aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
        Deltks[[k]] <- aa %*% diag(dd*Dks[k]/max(dd)) %*% t(aa) # incorporate snr
      }
    }
    
    # generate Xs
    Xs <- list()
    for (k in 1:K) {
      # generate data from matrix normal dist
      Sig_eig <- eigen(Sig); Sig_V <- Sig_eig$vectors; Sig_D <- Sig_eig$values
      Sig_sqrt <- Sig_V %*% sqrt(diag(Sig_D)) %*% t(Sig_V)
      Delt_eig <- eigen(Deltks[[k]]); Delt_V <- Delt_eig$vectors; Delt_D <- Delt_eig$values
      Delt_sqrt <- Delt_V %*% sqrt(diag(Delt_D)) %*% t(Delt_V)
      Xs[[k]] <- Sig_sqrt %*% matrix(rnorm(n*pks[k],0,1), nrow = n, ncol = pks[k]) %*% Delt_sqrt
      
      # add errors
      if (identical(err_type, 'laplace')) {
        Xs[[k]] <- Xs[[k]] + sqrt(err_params)*matrix(rlaplace(n*pks[k]), nrow = n, ncol = pks[k])
      }else if (identical(err_type, 't')) {
        Xs[[k]] <- Xs[[k]] + matrix(rt(n = n*pks[k], df = err_params), nrow = n, ncol = pks[k])
      }else if (!identical(err_type, 'none')) {
        stop('Unknown error type.')
      }
      
      # center data
      Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
    }
    
    truth <- list(Sig_true = Sig, Deltks_true = Deltks, y_true = as.factor(y))
    
    Xsims[[sim]] <- list(sim_data = Xs, truth = truth)
  }
  
  return(Xsims) # return centered data
}


# to vary number of dimensions of Sig (not a cluster model)
ipca_noncluster_model <- function(K = 3, dim_U = 2,
                                  n = 150, pks = c(300, 500, 400), pks_add = 0,
                                  snr = 3, Dsig = 50, Dks = c(900, 4500, 2400),
                                  cov_types = "base", cov_params = NA,
                                  err_type = "none", err_params = NA, nsim) {
  pks <- pks + pks_add
  
  Xsims <- list()
  
  for (sim in 1:nsim) {
    # simulate Sig
    U <- randortho(n=n, type="orthonormal") # simulate random orthogonal matrix (eigenvectors)
    d <- c(sort(runif(n=dim_U, min=5, max=75), decreasing=T), # make sure it has rank dim_U
           rep(1,n-dim_U)) 
    Sig <-U %*% diag(d) %*% t(U)
    if (!all.equal(Sig, t(Sig))) { # error checking
      stop("Sigma is not symmetric")
    }
    
    # simulate Delta_ks
    Deltks <- list()
    
    if (identical(cov_types, 'base')) {
      # toeplitz
      Deltx <- my_toeplitz(.9, pks[1])
      Delt_eigs <- eigen(Deltx)
      aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
      Deltks[[1]] <- aa %*% diag(dd*Dks[1]/max(dd)) %*% t(aa) # incorporate snr
      
      # real data cov
      dat <- log(read.csv("../data/TCGA_ov_mirna.csv", header = F))
      Delty <- my_realcov(dat, 25, pks[2])
      Delt_eigs <- eigen(Delty)
      aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
      Deltks[[2]] <- aa %*% diag(dd*Dks[2]/max(dd)) %*% t(aa) # incorporate snr
      
      # block diagonal
      Deltz <- my_blkdiag(5, c(.6, .4, .6, .2, .8), pks[3])
      Delt_eigs <- eigen(Deltz)
      aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
      Deltks[[3]] <- aa %*% diag(dd*Dks[3]/max(dd)) %*% t(aa) # incorporate snr
      
    }else { # if specifiy covariance types directly
      for (k in 1:K) {
        if (identical(cov_types[[k]], 'toeplitz')) {
          Delta <- my_toeplitz(cov_params[[k]][1], pks[k])
        }else if (identical(cov_types[[k]], 'blkdiag')) {
          Delta <- my_blkdiag(cov_params[[k]][1], cov_params[[k]][2], pks[k])
        }else if (identical(cov_types[[k]], 'realcov')) {
          Delta <- my_realcov(cov_params[[k]][1], cov_params[[k]][2], pks[k])
        }else {
          stop('Unknown covariance type.')
        }
        
        Delt_eigs <- eigen(Delta)
        aa <- Delt_eigs$vectors; dd <- Delt_eigs$values
        Deltks[[k]] <- aa %*% diag(dd*Dks[k]/max(dd)) %*% t(aa) # incorporate snr
      }
    }
    
    # generate Xs
    Xs <- list()
    Sig_sqrt <- U %*% sqrt(diag(d)) %*% t(U)
    for (k in 1:K) {
      # generate data from matrix normal dist
      Delt_eig <- eigen(Deltks[[k]]); Delt_V <- Delt_eig$vectors; Delt_D <- Delt_eig$values
      Delt_sqrt <- Delt_V %*% sqrt(diag(Delt_D)) %*% t(Delt_V)
      Xs[[k]] <- Sig_sqrt %*% matrix(rnorm(n*pks[k],0,1), nrow = n, ncol = pks[k]) %*% Delt_sqrt
      
      # add errors
      if (identical(err_type, 'laplace')) {
        Xs[[k]] <- Xs[[k]] + sqrt(err_params)*matrix(rlaplace(n*pks[k]), nrow = n, ncol = pks[k])
      }else if (identical(err_type, 't')) {
        Xs[[k]] <- Xs[[k]] + matrix(rt(n = n*pks[k], df = err_params), nrow = n, ncol = pks[k])
      }else if (!identical(err_type, 'none')) {
        stop('Unknown error type.')
      }

      # center data 
      Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
    }
    
    truth <- list(Sig_true = Sig, Deltks_true = Deltks)
    
    Xsims[[sim]] <- list(sim_data = Xs, truth = truth)
  }
  return(Xsims)
}
      


# jive model
jive_model <- function(K = 3, rankJ = 5, rankA = c(10,15,20),
                       n = 75, pks = c(100, 150, 200), 
                       noise = 1, model_type = "random", # model_type can be illustrative or random
                       nsim) {
  Xsims <- list()
  
  for (sim in 1:nsim) {
    if (model_type == "illustrative") { # small illustrative jive simulation
      K <- 2; n = 100; pks = c(50,50)
      
      J <- matrix(0, nrow = n, ncol = pks[1])
      j <- rnorm(n)
      J[,1:25] <- t(j)
      
      A1 <- matrix(rep(-2:2, each = 20), nrow = n, ncol = pks[1], byrow = F)
      a2 <- sample(-2:2, 100, replace = T)
      A2 <- matrix(a2, nrow = n, ncol = pks[2], byrow = F)
      As <- list(A1, A2)
      
      X1 <- J + A1 + matrix(rnorm(n*pks[1], mean = 0, sd = noise), nrow = n, ncol = pks[1])
      X2 <- J + A2 + matrix(rnorm(n*pks[2], mean = 0, sd = noise), nrow = n, ncol = pks[2])
      Xs <- list(X1, X2)
      
      J <- cbind(J,J)
      
    }else if (model_type == "random") {
      p <- sum(pks)
      
      S <- sim_from_rand_dist(n = n, p = rankJ)
      U <- sim_from_rand_dist(n = rankJ, p = p)
      J <- S %*% U
      
      Xs <- As <- list()
      idx <- 1
      for (k in 1:K) {
        rankAk <- rankA[k]
        Sk <- sim_from_rand_dist(n = n, p = rankAk)
        Wk <- sim_from_rand_dist(n = rankAk, p = pks[k])
        Ak <- Sk %*% Wk
        As[[k]] <- Ak
        
        col_idx <- idx:(idx + pks[k] - 1)
        Xs[[k]] <- J[,col_idx] + Ak + 
          matrix(rnorm(n = n*pks[k], mean = 0, sd = noise), nrow = n, ncol = pks[k])
        idx <- idx + pks[k]
      }
    }else {
      stop("Unknown model_type. model_type must be either illustrative or random")
    }
    
    Sig <- J %*% t(J)
    Deltks <- lapply(As, FUN = function(A) {return(t(A) %*% A)})                
    
    truth <- list(Sig_true = Sig, Deltks_true = Deltks, J_true = J, As_true = As, rankJ = rankJ, rankA = rankA)
    Xsims[[sim]] <- list(sim_data = Xs, truth = truth)
  }
  
  return(Xsims)
}


real_data_sim <- function(dat) { # for rosmap data
  Xsims <- list()
  Xs <- dat
  for (k in 1:length(dat)) {
    # center data
    Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
  }
  
  Xsims[[1]] <- list(sim_data = Xs)
  
  return(Xsims) # return centered data
}
