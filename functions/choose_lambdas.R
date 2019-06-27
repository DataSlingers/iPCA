choose_lambdas <- function(dat, q, lams, lam_grid, prop_miss = 0.05, trcma = T, maxit = 20,
                           greedy.search = T, maxit.search = 1, save.filepath, save.fileext = "", seed = 100) {
  # choose_lambdas: function to choose penalty parameters via ECM-type algorithm
  #
  # Inputs:
  #  -dat = list of K data matrices (centered or uncentered)
  #  -q = type of penalty, either 1, addfrob, or multfrob
  #  -lams = (optional) vector of possible lams to use for construction of lam_grid
  #  -lam_grid = grid of penalty parameters to try; each row is a set of penalty parameters to try;
  #              lamS (if needed) must be in the first column, followed by lamDs
  #  -prop_miss = proportion of data missing from each Xk (default = 0.05)
  #  -trcma = T/F; whether or not to use TRCMAimpute_iPCA method (instead of RMNimpute_iPCA)
  #  -maixt = maximum number of iterations for the imputation method
  #  -greedy.search = logical T/F; whether or not to do a greedy search to choose penalty parameters
  #  -maxit.search = number of iterations for the greedy search
  #  -save.filepath = path where to save files (optional)
  #  -save.fileext = extension name of files (optional)
  #  -seed = numeric for set.seed (optional)
  #  
  # Outputs:
  #  -best_lambdas = optimal set of penalty parameters
  #  -errs = frobenius error from the scattered missing data
  #  -errs_tab = data.frame with errors per dataset
  #  -Xmiss = X with scattered missing elements
  
  # if (!missing(seed)) {
    set.seed(seed)
  # }
  
  Xs <- dat
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  K <- length(Xs)
  
  # leave out scattered missing elements in each of the Xks
  Xmiss <- Xs
  for (k in 1:K) {
    Xk <- Xs[[k]]
    miss_idx <- sample(1:length(Xk), prop_miss*length(Xk), replace = F)
    Xmiss[[k]][miss_idx] <- NA
  }
  
  if (!missing(save.filepath)) {
    saveRDS(Xmiss, file = paste0(save.filepath, "Xmiss", save.fileext, ".rds"))
    saveRDS(lam_grid, file = paste0(save.filepath, "lam_grid", save.fileext, ".rds"))
  }
  
  if (!missing(lam_grid)) {
    err <- rep(NA,nrow(lam_grid))
    errs <- matrix(NA, nrow = nrow(lam_grid), ncol = K)
    Xdiff <- list()
  }else {
    err <- NA
    errs <- matrix(NA, nrow = 1, ncol = K)
    Xdiff <- list()
  }

  if (!greedy.search) {
    if (!missing(lams)) {
      nreps <- ifelse(q == "multfrob", K, K+1)
      lam_grid <- expand.grid(rep(list(lams),nreps))
    }
    for (l in 1:nrow(lam_grid)) {
      print(l)
      lambdas <- as.matrix(lam_grid[l,])
      if (q == "multfrob" & ncol(lam_grid) == K) {
        lamS <- NULL
        lamDs <- lambdas
      }else if ((q == 1 | q == "1_off" | q == "addfrob") & ncol(lam_grid) == K+1) {
        lamS <- lambdas[1]
        lamDs <- lambdas[-1]
      }else {
        stop("q must be either 1, 1_off, addfrob, or multfrob, and check the number of penalty parameters.")
      }
      
      if (trcma) {
        if (l == 1) {
          ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
          rmvn.ed1s <- ans$rmvn.ed1s
        }else {
          ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, rmvn.ed1s = rmvn.ed1s, seed = seed, maxit = maxit)
        }
      }else {
        ans <- RMNimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
      }
      
      Xhat <- ans$xhat
      
      Xdiff[[l]] <- mapply(X = Xs, Y = Xhat, FUN = function(X,Y) {return(X-Y)},
                           SIMPLIFY = F)
      if (!missing(save.filepath)) {
        saveRDS(Xdiff, file = paste0(save.filepath, "Xdiff", save.fileext, ".rds"))
      }
      
      Xnum <- mapply(X = Xdiff[[l]], I = Xmiss, FUN = function(X, I) { return(X[is.na(I)]) }, SIMPLIFY = F) # numerator
      Xden <- mapply(X = Xs, I = Xmiss, 
                     FUN = function(X, I) {
                       Xc <- scale(X, center = T, scale = F)
                       return(Xc[is.na(I)]) 
                       }, SIMPLIFY = F) # denominator
      
      errs[l,] <- mapply(X = Xnum, Y = Xden,
                         FUN = function(X, Y) {
                           return(as.numeric((sum(X^2)^2) / (sum(Y^2)^2)))
                         })
      
      err[l] <- sum(errs[l,])
      
      if (!missing(save.filepath)) {
        saveRDS(errs, file = paste0(save.filepath, "lambda_errors_table", save.fileext, ".rds"))
        saveRDS(err, file = paste0(save.filepath, "lambda_errors", save.fileext, ".rds"))
      }
    }
    
    best_lambdas <- lam_grid[which.min(err),]
    errs_tab <- errs
    errs <- err

  }else { # do greedy search
    if (!missing(lam_grid)) {
      lams <- list()
      best_lambdas <- rep(NA,ncol(lam_grid))
      for (col in 1:ncol(lam_grid)) {
        lams[[col]] <- unique(lam_grid[,col]) # get unique penalty parameters
        best_lambdas[col] <- lams[[col]][1] # initialize
      }
    }else {
      nreps <- ifelse(q == "multfrob", K, K+1)
      lams <- rep(list(lams),nreps)
      best_lambdas <- rep(NA, nreps)
      for (col in 1:nreps) {
        best_lambdas[col] <- lams[[col]][1] # initialize
      }
    }
    
    # fix all but one penalty parameter
    search.iter <- 1; ptr <- 1
    already_searched <- NULL
    while (search.iter <= maxit.search) {
      old_best_lambdas <- best_lambdas
      for (idx in 1:length(best_lambdas)) {
        err <- rep(NA,length(lams[[idx]]))
        for (l in 1:length(lams[[idx]])) {
          print(paste0("Search.iter = ", search.iter, " // idx = ", idx, " // l = ",l))
          lambdas <- best_lambdas; lambdas[idx] <- lams[[idx]][l]
            
          if (q == "multfrob" & length(best_lambdas) == K) {
            lamS <- NULL
            lamDs <- lambdas
          }else if ((q == 1 | q == "1_off" | q == "addfrob") & length(best_lambdas) == K+1) {
            lamS <- lambdas[1]
            lamDs <- lambdas[-1]
          }else {
            stop("q must be either 1, addfrob, or multfrob, and check the number of penalty parameters.")
          }
          
          if (is.null(already_searched)) {
            marker <- T
          }else {
            marker <- sum(apply(X = as.matrix(already_searched), 1, FUN = function(X) {return(identical(X[-1], lambdas))}) == 0)
          }
          
          if (marker) { # if havent already looked at this choice of lambda
            if (trcma) {
              if (search.iter == 1 & idx == 1 & l == 1) {
                ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
                rmvn.ed1s <- ans$rmvn.ed1s
              }else {
                ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, rmvn.ed1s = rmvn.ed1s, seed = seed, maxit = maxit)
              }
            }else {
              ans <- RMNimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
            }
            
            Xhat <- ans$xhat
            
            Xdiff[[ptr]] <- mapply(X = Xs, Y = Xhat, FUN = function(X,Y) {return(X-Y)}, 
                                   SIMPLIFY = F)
            
#             if (!missing(save.filepath)) {
#               saveRDS(Xdiff, file = paste0(save.filepath, "Xdiff", save.fileext, ".rds"))
#             }
            
            Xnum <- mapply(X = Xdiff[[ptr]], I = Xmiss, FUN = function(X, I) { return(X[is.na(I)]) }, SIMPLIFY = F) # numerator
            Xden <- mapply(X = Xs, I = Xmiss, 
                           FUN = function(X, I) { 
                             Xc <- scale(X, center = T, scale = F)
                             return(Xc[is.na(I)]) 
                             }, SIMPLIFY = F) # denominator
            current_errs <- mapply(X = Xnum, Y = Xden,
                                   FUN = function(X, Y) {
                                     return(as.numeric((sum(X^2)^2) / (sum(Y^2)^2)))
                                   })
            errs <- rbind(errs, current_errs)
            err[l] <- sum(current_errs)
            
            already_searched <- rbind(already_searched, c(err[l], lambdas))
            ptr <- ptr + 1
          }else {
            err_idx <- which(apply(X = as.matrix(already_searched), 1, FUN = function(X) {return(identical(X[-1], lambdas))}))
            err[l] <- already_searched[err_idx, 1]
          }
          
#           if (!missing(save.filepath)) {
#             saveRDS(errs, file = paste0(save.filepath, "lambda_errors_table", save.fileext, ".rds"))
#           }
          
          if (l > 1) {
            if (err[l] > err[l-1]) {
              break
            }
          }
        }
        
        # update best lambdas
        best_lambdas[idx] <- lams[[idx]][which.min(err)]
        print(best_lambdas)
      }
      
      search.iter <- search.iter + 1
      
      if (identical(old_best_lambdas, best_lambdas)) { # no change
        search.iter <- maxit.search + 1 # exit while loop
      }
    }
    
    lam_grid <- already_searched[,-1]
    errs_tab <- errs
    errs <- rowSums(errs)

  }
  
  return(list(best_lambdas = best_lambdas, errs = errs, errs_tab = errs_tab, lam_grid = lam_grid, Xmiss = Xmiss))
}
