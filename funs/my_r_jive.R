# my version of the jive() function from library(r.jive)
# only difference is that my version can change the "nperms" parameter
library(SpatioTemporal)
my_r_jive <- function (data, rankJ = 1, rankA = rep(1, length(data)), method = "perm", 
          dnames = names(data), conv = "default", maxiter = 1000, nperms,
          scale = TRUE, center = TRUE, orthIndiv = TRUE, est = TRUE, 
          showProgress = TRUE) 
{
  l <- length(data)
  n <- c()
  d <- c()
  for (i in 1:length(data)) {
    n[i] <- nrow(data[[i]]) * ncol(data[[i]])
    d[i] <- nrow(data[[i]])
  }
  for (i in 1:l) {
    temp <- SVDmiss(data[[i]], ncomp = min(ncol(data[[i]]), 
                                           nrow(data[[i]])))[[1]]
    data[[i]] <- temp$u %*% diag(x = temp$d) %*% t(temp$v)
  }
  centerValues <- list()
  scaleValues <- c()
  for (i in 1:l) {
    if (center) {
      centerValues[[i]] <- apply(data[[i]], 1, mean, na.rm = T)
      data[[i]] <- data[[i]] - matrix(rep(centerValues[[i]], 
                                          ncol(data[[i]])), nrow = nrow(data[[i]]))
    }
    if (!center) {
      centerValues[[i]] <- rep(0, d[i])
    }
    if (scale) {
      scaleValues[i] <- norm(data[[i]], type = "f") * 
        sqrt(sum(n))
      data[[i]] <- data[[i]]/scaleValues[i]
    }
    if (!scale) {
      scaleValues[i] <- 1
    }
  }
  if (conv == "default") {
    conv = 10^(-6) * norm(do.call(rbind, data), type = "f")
  }
  orig <- data
  if (method == "given") {
    if (est) {
      u <- list()
      for (i in 1:l) {
        if (nrow(data[[i]]) > ncol(data[[i]])) {
          temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                             nv = ncol(data[[i]]))
          data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                            nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
          u[[i]] <- temp$u
        }
        else {
          u[[i]] <- diag(1, nrow(data[[i]]))
        }
      }
    }
    if (showProgress) {
      cat("Running JIVE algorithm for ranks:\njoint rank:", 
          rankJ, ", individual ranks:", rankA, "\n")
    }
    temp <- jive.iter(data, rankJ, rankA, conv = conv, maxiter = maxiter, 
                      orthIndiv = orthIndiv, showProgress = showProgress)
    joint <- temp$joint
    individual <- temp$individual
    if (est) {
      for (i in 1:l) {
        joint[[i]] <- u[[i]] %*% joint[[i]]
        individual[[i]] <- u[[i]] %*% individual[[i]]
      }
    }
  }
  else if (method == "perm") {
    temp <- jive.perm(data, est = est, conv = conv, maxiter = maxiter, nperms = nperms,
                      orthIndiv = orthIndiv, showProgress = showProgress)
    joint <- temp$joint
    individual <- temp$individual
    rankJ <- temp$rankJ
    rankA <- temp$rankA
    converged <- temp$converged
  }
  else if (method == "bic") {
    if (est) {
      u <- list()
      for (i in 1:l) {
        if (nrow(data[[i]]) > ncol(data[[i]])) {
          temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                             nv = ncol(data[[i]]))
          data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                            nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
          u[[i]] <- temp$u
        }
        else {
          u[[i]] <- diag(1, nrow(data[[i]]))
        }
      }
    }
    temp <- bic.jive(data, n, d, conv = conv, maxiter = maxiter, 
                     orthIndiv = orthIndiv, showProgress = showProgress)
    joint <- temp$joint
    individual <- temp$individual
    rankJ <- temp$rankJ
    rankA <- temp$rankA
    bic.table <- temp$bic.table
    if (est) {
      for (i in 1:l) {
        joint[[i]] <- u[[i]] %*% joint[[i]]
        individual[[i]] <- u[[i]] %*% individual[[i]]
      }
    }
  }
  if (is.null(dnames)) {
    for (i in 1:l) {
      names(orig)[[i]] <- paste("Source", i, sep = "_")
    }
  }
  else {
    names(orig) <- dnames
  }
  result <- list(data = orig, joint = joint, individual = individual, 
                 rankJ = rankJ, rankA = rankA, method = method)
  if (method == "bic") {
    result$bic.table <- bic.table
  }
  if (method == "perm") {
    result$converged <- converged
  }
  result$scale <- list(center, scale, centerValues, scaleValues)
  names(result$scale) <- c("Center", "Scale", "Center Values", 
                           "Scale Values")
  class(result) <- "jive"
  return(result)
}
