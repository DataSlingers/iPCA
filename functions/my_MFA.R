# MFA: runs multiple factor analysis (J. Pages, 1994) i.e. multiblock PCA
#
# inputs:
#   dat = K-list of n x pk CENTERED data matrices
#   npc = number of PCs to compute (default = 10)
#
# outputs:
#   Sig = n x n row covariance matrix estimate
#   Delts = K-list of pk x pk column covariance matrix estimates
#   var_explained = eig from FactoMinR::MFA function 

library(FactoMineR)

my_MFA <- function(dat, npc = 10) {
  
  Xs <- dat
  
  pks <- mapply(Xs, FUN = ncol)
  K <- length(Xs)
  
  # center data
  for (k in 1:K) {
    Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
  }
  
  Xs_df <- data.frame(do.call(cbind, Xs))
  
  res <- MFA(base = Xs_df, group = pks, type = rep("c",K), ncp = npc, graph = F)
  U <- res$ind$coord; colnames(U) <- paste0("X",1:npc)
  Vs <- list(); V <- res$quanti.var$coord; ind <- 1
  for (k in 1:K) {
    Vs[[k]] <- V[ind:(ind + pks[k] - 1),]; colnames(Vs[[k]]) <- paste0("X",1:npc)
    ind <- ind + pks[k]
  }
  var_explained <- res$eig
  
  return(list(U = U, Vs = Vs, var_explained = var_explained))
}



