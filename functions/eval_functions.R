## metrics

subspace_recovery <- function(Sig, Sigh, Uh, dim_U) {
  # Sig = true Sigma
  # Sigh = estimated Sigma
  # Uh = estimated eigenvectors of Sigma (only required if Sigh is missing)
  # dim_U = dimension of true subspace of Sigma
  # must specify either Sigh or Uh
  
  if (!missing(Sigh)) {
    Sigh_eigs <- eigs(Sigh, dim_U, which = "LM")
    Uh <- Sigh_eigs$vectors # eigenvectors of Sigh
  }else {
    Uh <- Uh[,1:dim_U]
  }
  
  Sig_eigs <- eigs(Sig, dim_U, which = "LM")
  U <- Sig_eigs$vectors # truth
  
  sr_metric <- base::norm(U %*% solve(t(U) %*% U) %*% t(U) - Uh %*% solve(t(Uh) %*% Uh) %*% t(Uh), "F")^2
  return(1/dim_U * sr_metric)
}

# for the following metric, see "Signal-plus-noise matrix models: eigenvector deviations and fluctuations"
# smaller subspace_angle is better
# note: arccos undefined for values outside [-1,1] so doens't work for non-orthogonal matrices
# canonical angles with spectral norm
subspace_angle_spect <- function(Sig, Sigh, Uh, dim_U) {
  # Sig = true Sigma
  # Sigh = estimated Sigma
  # Uh = estimated eigenvectors of Sigma (only required if Sigh is missing)
  # dim_U = dimension of true subspace of Sigma
  # must specify either Sigh or Uh
  
  if (!missing(Sigh)) {
    Sigh_eigs <- eigs(Sigh, dim_U, which = "LM")
    Uh <- Sigh_eigs$vectors # eigenvectors of Sigh
  }else {
    Uh <- Uh[,1:dim_U]
  }
  
  Sig_eigs <- eigs(Sig, dim_U, which = "LM")
  U <- Sig_eigs$vectors # truth
  
  u_svd <- svd(t(Uh) %*% U)
  canonical_angles <- diag(acos(diag(u_svd$d)))
  
  angle_metric <- base::norm(sin(canonical_angles), '2')
  return(1/dim_U * angle_metric)
}

# canonical angles with frobenius norm
subspace_angle_frob <- function(Sig, Sigh, Uh, dim_U) {
  # Sig = true Sigma
  # Sigh = estimated Sigma
  # Uh = estimated eigenvectors of Sigma (only required if Sigh is missing)
  # dim_U = dimension of true subspace of Sigma
  # must specify either Sigh or Uh
  
  if (!missing(Sigh)) {
    Sigh_eigs <- eigs(Sigh, dim_U, which = "LM")
    Uh <- Sigh_eigs$vectors # eigenvectors of Sigh
  }else {
    Uh <- Uh[,1:dim_U]
  }
  
  Sig_eigs <- eigs(Sig, dim_U, which = "LM")
  U <- Sig_eigs$vectors # truth
  
  u_svd <- svd(t(Uh) %*% U)
  canonical_angles <- diag(acos(diag(u_svd$d)))
  
  angle_metric <- base::norm(sin(canonical_angles), 'F')
  return(1/dim_U * angle_metric)
}
