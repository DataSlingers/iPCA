sim_rand_cov <- function(num, p) {
  # simulate random covariance matrix
  
  if (missing(num)) {
    num <- sample(1:4, 1)
  }
  
  if (num == 1) {
    cov <- my_toeplitz(runif(n=1, min=-.9, max=.9), p)
  }else if (num == 2) {
    n_blks <- sample(3:10,1)
    cov <- my_blkdiag(n_blks=n_blks, rhos=runif(n=n_blks, min=0, max=.9), p=p)
  }else if (num == 3) {
    U <- randortho(n=p, type="orthonormal")
    r <- sample(5:min(50,p),1) # rank
    d <- c(runif(n=r, min=5, max=75), rep(1,p-r))
    cov <- U %*% diag(d) %*% t(U)
  }else if (num == 4) {
    dat <- log(read.csv("../data/TCGA_ov_mirna.csv", header = F))
    cov <- my_realcov(dat, 25, p)
  }

  return(cov)
}