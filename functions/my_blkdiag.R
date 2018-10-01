my_blkdiag <- function(n_blks, rhos, p) {
  
  # my_blkdiag
  # 
  # inputs:
  # -n_blks = integer; number of blocks
  # -rhos = n_blks-vectors with correlations of each block
  # -p = number of rows (columns) in Delta
  # 
  # outputs:
  # -Delta = pxp SPARSE block diagonal matrix
  
  blks <- list()
  blk_p <- ceiling(p/n_blks)
  for (blk in 1:n_blks) {
    blks[[blk]] <- matrix(rhos[blk], nrow = blk_p, ncol = blk_p) + 
      diag(x = 1-rhos[blk], nrow = blk_p, ncol = blk_p)
  }
  Delta <- bdiag(blks)
  Delta <- Delta[1:p, 1:p]
  return(Delta)
}