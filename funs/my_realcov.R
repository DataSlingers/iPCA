my_realcov <- function(dat, top_p, p) {
  
  # my_realcov
  # 
  # inputs:
  # -dat = real data matrix
  # -top_p = integer; number of columns with highest variance to keep
  # -p = number of rows (columns) in Delta
  # 
  # outputs:
  # -Delta = pxp covariance matrix from real data sim
      
  ordered <- sort(x = apply(dat, 2, var), decreasing = T, index.return = T)
  dat_e <- dat[,-ordered$ix[1:top_p]]
  tmp <- cbind(dat[,ordered$ix[1:top_p]], 
               dat_e[,sample(1:ncol(dat_e), p-top_p, replace = F)])
  tmp <- as.matrix(tmp)
  Delta <- t(tmp) %*% tmp / p + diag(p)*.1
  return(Delta)
}