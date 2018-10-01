my_toeplitz <- function(rho, p) {
  
  # my_toeplitz
  # 
  # inputs:
  #  -rho = scalar between 0 and 1
  #  -p = number of rows (columns) in Delta
  # output:
  #  -Delta = pxp toeplitz matrix with ij-entry = rho^(abs(i-j))
      
  row1 <- rho^(0:(p-1))
  Delta <- toeplitz(row1)
  return(Delta)
}