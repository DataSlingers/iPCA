sim_from_rand_dist <- function(num, n, p) {
  # simulate from a random distribution and output n x pk matrix
  
  if (missing(num)) {
    num <- sample(1:4, 1)
  }
  
  if (num == 1) {
    dist <- rnorm(n*p)
  }else if (num == 2) {
    dist <- sample(-2:2, n*p, replace = T)
  }else if (num == 3) {
    dist <- runif(min=0, max=1, n=n*p)
  }else if (num == 4) {
    dist <- rexp(n=n*p, rate=1)
  }
  out <- matrix(dist, nrow = n, ncol = p)
  return(out)
}