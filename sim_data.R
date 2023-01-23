library(mvtnorm)
library(fields)

simSpatialData <- function(n, X = NULL, range = c(0, 10), dims = 2, 
                           theta = 5, sigma2 = 3, beta = c(2,2), tau2 = 2,
                           covariance = "exponential") {
  
  # Sample the locations and put them into an n-by-dims matrix
  locations <- runif(n * dims, range[1], range[2])
  s <- matrix(locations, nrow = n, ncol = dims)
  
  # Compute the covariance matrix (symmetric)
  if (covariance == "exponential") {
    D <- rdist(s)
  } else if (covariance == "exp_squared") {
    D <- rdist(s)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  C <- sigma2 * exp(- theta * D)
  
  # Sample W
  W <- t(rmvnorm(1, sigma = C))
  
  # If X is not supplied, set it to a matrix containing a column of ones (intercept) 
  # and one standard normal predictor
  if (is.null(X)) {
    Xint <- rep(1, n)
    Xpred <- rnorm(n)
    X <- matrix(Xint, Xpred, nrow = n, ncol = 2)
  }
  
  # Sample epsilon
  eps <- rnorm(n, 0, 1 / tau2)
  
  # Generate Y
  Y <- X %*% beta + W + eps
  
  # Return data
  return(list(Y = as.vector(Y), W = as.vector(W), S = s))
}






