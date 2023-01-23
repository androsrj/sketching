library(mvtnorm)

### LOG PRIORS ### 

# Sigma2 (inverse gamma)
logPriorSigma2 <- function(sigma2, a = 1, b = 2) {
  log(b^a) - lgamma(a) - (a - 1) * sigma2 - b / sigma2
}

# Tau2 (inverse gamma)
logPriorTau2 <- function(tau2, a = 1, b = 2) {
  log(b^a) - lgamma(a) - (a - 1) * tau2 - b / tau2
}

# Theta (uniform for now, could try discrete later)
logPriorTheta <- function(theta, a = 1, b = 3) {
  dunif(theta, a, b, log = TRUE)
}

# Beta (standard MV normal)
logPriorBeta <- function(beta) {
  p <- length(beta)
  dmvnorm(beta, mean = rep(0, p), sigma = diag(p), log = TRUE)
}

# Inverse transformation for theta, where trTheta = log((theta - a) / (b - theta))
fInv <- function(trTheta, a = 1, b = 3) {
  (b * exp(trTheta) + a) / (1 + exp(trTheta))
}

# Log-Jacobian for theta, (log-derivative of fInv function above)
jac <- function(trTheta, a = 1, b = 3) {
  # log( (b - a) * exp(trTheta) / (1 + exp(trTheta))^2 ) $ or simplify, as below
  log(b - a) + trTheta - 2 * log(1 + exp(trTheta))
}

### LOG LIKELIHOOD ###

logLik <- function(sigma2, tau2, theta, beta) {
  p <- length(beta)
  C <- sigma2 * exp(- theta * Dcov)
  Cstar <- sigma2 * exp(- theta * Dstar)
  dmvnorm(as.vector(newY), 
          as.vector(newX %*% beta), 
          phi %*% C %*% solve(Cstar) %*% t(C) %*% t(phi) + tau2 * diag(n_star))
}

### LOG POSTERIOR ###
logPost <- function(trSigma2, trTau2, trTheta, beta) {
  logLik(exp(trSigma2), exp(trTau2), fInv(trTheta), beta) + # Likelihood
    logPriorBeta(beta) + logPriorTheta(fInv(trTheta)) + # Priors
    logPriorTau2(exp(trTau2)) + logPriorSigma2(exp(trSigma2)) + # Priors
    exp(trTau2) + exp(trSigma2) + jac(trTheta) # Jacobians
}






