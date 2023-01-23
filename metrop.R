source("mcmc_functions.R")
library(mcmc)

### LOG POSTERIOR ###
obj <- function(params) {
  sigma2 <- params[1]
  tau2 <- params[2]
  theta <- params[3]
  beta <- params[4:length(params)]
  logLik(sigma2, tau2, theta, beta) + # Likelihood
    logPriorBeta(beta) + logPriorTheta(theta) + # Priors
    logPriorTau2(tau2) + logPriorSigma2(sigma2) # Priors
}

p <- 5
m <- metrop(obj = obj, initial = c(0.5, 0.5, 2, rep(0, p)), nbatch = 2, blen = 100)
m$batch
m$final
m$accept

