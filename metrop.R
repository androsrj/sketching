source("mcmc_functions.R")
library(mcmc)

### LOG POSTERIOR ###
obj <- function(params) {
  trSigma2 <- params[1]
  trTau2 <- params[2]
  trTheta <- params[3]
  beta <- params[4:length(params)]
  logLik(exp(trSigma2), exp(trTau2), fInv(trTheta), beta) + # Likelihood
    logPriorBeta(beta) + logPriorTheta(fInv(trTheta)) + # Priors
    logPriorTau2(exp(trTau2)) + logPriorSigma2(exp(trSigma2)) + # Priors
    trTau2 + trSigma2 + jac(trTheta) # Jacobians
}

p <- 5
met <- metrop(obj = obj, initial = c(log(1), log(1), fInv(3), trueBeta), nbatch = 4, blen=10000)
met$batch
met$final
met$accept

