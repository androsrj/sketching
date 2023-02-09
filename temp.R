source("mcmc_functions.R")

trueBeta
trueSigma2 <- 2
trueTau2 <- 7
trueTheta <- 3

f <- function(x, a=1, b=5) {
  log((a-x)/(x-b))
}

sigma2Vals <- seq(0.5, 5, by=0.1)
sigma2Density <- sapply(sigma2Vals, \(x) logPost(log(x), log(trueTau2), trueTheta, trueBeta))
plot(sigma2Vals, sigma2Density)
sigma2Vals[which.max(sigma2Density)]

tau2Vals <- seq(0.5, 5, by=0.1)
tau2Density <- sapply(tau2Vals, \(x) logPost(log(trueSigma2), log(x), trueTheta, trueBeta))
plot(tau2Vals, tau2Density)
tau2Vals[which.max(tau2Density)]

thetaVals <- seq(1.1, 4.9, by=0.1)
thetaDensity <- sapply(thetaVals, \(x) logPost(log(trueSigma2), log(trueTau2), f(x), trueBeta))
plot(thetaVals, thetaDensity)
thetaVals[which.max(thetaDensity)]


