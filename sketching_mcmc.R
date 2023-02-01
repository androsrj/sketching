library(mvtnorm)
library(fields)

# SOURCES
source("sim_data.R") # Spatial data simulation
source("pred_process.R") # Predictive process BFE
source("mcmc_functions.R") # Log likelihood/priors/posterior

## Define X and beta if desired
nTrain <- 500
nTest <- round(0.01 * nTrain)
n <- nTrain + nTest
p <- 5
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
trueBeta <- runif(p, 0, 10)

# Define distribution bounds for uniform prior on theta
a <- 1
b <- 5

# Simulate Y_star and W_star (having nStar observations)
fullData <- simSpatialData(n = n, X = X, beta = trueBeta)

# Extract training data
xTrain <- X[1:nTrain, ]
yTrain <- fullData$Y[1:nTrain]
wTrain <- fullData$W[1:nTrain]
sTrain <- fullData$S[1:nTrain]

# Extract testing data
xTest <- X[(nTrain+1):n, ]
yTest <- fullData$Y[(nTrain+1):n]
wTest <- fullData$W[(nTrain+1):n]
sTest <- fullData$S[(nTrain+1):n]

#subsetData <- simSpatialData(n = n + n_pred)
#Wstar <- subsetData$W
#Ystar <- subsetData$Y
#Sstar <- subsetData$S

# Fit predictive process BFE on training data
nKnots <- 20
#BFE <- predProcess(Wstar, Ystar, Sstar, n = n, X = X, beta = trueBeta)
#y <- fullData$Y
#condExpW <- fullData$condExpW
#S <- fullData$S

# Dimensions for subset data and distance matrix for full/subset data
DTrain <- rdist(sTrain)
DTest <- rdist(sTest)
DCov <- rdist(sTrain, sTest)

# Generate phi and compress data
m <- round(0.01 * nTrain)
phi <- matrix(rnorm(m * nTrain, 0, nTrain), nrow = m, ncol = nTrain)
newY <- phi %*% yTrain
newX <- phi %*% xTrain

# MCMC chain properties
nBurn <- 1000 # 3 to 4 thousand ideally
nThin <- 2
nIter <- nBurn + 10000 # 15 to 20 thousand ideally
sd <- 2 # (for proposal distributions, will need to tune)

trSigma2 <- trTau2 <- trTheta <- numeric(nIter) # Transformed parameters
beta <- matrix(0, nrow = p, ncol = nIter) # Beta
acceptSigma2 <- acceptTau2 <- acceptTheta <- 0 # Track acceptance rates

# Initial values of transformed parameters (except for beta, not transformed)
# (could definitely be different than what I have here)
trSigma2[1] <- log(1)
trTau2[1] <- log(1)
trTheta[1] <- log((2 - a) / (b - 2))
beta[ , 1] <- rep(0, p)

# Run Gibbs/Metropolis for one chain
for (i in 2:nIter) {
  
  ### Metropolis update (sigma2) ###
  
  propTrSigma2 <- rnorm(1, mean = trSigma2[i - 1], sd = sd)
  MHratio <- logPost(propTrSigma2, trTau2[i - 1], trTheta[i - 1], beta[ , i - 1]) - 
    logPost(trSigma2[i - 1], trTau2[i - 1], trTheta[i - 1], beta[ , i - 1])
  if(runif(1) < exp(MHratio)) {
    trSigma2[i] <- propTrSigma2
    acceptSigma2 <- acceptSigma2 + 1
  } else {
    trSigma2[i] <- trSigma2[i - 1]
  }
  
  ### Metropolis update (tau2) ###
  
  propTrTau2 <- rnorm(1, mean = trTau2[i - 1], sd = sd)
  MHratio <- logPost(trSigma2[i], propTrTau2, trTheta[i - 1], beta[ , i - 1]) - 
    logPost(trSigma2[i], trTau2[i - 1], trTheta[i - 1], beta[ , i - 1])
  if(runif(1) < exp(MHratio)) { 
    trTau2[i] <- propTrTau2
    acceptTau2 <- acceptTau2 + 1
  } else {
    trTau2[i] <- trTau2[i - 1]
  }
  
  ### Metropolis update (theta) ###
  
  propTrTheta <- rnorm(1, mean = trTheta[i - 1], sd = sd)
  MHratio <- logPost(trSigma2[i], trTau2[i], propTrTheta, beta[ , i - 1]) - 
    logPost(trSigma2[i], trTau2[i], trTheta[i - 1], beta[ , i - 1])
  if(runif(1) < exp(MHratio)) {
    trTheta[i] <- propTrTheta
    acceptTheta <- acceptTheta + 1 
  } else {
    trTheta[i] <- trTheta[i - 1]
  }
  
  ### Gibbs update (beta) ###
  
  # Recalculate C and Cstar based on most recent parameter values
  tempSigma2 <- exp(trSigma2[i])
  tempTau2 <- exp(trTau2[i])
  tempTheta <- fInv(trTheta[i])
  C <- tempSigma2 * exp(- tempTheta * DCov)
  Cstar <- tempSigma2 * exp(- tempTheta * DTest)
  Sigma <- phi %*% C %*% solve(Cstar) %*% t(C) %*% t(phi) + tempTau2 * diag(m)
  SigmaInv <- solve(Sigma)
  SigmaBeta <- (n / m) * t(newX) %*% SigmaInv %*% newX + diag(p)
  SigmaBetaInv <- solve(SigmaBeta)
  meanBeta <- (n / m) * SigmaBetaInv %*% t(newX) %*% SigmaInv %*% newY
  beta[ , i] <- t(rmvnorm(1, meanBeta, SigmaBetaInv))
}

# Acceptance rates (for Metropolis-sampled parameters)
acceptSigma2 / nIter
acceptTau2 / nIter
acceptTheta / nIter

# Remove burn-in and perform thinning
index <- seq(nBurn + 1, nIter, by = nThin)
trSigma2 <- trSigma2[index]
trTau2 <- trTau2[index]
trTheta <- trTheta[index]
beta <- beta[ , index]
nSamples <- length(index)

# Back-transform
sigma2 <- exp(trSigma2)
tau2 <- exp(trTau2)
theta <- fInv(trTheta)

# Trace plots
plot(1:nSamples, sigma2, type = 'l', main = "Trace Plot for sigma2")
plot(1:nSamples, tau2, type = 'l', main = "Trace Plot for tau2")
plot(1:nSamples, theta, type = 'l', main = "Trace Plot for theta")
plot(1:nSamples, beta[1, ], type = 'l', main = "Trace Plot for beta_1")
plot(1:nSamples, beta[p, ], type = 'l', main = "Trace Plot for beta_p")

# Posterior mean estimates
mean(sigma2)
mean(tau2)
mean(theta)
apply(beta, 1, mean)

# 95% credible intervals
quantile(sigma2, c(0.025, 0.975))
quantile(tau2, c(0.025, 0.975))
quantile(theta, c(0.025, 0.975))
apply(beta, 1, quantile, c(0.025, 0.975))





