# This is a general sandbox / experimental script for testing 
# various parts of the MCMC sampler. It does not constitute any 
# part of the end product.

source("mcmc_functions.R")
source("sim_data.R")

## Define X and beta
nTrain <- 500
nTest <- round(0.01 * nTrain)
n <- nTrain + nTest
p <- 5
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
trueBeta <- runif(p, 0, 10)

# Define distribution bounds for uniform prior on theta
a <- 0.5
b <- 10

# Define true values for other parameters
trueSigma2 <- 2
trueTau2 <- 2
trueTheta <- 7

# Simulate Y_star and W_star (having nStar observations)
fullData <- simSpatialData(n = n, X = X, 
                           beta = trueBeta, 
                           theta = trueTheta, 
                           sigma2 = trueSigma2, 
                           tau2 = trueTau2)

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

# Dimensions for subset data and distance matrix for full/subset data
DTrain <- rdist(sTrain)
DTest <- rdist(sTest)
DCov <- rdist(sTrain, sTest)

# Generate phi and compress data
m <- round(0.02 * nTrain)
phi <- matrix(rnorm(m * nTrain, 0, sqrt(nTrain)), nrow = m, ncol = nTrain)
newY <- phi %*% yTrain
newX <- phi %*% xTrain

# Exploratory plots (log likelihood)
sigma2Vals <- seq(0.5, 20, by=0.1)
sigma2Density <- sapply(sigma2Vals, \(x) logLik(x, trueTau2, trueTheta, trueBeta))
plot(sigma2Vals, sigma2Density)
sigma2Vals[which.max(sigma2Density)]

tau2Vals <- seq(0.5, 20, by=0.1)
tau2Density <- sapply(tau2Vals, \(x) logLik(trueSigma2, x, trueTheta, trueBeta))
plot(tau2Vals, tau2Density)
tau2Vals[which.max(tau2Density)]

thetaVals <- seq(0.6, 9.9, by=0.1)
thetaDensity <- sapply(thetaVals, \(x) logLik(trueSigma2, trueTau2, x, trueBeta))
plot(thetaVals, thetaDensity)
thetaVals[which.max(thetaDensity)]

# Exploratory plots (log posterior)
sigma2Vals <- seq(0.5, 20, by=0.1)
sigma2Density <- sapply(sigma2Vals, \(x) logPost(log(x), log(trueTau2), f(trueTheta), trueBeta))
plot(sigma2Vals, sigma2Density)
sigma2Vals[which.max(sigma2Density)]

tau2Vals <- seq(0.5, 20, by=0.1)
tau2Density <- sapply(tau2Vals, \(x) logPost(log(trueSigma2), log(x), f(trueTheta), trueBeta))
plot(tau2Vals, tau2Density)
tau2Vals[which.max(tau2Density)]

thetaVals <- seq(0.6, 9.9, by=0.1)
thetaDensity <- sapply(thetaVals, \(x) logPost(log(trueSigma2), log(trueTau2), f(x), trueBeta))
plot(thetaVals, thetaDensity)
thetaVals[which.max(thetaDensity)]


