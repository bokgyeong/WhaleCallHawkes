rm(list = ls())
library(Rcpp); library(RcppArmadillo); library(foreach); library(readr); library(tidyverse)

fold = 'sim/'
path.data = paste0(fold, 'data/')
ifelse(!dir.exists(path.data), dir.create(path.data, recursive = T), FALSE)

datai = 'simLGCPSE'
filename = paste0(path.data, datai, '.RData')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# =============================================================================-
# True parameters from the model fit ----
# =============================================================================-

## space and time window ----
load(paste0(path.data, 'ccb_env.RData'))

maxT = 5 * 24 * 60 # unit = min
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs in km
rho_w = 60 # effective range is rho_w * 3

soundspeed = 1.5 * 60 # km per min

## covariates ----
noise = data.frame(ts = knts) %>%
  left_join(noise)
noise[,-1] = scale(noise[,-1])

K = nrow(distmat)
Xm = list()
for(i in 1:K){
  Xm[[i]] = cbind(noise[,i+1],
                  sin(2*pi*(knts + 30)/(8*60)),
                  cos(2*pi*(knts + 30)/(8*60)),
                  sin(2*pi*(knts + 30)/(12*60)),
                  cos(2*pi*(knts + 30)/(12*60)),
                  sin(2*pi*(knts + 30)/(24*60)),
                  cos(2*pi*(knts + 30)/(24*60))) # should be greater than sback*4
}
p = ncol(Xm[[1]])
# par(mfrow = c(2, 2))
# for(i in 2:ncol(Xm[[1]])){
#   plot(knts, Xm[[1]][,i], type = 'l')
# }


## true parameter values ----
load(paste0(path.data, 'ccbfitLGCPSE.RData'))

# burn = 10000
# tau = round(colMeans(postTau[-(1:burn),]), 2)
# beta0 = round(colMeans(postBeta0[-(1:burn),]), 2) + 0.8
# beta = t(sapply(1:length(postBeta), function(i) round(colMeans(postBeta[[i]][-(1:burn),]), 2)))
# delta = exp(round(colMeans(postDelta[-(1:burn),]), 2))
# zeta = round(colMeans(postZeta[-(1:burn),]), 2)
# eta = round(mean(postEta[-(1:burn)]), 3)
# phi = round(mean(postPhi[-(1:burn)]), 3)


trpar = list(tau = tau, beta0 = beta0, beta = beta, delta = delta, 
             zeta = zeta, eta = eta, phi = phi, rho_beta = rho_beta, rho_w = rho_w)




# =============================================================================-
# Simulate data
# =============================================================================-

# load(filename)

source(path.r)
sourceCpp(path.cpp)


# Compute Wm ----
tdiffm = matrix(0, m+1, m+1)
for(i in 1:(m+1)){
  for(j in 1:i){
    tdiffm[i,j] = tdiffm[j,i] = abs(knts[i] - knts[j])
  }
}
Sigmam = exp(-tdiffm / rho_w)
invSigmam = solve(Sigmam)
cholSigmam = chol(Sigmam)
Wm = t(cholSigmam) %*% rnorm(m+1)
trpar$Wm = Wm
lam0m = sapply(1:K, function(i) as.vector( exp(trpar$beta0[i] + Xm[[i]] %*% trpar$beta[,i] + trpar$delta[i] * trpar$Wm) ))
aaa = foreach(k = 1:K, .combine = c) %do% {
  compIntLam0k(maxT, beta0[k], Xm[[k]], beta[,k], delta[k], Wm)
}

compIntLam0(maxT, beta0, Xm, beta, delta, Wm) 
sum(zeta/eta * aaa)
compIntLam0(maxT, beta0, Xm, beta, delta, Wm) + sum(zeta/eta * aaa)


# Simulate data ----
data = simulateMultiselfex(
  Xm, distmat, maxT, lam0m, knts, trpar$beta0, trpar$beta, 
  trpar$Wm, trpar$delta, trpar$zeta, trpar$eta, trpar$phi,
  soundspeed)

ts = data$ts[,1]
marks = data$ts[,2] - 1
branching = data$branching
n = length(ts)
sapply(1:K, function(k) sum(branching[marks == (k-1)] == 0))
sapply(1:K, function(k) sum(marks[branching[branching > 0]] == (k-1)))
sapply(1:K, function(k) sum(marks[branching[branching > 0]] == (k-1)) / sum(marks == (k-1)))



# Approximation ----
indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)
indmarks = sapply(1:K, function(i) which(marks == (i-1))-1)
lam0 = compLam0(ts, marks, maxT, trpar$beta0, Xm, trpar$beta, trpar$delta, trpar$Wm, indlam0, knts)
intLam0 = compIntLam0(maxT, trpar$beta0, Xm, trpar$beta, trpar$delta, trpar$Wm) 
# par(mfrow = c(2,3))
# for(i in 1:K){
#   plot(knts, lam0m[,i], type = 'l')
#   points(ts[indmarks[[i]] + 1], lam0[indmarks[[i]] + 1], col = 2)
# }
# sum(data$branching == 0)
# compIntLam0(maxT, trpar$beta0, Xm, trpar$beta, trpar$delta, trpar$Wm) 
# compIntLam0_parallel(maxT, trpar$beta0, Xm, trpar$beta, trpar$delta, trpar$Wm, 2) 



# Save data ----
save(sback, maxT, knts, m, rho_beta, rho_w, Xm, p, distmat, trpar, lam0m, data, indlam0, soundspeed,
     file = filename)

