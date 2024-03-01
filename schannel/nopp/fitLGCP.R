rm(list = ls())
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')

load('nopp/data/nopp.RData')
filename = paste0('nopp/fit/noppLGCP.RData')

# ===========================================================================-
# Set up ----
# ===========================================================================-
ts = data$ts
maxT = ceiling(max(ts))
sback = 20 # unit = min
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho = 60 # effective range is 3 hour


## covariates ----
noise = data.frame(ts = knts) %>%
  left_join(noise)
noiseVar = as.vector(scale(noise$noise))

noise$UTC[1]

Xm = cbind(1,
           noiseVar,
           sin(2*pi*(knts + 15*60 + 1)/(8*60)),
           cos(2*pi*(knts + 15*60 + 1)/(8*60)),
           sin(2*pi*(knts + 15*60 + 1)/(12*60)),
           cos(2*pi*(knts + 15*60 + 1)/(12*60)),
           sin(2*pi*(knts + 15*60 + 1)/(24*60)),
           cos(2*pi*(knts + 15*60 + 1)/(24*60))) # should be greater than sback*4 
p = ncol(Xm)
# par(mfrow = c(2,3))
# for(i in 2:p){
#   plot(knts, Xm[,i], type = 'l')
# }


## Initial values for parameters ----
beta = rnorm(p)
delta = log(1)


## Compute Wm ----
tdiffm = matrix(0, m+1, m+1)
for(i in 1:(m+1)){
  for(j in 1:i){
    tdiffm[i,j] = tdiffm[j,i] = abs(knts[i] - knts[j])
  }
}
Sigmam = exp(- tdiffm / rho)
invSigmam = solve(Sigmam)
cholSigmam = chol(Sigmam)
Wm = t(cholSigmam) %*% rnorm(m+1)
lam0m = as.vector( exp(Xm %*% beta + exp(delta) * Wm) )


## Approximation ----
n = length(ts)
indlam0 = sapply(1:n, function(i) which(knts >= ts[i])[1] - 1 - 1)
lam0 = compLam0(ts, maxT, lam0m, indlam0, knts)
intLam0 = compIntLam0(maxT, lam0m)
# par(mfrow = c(1, 1))
# plot(knts, lam0m, type = 'l')
# points(ts, lam0, col = 2)


## Others ---
betaInd = 1:p
deltaInd = p+1

sigma2 = rep(0.2^2, 2)
adapIter = rep(1, 2)
COVbeta = diag(p)
COVdelta = 1


niters = 100000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
outers = c(0, seq(1000, niters, by = 1000))

start = 1; postSamples = postWm = c(); Accprob = 0; rtime = 0
# load(filename); start = which(outers == nrow(postSamples))

# i = 2; outers = seq(0, 2000, by = 100); adaptInterval = 50

sourceCpp('src/RcppFtns.cpp')

for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = fitLGCP(
    outeri, ts, Xm, maxT, knts, tdiffm, beta, delta, rho, 
    Wm, indlam0, sigma2, COVbeta, COVdelta, updateCOV, adaptInterval, 
    adaptFactorExponent, adapIter) 
  rtime = rtime + proc.time()[3] - ptm
  
  postSamples = rbind(postSamples, dummy$postSamples)
  postWm = rbind(postWm, dummy$postWm)
  Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
  
  nSamples = nrow(postSamples)
  beta = postSamples[nSamples, betaInd]
  delta = postSamples[nSamples, deltaInd]
  
  sigma2 = dummy$sigma2
  adapIter = dummy$adapIter
  Wm = dummy$Wm
  COVbeta = dummy$COVbeta
  COVdelta = dummy$COVdelta
  
  save(rtime, postSamples, postWm, Accprob, beta, delta,
       sigma2, adapIter, Wm, COVbeta, COVdelta,
       updateCOV, adaptInterval, adaptFactorExponent, adapIter,
       sback, knts, rho, Xm,
       file = filename)
}


