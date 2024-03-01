rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/fit/ccbfitLGCPSE.RData'))
filename = paste0('ccb/loglik/ccbfitLGCPSEloglik.RData')

p = nrow(beta)
K = nrow(distmat)

burn = 10000

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postDelta = postDelta[-(1:burn),]
postZeta = postZeta[-(1:burn),]
postEta = postEta[-(1:burn)]
postPhi = postPhi[-(1:burn)]

niters = nrow(postBeta0)

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}



# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts # unit is minutes
marks = data$marks
maxT = ceiling(max(ts))
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
rho_w = 60 # effective range is rho_w * 3

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

indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)

aaa = unique(c(0, seq(200, niters, by = 200), niters))

start = 1; postLogLik = c()
# load(filename); start = which(length(postLogLik) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = foreach(iter = (aaa[aa-1]+1):(aaa[aa]), .combine = c, .packages = 'Rcpp') %do% {
    lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( postBeta0[iter,k] + Xm[[k]] %*% postBetanew[[k]][iter,] + exp(postDelta[iter,k]) * postWm[iter,] ) }
    compLogLiki(ts, marks, distmat, maxT, lam0m, indlam0, knts, postZeta[iter,], postEta[iter], postPhi[iter])
  }
  
  print(paste0('Computed ', aaa[aa], 'th iteration'))
  
  postLogLik = c(postLogLik, dummy)
  save(postLogLik, file = filename)
}
save(postLogLik, file = filename)




