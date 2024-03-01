rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)


load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/fit/ccbfitNHPPSE.RData'))
filename = paste0('ccb/rtct/ccbfitNHPPSErtct.RData')

p = nrow(beta)
K = nrow(distmat)

burn = 10000

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
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


aaa = unique(c(0, seq(100, length(ts), by = 100), length(ts)))

start = 1; postCompen = c()
# load(filename); start = which(nrow(postCompen) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(i in (aaa[aa-1]+1):(aaa[aa])) {
    
    indlam0i = which(knts >= ts[i])[1] - 1 - 1
    
    if(i == 1){
      intlampre = matrix(0, nrow = niters, ncol = K)
    } else {
      intlampre = intlami
    }
    
    intlam0i = rtctIntLam0ik(ts[i], indlam0i, maxT, knts, Xm, postBeta0, postBetanew, matrix(0, nrow = niters, ncol = K), matrix(0, nrow = niters, ncol = nrow(Xm[[1]])))
    intTrigi = rtctAlphaSumIntHik(ts[i], ts, marks, distmat, postZeta, postEta, postPhi)
    intlami = intlam0i + intTrigi
    
    dummy = rbind(dummy,
                  c(quantile(rowSums(intlami), probs = c(0.025, 0.5, 0.975)), 
                    quantile(rowSums(intlami) - rowSums(intlampre), probs = c(0.025, 0.5, 0.975)),
                    sapply(1:K, function(k) quantile(intlami[,k], probs = c(0.025, 0.5, 0.975))),
                    sapply(1:K, function(k) quantile(intlami[,k] - intlampre[,k], probs = c(0.025, 0.5, 0.975)))))
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  postCompen = rbind(postCompen, dummy)
  save(intlami, ts, postCompen, file = filename)
}
save(intlami, ts, postCompen, file = filename)


