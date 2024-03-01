rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/fit/ccbfitLGCPSE.RData'))
filename = paste0('ccb/lam/ccbfitLGCPSElam.RData')

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

tsnew = seq(min(ts), max(ts), length.out = 2000)

aaa = unique(c(0, seq(100, length(tsnew), by = 100), length(tsnew)))

start = 1; postBack = postSE = postLam = c()
# load(filename); start = which(nrow(postLam) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummyBack = dummySE = dummyLam = c()
  for(i in (aaa[aa-1]+1):(aaa[aa])) {
    
    indlam0i = which(knts >= tsnew[i])[1] - 1 - 1
    
    lami = compLamt(tsnew[i], ts, marks, distmat, knts, maxT, indlam0i, Xm, postBeta0, 
                    postBetanew, exp(postDelta), postWm, postZeta, postEta, postPhi)
    
    dummyBack = rbind(dummyBack, as.vector(lami$resBack))
    dummySE = rbind(dummySE, as.vector(lami$resSE))
    dummyLam = rbind(dummyLam, as.vector(lami$resLam))
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  postBack = rbind(postBack, dummyBack)
  postSE = rbind(postSE, dummySE)
  postLam = rbind(postLam, dummyLam)
  
  save(tsnew, postBack, postSE, postLam, file = filename)
}
save(tsnew, postBack, postSE, postLam, file = filename)



