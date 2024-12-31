rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)


fold = 'real/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.rtct = paste0(fold, 'rtct/')

ifelse(!dir.exists(path.rtct), dir.create(path.rtct, recursive = T), FALSE)


datai = 'ccb'
fiti = 'NHPPSE'

filename = paste0(path.rtct, datai, '_', fiti, '_rtct.RData')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# =============================================================================-
# Load results ----
# =============================================================================-
load(paste0(path.data, datai, '.RData'))
load(paste0(path.fit, datai, '_', fiti, '.RData'))

p = nrow(beta)
K = nrow(distmat)

burn = 75000

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
soundspeed = 1.5 * 60 # km per min

maxT = ceiling(max(ts))
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3


aaa = unique(c(0, seq(100, length(ts), by = 100), length(ts)))

start = 1; postCompen = foreach(k = 1:(K+1)) %do% {c()}
# load(filename); start = which(ncol(postCompen[[1]]) == aaa)

sourceCpp(path.cpp)

for(aa in (start+1):length(aaa) ){
  
  for(i in (aaa[aa-1]+1):(aaa[aa])) {
    
    indlam0i = which(knts >= ts[i])[1] - 1 - 1
    
    if(i == 1){
      intlampre = matrix(0, nrow = niters, ncol = K)
    } else {
      intlampre = intlami
    }
    
    intlam0i = rtctIntLam0ik(ts[i], indlam0i, maxT, knts, Xm, postBeta0, postBetanew, matrix(0, nrow = niters, ncol = K), matrix(0, nrow = niters, ncol = nrow(Xm[[1]])))
    intTrigi = rtctAlphaSumIntHik(ts[i], ts, marks, distmat, soundspeed, postZeta, postEta, postPhi)
    intlami = intlam0i + intTrigi
    
    dstarSum = rowSums(intlami) - rowSums(intlampre)
    postCompen[[1]] = cbind(postCompen[[1]], dstarSum)
    
    dstar = intlami - intlampre
    for(ll in 1:K){
      postCompen[[ll+1]] = cbind(postCompen[[ll+1]], dstar[,ll])
    }
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  save(intlami, ts, postCompen, file = filename)
}
save(intlami, ts, postCompen, file = filename)


