rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

path.data = paste0('data/')
path.fit = paste0('fit/')
path.rtct = paste0('rtct/')
ifelse(!dir.exists(path.rtct), dir.create(path.rtct, recursive = T), FALSE)


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 1
datai = paste0('sim', datasets[runID])
fiti = 'LGCP'

filename = paste0(path.rtct, datai, '_', fiti, '_rtct.RData')


path.r = paste0('src/RFtns.R')
path.cpp = paste0('src/RcppFtns.cpp')
# ===========================================================================-
# Set up ----
# ===========================================================================-
load(paste0(path.data, datai, '.RData'))
load(paste0(path.fit, datai, '_', fiti, '.RData'))

p = nrow(beta)
K = nrow(distmat)

burn = 10000

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postLogdelta = postLogdelta[-(1:burn),]

niters = nrow(postBeta0)

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}



# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts[,1] # unit is minutes
marks = data$ts[,2]-1
soundspeed = 1.5 * 60 # km per min
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3

# par(mfrow = c(2, 2))
# for(i in 1:ncol(Xm[[1]])){
#   plot(knts, Xm[[1]][,i], type = 'l')
# }

aaa = unique(c(0, seq(100, length(ts), by = 100), length(ts)))


start = 1; postCompen = foreach(k = 1:(K+1)) %do% {c()}
# load(filename); start = which(ncol(postCompen[[1]]) == aaa)


sourceCpp(path.cpp)

for(aa in (start+1):length(aaa) ){
  
  for(i in (aaa[aa-1]+1):(aaa[aa])){
    
    indlam0i = which(knts >= ts[i])[1] - 1 - 1
    
    if(i == 1){
      intlampre = matrix(0, nrow = niters, ncol = K)
    } else {
      intlampre = intlami
    }
    
    intlam0i = rtctIntLam0ik(ts[i], indlam0i, maxT, knts, Xm, postBeta0, postBetanew, exp(postLogdelta), postWm)
    intlami = intlam0i
    
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


