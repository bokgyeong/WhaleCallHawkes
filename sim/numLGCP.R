rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

fold = 'sim/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.num = paste0(fold, 'num/')
ifelse(!dir.exists(path.num), dir.create(path.num, recursive = T), FALSE)


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 1
datai = paste0('sim', datasets[runID])
fiti = 'LGCP'

filename = paste0(path.num, datai, '_', fiti, '_num.RData')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# ===========================================================================-
# Set up ----
# ===========================================================================-
load(paste0(path.data, datai, '.RData'))
load(paste0(path.fit, datai, '_', fiti, '.RData'))

p = nrow(beta)
K = nrow(distmat)

# burn = 10000
burn = 100

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
# rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3
# par(mfrow = c(2, 2))
# for(i in 1:ncol(Xm[[1]])){
#   plot(knts, Xm[[1]][,i], type = 'l')
# }

indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)

# aaa = unique(c(0, seq(1000, niters, by = 1000), niters))
aaa = c(0,niters)

start = 1; postNum = c()
# load(filename); start = which(nrow(postNum) == aaa)


sourceCpp(path.cpp)

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(iter in (aaa[aa-1]+1):(aaa[aa])){
    
    numBack = sapply(1:K, function(k) compIntLam0k(maxT, postBeta0[iter,k], Xm[[k]], postBetanew[[k]][iter,], exp(postLogdelta[iter,k]), postWm[iter,]))
    
    dummy = rbind(dummy, numBack)
  }
  
  print(paste0('Computed ', aaa[aa], 'th iteration'))
  
  postNum = rbind(postNum, dummy)
  save(postNum, file = filename)
}
save(postNum, file = filename)




