rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(runID)

load(paste0('sim/data/sim', datasets[runID], '.RData'))
load(paste0('sim/fit/sim', datasets[runID], 'fitLGCP.RData'))
filename = paste0('sim/rtct/sim', datasets[runID], 'fitLGCPrtct.RData')

p = nrow(beta)
K = nrow(distmat)

burn = 10000

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postDelta = postDelta[-(1:burn),]

niters = nrow(postBeta0)

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}



# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts[,1] # unit is minutes
marks = data$ts[,2]-1
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
rho_w = 60 # effective range is rho_w * 3
p = ncol(Xm[[1]])
K = nrow(distmat)
# par(mfrow = c(2, 2))
# for(i in 1:ncol(Xm[[1]])){
#   plot(knts, Xm[[1]][,i], type = 'l')
# }

aaa = unique(c(0, seq(20, length(ts), by = 20), length(ts)))

start = 1; postCompen = c()
# load(filename); start = which(nrow(postCompen) == aaa)


sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(i in (aaa[aa-1]+1):(aaa[aa])){
    
    indlam0i = which(knts >= ts[i])[1] - 1 - 1
    
    if(i == 1){
      intlampre = matrix(0, nrow = niters, ncol = K)
    } else {
      intlampre = intlami
    }
    
    intlam0i = rtctIntLam0ik(ts[i], indlam0i, maxT, knts, Xm, postBeta0, postBetanew, exp(postDelta), postWm)
    intlami = intlam0i
    
    dummy = rbind(dummy, c(quantile(rowSums(intlami), probs = c(0.025, 0.5, 0.975)), 
                           quantile(rowSums(intlami) - rowSums(intlampre), probs = c(0.025, 0.5, 0.975)),
                           sapply(1:K, function(k) quantile(intlami[,k], probs = c(0.025, 0.5, 0.975))),
                           sapply(1:K, function(k) quantile(intlami[,k] - intlampre[,k], probs = c(0.025, 0.5, 0.975)))))
    
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  postCompen = rbind(postCompen, dummy)
  save(intlami, ts, postCompen, file = filename)
}
save(intlami, ts, postCompen, file = filename)


