rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)

load(paste0('sim/data/sim', datasets[run_ID], '.RData'))
load(paste0('sim/fit/sim', datasets[run_ID], 'LGCP.RData'))
filename = paste0('sim/rtct/sim', datasets[run_ID], 'LGCPrtct.RData')

p = length(beta)

betaInd = 1:p
deltaInd = p+1

burn = 10000
postSamples = postSamples[-(1:burn),]
postWm = postWm[-(1:burn),]

niters = nrow(postSamples)

c(deltaInd, ncol(postSamples))


# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts
sback = 20 # unit = min
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho = 60 # effective range is 3 hour
p = ncol(Xm)

aaa = unique(c(0, seq(100, length(ts), by = 100), length(ts)))

start = 1; postCompen = c()
# load(filename); start = which(nrow(postCompen) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(i in (aaa[aa-1]+1):(aaa[aa]) ) {
    
    indlam0i = which(knts >= ts[i])[1] - 1 - 1
    
    if(i == 1){
      intlampre = 0
    } else {
      intlampre = intlami
    }
    
    intlam0i = rtctIntLam0i(ts[i], Xm, maxT, knts, indlam0i, postSamples[,betaInd], exp(postSamples[,deltaInd]) * postWm)
    intlami = intlam0i
    di = intlami - intlampre
    
    dummy = rbind(dummy, c(quantile(intlami, probs = c(0.025, 0.5, 0.975)), quantile(di, probs = c(0.025, 0.5, 0.975))))
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  postCompen = rbind(postCompen, dummy)
  save(intlami, ts, postCompen, file = filename)
}
save(intlami, ts, postCompen, file = filename)



