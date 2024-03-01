rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)

load(paste0('sim/data/sim', datasets[run_ID], '.RData'))
load(paste0('sim/fit/sim', datasets[run_ID], 'NHPP.RData'))
filename = paste0('sim/num/sim', datasets[run_ID], 'NHPPnum.RData')

p = length(beta)

betaInd = 1:p

burn = 10000
postSamples = postSamples[-(1:burn),]

niters = nrow(postSamples)


# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts
sback = 20 # unit = min
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho = 60 # effective range is 3 hour
p = ncol(Xm)


indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)

aaa = unique(c(0, seq(1000, niters, by = 1000), niters))

start = 1; postNum = c()
# load(filename); start = which(nrow(postNum) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(iter in (aaa[aa-1]+1):(aaa[aa])){
    
    lam0m = exp( Xm %*% postSamples[iter,betaInd] );
    numBack = compIntLam0(maxT, lam0m)
    numSE = 0
    
    dummy = rbind(dummy, c(numBack, numSE))
  }
  
  print(paste0('Computed ', aaa[aa], 'th iteration'))
  
  postNum = rbind(postNum, dummy)
  save(postNum, file = filename)
}
save(postNum, file = filename)

