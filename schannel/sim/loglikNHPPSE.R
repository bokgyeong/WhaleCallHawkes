rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)

load(paste0('sim/data/sim', datasets[run_ID], '.RData'))
load(paste0('sim/fit/sim', datasets[run_ID], 'NHPPSE.RData'))
filename = paste0('sim/loglik/sim', datasets[run_ID], 'NHPPSEloglik.RData')

p = length(beta)

betaInd = 1:p
alphaInd = p+1
etaInd = p+2

burn = 10000
postSamples = postSamples[-(1:burn),]
niters = nrow(postSamples)

c(etaInd, ncol(postSamples))



# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts
m = length(knts) - 1
p = ncol(Xm)

indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)

aaa = unique(c(0, seq(1000, niters, by = 1000), niters))

start = 1; postLogLik = c()
# load(filename); start = which(length(postLogLik) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = foreach(iter = (aaa[aa-1]+1):(aaa[aa]), .combine = c, .packages = 'Rcpp') %do% {
    
    lam0m = exp( Xm %*% postSamples[iter,betaInd] );
    compLogLiki(ts, maxT, lam0m, indlam0, knts, postSamples[iter,alphaInd], postSamples[iter,etaInd])
  }
  
  print(paste0('Computed ', aaa[aa], 'th iteration'))
  
  postLogLik = c(postLogLik, dummy)
  save(postLogLik, file = filename)
}
save(postLogLik, file = filename)




