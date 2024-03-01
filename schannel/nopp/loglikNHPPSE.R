rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

load('nopp/data/nopp.RData')
load(paste0('nopp/fit/noppNHPPSE.RData'))
filename = paste0('nopp/loglik/noppNHPPSEloglik.RData')

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
maxT = ceiling(max(ts))
m = length(knts) - 1
p = ncol(Xm)


indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)


aaa = unique(c(0, seq(1000, niters, by = 1000), niters))

start = 1; postLogLik = c()
# load(filename); start = which(length(postLogLik) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(iter in (aaa[aa-1]+1):(aaa[aa])){
    
    lam0m = exp( Xm %*% postSamples[iter,betaInd] );
    dummy = c(dummy, compLogLiki(ts, maxT, lam0m, indlam0, knts, postSamples[iter,alphaInd], postSamples[iter,etaInd]))
  }
  
  print(paste0('Computed ', aaa[aa], 'th iteration'))
  
  postLogLik = c(postLogLik, dummy)
  save(postLogLik, file = filename)
}
save(postLogLik, file = filename)



