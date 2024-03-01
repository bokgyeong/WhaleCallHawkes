rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

load('nopp/data/nopp.RData')
load(paste0('nopp/fit/noppLGCPSE.RData'))
filename = paste0('nopp/rtct/noppLGCPSErtct.RData')

p = length(beta)

betaInd = 1:p
deltaInd = p+1
alphaInd = p+2
etaInd = p+3


burn = 10000
postSamples = postSamples[-(1:burn),]
postWm = postWm[-(1:burn),]

niters = nrow(postSamples)



# =============================================================================-
# Posterior intensity ----
# =============================================================================-
ts = data$ts
maxT = ceiling(max(ts))
m = length(knts) - 1
p = ncol(Xm)


aaa = unique(c(0, seq(10, length(ts), by = 10), length(ts)))


start = 1; postCompen = c()
# load(filename); start = which(nrow(postCompen) == aaa)

sourceCpp('src/RcppFtns.cpp')

for(aa in (start+1):length(aaa) ){
  
  dummy = c()
  for(i in (aaa[aa-1]+1):(aaa[aa])){
    
    indlam0i = which(knts >= ts[i])[1] - 1 - 1
    
    if(i == 1){
      intlampre = 0
    } else {
      intlampre = intlami
    }
    
    intlam0i = rtctIntLam0i(ts[i], Xm, maxT, knts, indlam0i, matrix(postSamples[,betaInd], nrow = niters), exp(postSamples[,deltaInd]) * postWm)
    intTrigi = rtctSumIntHi(ts[i], ts, postSamples[,alphaInd], postSamples[,etaInd])
    intlami = intlam0i + intTrigi
    di = intlami - intlampre
    
    dummy = rbind(dummy, c(quantile(intlami, probs = c(0.025, 0.5, 0.975)), quantile(di, probs = c(0.025, 0.5, 0.975))))
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  postCompen = rbind(postCompen, dummy)
  save(intlami, ts, postCompen, file = filename)
}
save(intlami, ts, postCompen, file = filename)


