rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

load('nopp/data/nopp.RData')
load(paste0('nopp/fit/noppLGCP.RData'))
filename = paste0('nopp/rtct/noppLGCPrtct.RData')


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
maxT = ceiling(max(ts))
m = length(knts) - 1
p = ncol(Xm)


aaa = unique(c(0, seq(10, length(ts), by = 10), length(ts)))


start = 1; postCompen = c(); intlami = 0
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



