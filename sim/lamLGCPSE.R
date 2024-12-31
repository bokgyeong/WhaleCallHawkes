rm(list = ls())
library(foreach)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

fold = 'sim/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.lam = paste0(fold, 'lam/')
ifelse(!dir.exists(path.lam), dir.create(path.lam, recursive = T), FALSE)


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 1
datai = paste0('sim', datasets[runID])
fiti = 'LGCPSE'

filename = paste0(path.lam, datai, '_', fiti, '_lam.RData')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# =============================================================================-
# Load results ----
# =============================================================================-
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
ts = data$ts[,1] # unit is minutes
marks = data$ts[,2]-1
soundspeed = 1.5 * 60 # km per min
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
# rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3

tsnew = seq(min(ts), max(ts), length.out = 2000)

aaa = unique(c(0, seq(100, length(tsnew), by = 100), length(tsnew)))

start = 1; postBack = postSE = postLam = c()
# load(filename); start = which(nrow(postLam) == aaa)

sourceCpp(path.cpp)

for(aa in (start+1):length(aaa) ){
  
  dummyBack = dummySE = dummyLam = c()
  for(i in (aaa[aa-1]+1):(aaa[aa])) {
    
    indlam0i = which(knts >= tsnew[i])[1] - 1 - 1
    
    lami = compLamt(tsnew[i], ts, marks, distmat, soundspeed, knts, maxT, indlam0i, Xm, postBeta0, 
                    postBetanew, exp(postLogdelta), postWm, postZeta, postEta, postPhi)
    
    dummyBack = rbind(dummyBack, as.vector(lami$resBack))
    dummySE = rbind(dummySE, as.vector(lami$resSE))
    dummyLam = rbind(dummyLam, as.vector(lami$resLam))
  }
  
  print(paste0('Completed by ', aaa[aa], 'th time event'))
  
  postBack = rbind(postBack, dummyBack)
  postSE = rbind(postSE, dummySE)
  postLam = rbind(postLam, dummyLam)
  
  save(tsnew, postBack, postSE, postLam, file = filename)
}
save(tsnew, postBack, postSE, postLam, file = filename)



