rm(list = ls())
library(Rcpp); library(RcppArmadillo); library(foreach); library(tidyverse)

fold = 'sim/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
ifelse(!dir.exists(path.fit), dir.create(path.fit, recursive = T), FALSE)


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 1
datai = paste0('sim', datasets[runID])
fiti = 'NHPP'

filename = paste0(path.fit, datai, '_', fiti, '.RData')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# ===========================================================================-
# Set up ----
# ===========================================================================-
load(paste0(path.data, datai, '.RData'))


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



## Initial values for parameters ----
beta0 = rnorm(K)
betaTilde0 = mean(beta0)
beta = matrix(rnorm(p*K), ncol = K)
betaTilde = rowMeans(beta)
tau = rgamma(1+p, shape = 10, scale = 0.1)
indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)



## Others ----
sigma2 = rep(0.1, 1+nrow(beta))
COVbeta0 = diag(length(beta0))
COVbeta = sapply(1:nrow(beta), function(i) diag(ncol(beta)), simplify = F)



## parameters of invgamma priors
shape_tau = 2
rate_tau = 1


niters = 100000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = rep(1, length(sigma2))
outers = unique(c(0, seq(1000, niters, by = 1000), niters))
# outers = c(0, niters)


start = 1; postBeta0 = postTilde = postTau = c(); postBeta = sapply(1:p, function(i) c(), simplify = F); Accprob = 0; rtime = 0
# load(filename); start = which(outers == nrow(postBeta0))

# i = 2; outers = seq(0, 100, by = 50); adaptInterval = 25

sourceCpp(path.cpp)


for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = fitNHPP(
    outeri, ts, marks, Xm, maxT, rho_beta, distmat, knts, 
    betaTilde0, beta0, betaTilde, beta, tau, indlam0, shape_tau, 
    rate_tau, sigma2, COVbeta0, COVbeta, updateCOV, adaptInterval, 
    adaptFactorExponent, adapIter, 1) 
  rtime = rtime + proc.time()[3] - ptm
  
  Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
  postBeta0 = rbind(postBeta0, dummy$postBeta0)
  for(i in 1:p){
    postBeta[[i]] = rbind(postBeta[[i]], dummy$postBeta[[i]])
  }
  postTilde = rbind(postTilde, dummy$postTilde)
  postTau = rbind(postTau, dummy$postTau)
  
  nSamples = nrow(postBeta0)
  beta0 = postBeta0[nSamples,]
  betaTilde0 = postTilde[nSamples,1]
  beta = t(sapply(1:p, function(i) postBeta[[i]][nSamples,]))
  betaTilde = postTilde[nSamples,2:(1+p)]
  tau = postTau[nSamples,]
  
  sigma2 = dummy$sigma2
  adapIter = dummy$adapIter
  COVbeta0 = dummy$COVbeta0
  COVbeta = dummy$COVbeta
  
  save(sback, rtime, postBeta0, postBeta, postTilde, postTau,
       Accprob, betaTilde0, beta0, betaTilde, beta, tau, rho_beta,
       shape_tau, rate_tau, sigma2, adapIter, COVbeta0, COVbeta, 
       updateCOV, adaptInterval, adaptFactorExponent, adapIter,
       Xm, rho_beta, rho_w,
       file = filename)
}



