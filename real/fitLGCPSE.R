rm(list = ls())
library(Rcpp); library(RcppArmadillo); library(foreach); library(tidyverse)

fold = 'real/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
ifelse(!dir.exists(path.fit), dir.create(path.fit, recursive = T), FALSE)

datai = 'ccb'
fiti = 'LGCPSE'

filename = paste0(path.fit, datai, '_', fiti, '.RData')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# ===========================================================================-
# Set up ----
# ===========================================================================-
load(paste0(path.data, datai, '.RData'))

ts = data$ts # unit is minutes
marks = data$marks
soundspeed = 1.5 * 60 # km per min

maxT = ceiling(max(ts))
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
rho_w = 60 # effective range is rho_w * 3



## covariates ----
noise = data.frame(ts = knts) %>%
  left_join(noise)
noise[,-1] = scale(noise[,-1])

K = nrow(distmat)
Xm = list()
for(i in 1:K){
  Xm[[i]] = cbind(
    noise[,i+1],
    sin(3 * 2 * pi * (knts + 15*60 + 1) / (24 * 60)), # 8h
    cos(3 * 2 * pi * (knts + 15*60 + 1) / (24 * 60)),
    sin(2 * 2 * pi * (knts + 15*60 + 1) / (24 * 60)), # 12h
    cos(2 * 2 * pi * (knts + 15*60 + 1) / (24 * 60)),
    sin(1 * 2 * pi * (knts + 15*60 + 1) / (24 * 60)), # 24h
    cos(1 * 2 * pi * (knts + 15*60 + 1) / (24 * 60))  # should be greater than sback*4 
  ) # should be greater than sback*4
}

p = ncol(Xm[[1]])
# par(mfrow = c(2, 2))
# for(i in 1:ncol(Xm[[1]])){
#   plot(knts, Xm[[1]][,i], type = 'l')
# }




## Initial values for parameters ----
beta0 = rnorm(K)
betaTilde0 = mean(beta0)
beta = matrix(rnorm(p*K), ncol = K)
betaTilde = rowMeans(beta)
tau = rgamma(1+p+1, shape = 10, scale = 0.1)
logdelta = rnorm(K)
logdeltaTilde = mean(logdelta)
zeta = rgamma(1, shape = 10, scale = 0.1)
eta = 1
phi = 0.2


## Compute Wm ----
tdiffm = matrix(0, m+1, m+1)
for(i in 1:(m+1)){
  for(j in 1:i){
    tdiffm[i,j] = tdiffm[j,i] = abs(knts[i] - knts[j])
  }
}
Sigmam = exp(-tdiffm / rho_w)
invSigmam = solve(Sigmam)
cholSigmam = chol(Sigmam)
Wm = t(cholSigmam) %*% rnorm(m+1)
indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)


## Others ----
sigma2 = rep(0.1, 1+nrow(beta)+1+1+1)
COVbeta0 = diag(length(beta0))
COVbeta = sapply(1:nrow(beta), function(i) diag(ncol(beta)), simplify = F)
COVlogdelta = diag(length(logdelta))
COVeta = 1
COVphi = 1


## parameters of invgamma priors
shape_tau = 2
rate_tau = 1


## parameters of gamma prior
shape_zeta = rate_zeta = 1/1000


## parameters of uniform prior
lb_eta = 3 / 20
ub_eta = 3 / 0.003948651
lb_phi = 3 / max(distmat)
ub_phi = 3 / min(distmat[upper.tri(distmat)])


niters = 100000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = rep(1, length(sigma2))
outers = unique(c(0, seq(1000, niters, by = 1000), niters))
# outers = c(0, niters)


start = 1; Accprob = 0; rtime = 0
postBeta0 = postTau = postLogdelta = postTilde = postZeta = postEta = postPhi = postWm = postBranching = c()
postBeta = sapply(1:p, function(i) c(), simplify = F)
# load(filename); start = which(outers == nrow(postBeta0))

# i = 2; outers = seq(0, 100, by = 50); adaptInterval = 25

sourceCpp(path.cpp)


for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = fitLGCPSE(
    outeri, ts, marks, Xm, maxT, soundspeed, rho_beta, rho_w, 
    distmat, knts, tdiffm, betaTilde0, beta0, betaTilde, beta, 
    tau, Wm, logdeltaTilde, logdelta, indlam0, zeta, eta, phi, shape_tau, 
    rate_tau, shape_zeta, rate_zeta, lb_eta, ub_eta, lb_phi, 
    ub_phi, sigma2, COVbeta0, COVbeta, COVlogdelta, COVeta, COVphi, 
    updateCOV, adaptInterval, adaptFactorExponent, adapIter, 
    1) 
  rtime = rtime + proc.time()[3] - ptm
  
  Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
  postBranching = rbind(postBranching, dummy$postBranching)
  postBeta0 = rbind(postBeta0, dummy$postBeta0)
  postWm = rbind(postWm, dummy$postWm)
  for(i in 1:p){
    postBeta[[i]] = rbind(postBeta[[i]], dummy$postBeta[[i]])
  }
  postTau = rbind(postTau, dummy$postTau)
  postLogdelta = rbind(postLogdelta, dummy$postLogdelta)
  postTilde = rbind(postTilde, dummy$postTilde)
  postZeta = rbind(postZeta, dummy$postZeta)
  postEta = c(postEta, dummy$postEta)
  postPhi = c(postPhi, dummy$postPhi)
  
  nSamples = nrow(postBeta0)
  beta0 = postBeta0[nSamples,]
  betaTilde0 = postTilde[nSamples,1]
  beta = t(sapply(1:p, function(i) postBeta[[i]][nSamples,]))
  betaTilde = postTilde[nSamples,2:(1+p)]
  Wm = dummy$Wm
  tau = postTau[nSamples,]
  logdelta = postLogdelta[nSamples,]
  logdeltaTilde = postTilde[nSamples,2+p]
  zeta = postZeta[nSamples,]
  eta = postEta[nSamples]
  phi = postPhi[nSamples]
  
  sigma2 = dummy$sigma2
  adapIter = dummy$adapIter
  COVbeta0 = dummy$COVbeta0
  COVbeta = dummy$COVbeta
  COVlogdelta = dummy$COVlogdelta
  COVeta = dummy$COVeta
  COVphi = dummy$COVphi
  
  save(sback, rtime, postBeta0, postBeta, postLogdelta, postTilde, postWm, 
       postTau, postZeta, postEta, postPhi, postBranching, Accprob,
       betaTilde, beta0, betaTilde, beta, tau, logdeltaTilde, logdelta, 
       zeta, eta, phi, Wm, rho_beta, rho_w,
       shape_tau, rate_tau, shape_zeta, rate_zeta,
       lb_eta, ub_eta, lb_phi, ub_phi,
       sigma2, adapIter, COVbeta0, COVbeta, COVlogdelta, COVeta, COVphi,
       updateCOV, adaptInterval, adaptFactorExponent, adapIter,
       Xm, rho_beta, rho_w, soundspeed,
       file = filename)
}


