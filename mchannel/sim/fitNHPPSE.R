rm(list = ls())
library(Rcpp); library(RcppArmadillo); library(foreach)
# ncores = 20
ncores = 1


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(runID)

load(paste0('sim/data/sim', datasets[runID], '.RData'))
filename = paste0('sim/fit/sim', datasets[runID], 'fitNHPPSE.RData')

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')
# sourceCpp('src/test.cpp')
# ===========================================================================-
# Set up ----
# ===========================================================================-
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
zeta = rgamma(1, shape = 10, scale = 0.1)
eta = 1
phi = 0.2
indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)



## Others ----
sigma2 = rep(0.1, 1+nrow(beta)+1+1)
COVbeta0 = diag(length(beta0))
COVbeta = sapply(1:nrow(beta), function(i) diag(ncol(beta)), simplify = F)
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
outers = c(0, seq(1000, niters, by = 1000))

start = 1; postBeta0 = postTilde = postTau = postZeta = postEta = postPhi = postBranching = c(); postBeta = sapply(1:p, function(i) c(), simplify = F); Accprob = 0; rtime = 0
# load(filename); start = which(outers == nrow(postBeta0))

# i = 2; outers = seq(0, 100, by = 50); adaptInterval = 25
# sourceCpp('src/test.cpp')

sourceCpp('src/RcppFtns.cpp')


for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = fitNHPPSE(
    outeri, ts, marks, Xm, maxT, rho_beta, distmat, knts, 
    betaTilde0, beta0, betaTilde, beta, tau, indlam0, zeta, eta, 
    phi, shape_tau, rate_tau, shape_zeta, rate_zeta, lb_eta, 
    ub_eta, lb_phi, ub_phi, sigma2, COVbeta0, COVbeta, COVeta, 
    COVphi, updateCOV, adaptInterval, adaptFactorExponent, adapIter, 
    ncores)  
  rtime = rtime + proc.time()[3] - ptm
  
  Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
  postBranching = rbind(postBranching, dummy$postBranching)
  postBeta0 = rbind(postBeta0, dummy$postBeta0)
  for(i in 1:p){
    postBeta[[i]] = rbind(postBeta[[i]], dummy$postBeta[[i]])
  }
  postTilde = rbind(postTilde, dummy$postTilde)
  postTau = rbind(postTau, dummy$postTau)
  postZeta = rbind(postZeta, dummy$postZeta)
  postEta = c(postEta, dummy$postEta)
  postPhi = c(postPhi, dummy$postPhi)
  
  nSamples = nrow(postBeta0)
  beta0 = postBeta0[nSamples,]
  betaTilde0 = postTilde[nSamples,1]
  beta = t(sapply(1:p, function(i) postBeta[[i]][nSamples,]))
  betaTilde = postTilde[nSamples,2:(1+p)]
  tau = postTau[nSamples,]
  zeta = postZeta[nSamples,]
  eta = postEta[nSamples]
  phi = postPhi[nSamples]
  
  sigma2 = dummy$sigma2
  adapIter = dummy$adapIter
  COVbeta0 = dummy$COVbeta0
  COVbeta = dummy$COVbeta
  COVeta = dummy$COVeta
  COVphi = dummy$COVphi
  
  save(sback, rtime, postBeta0, postBeta, postTilde, postTau,
       postZeta, postEta, postPhi, postBranching, Accprob,
       betaTilde0, beta0, betaTilde, beta, tau, zeta, eta, phi, rho_beta,
       shape_tau, rate_tau, shape_zeta, rate_zeta,
       lb_eta, ub_eta, lb_phi, ub_phi,
       sigma2, adapIter, COVbeta0, COVbeta, COVeta, COVphi,
       updateCOV, adaptInterval, adaptFactorExponent, adapIter,
       file = filename)
}




