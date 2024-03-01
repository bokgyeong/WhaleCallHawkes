rm(list = ls())
library(Rcpp); library(RcppArmadillo)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)


load(paste0('sim/data/sim', datasets[run_ID], '.RData'))
filename = paste0('sim/fit/sim', datasets[run_ID], 'LGCPSE.RData')

# ===========================================================================-
# Set up ----
# ===========================================================================-
ts = data$ts
sback = 20 # unit = min
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho = 60 # effective range is 3 hour
p = ncol(Xm)

# par(mfrow = c(2,3))
# for(i in 2:p){
#   plot(knts, Xm[,i], type = 'l')  
# }


## Initial values for parameters ----
beta = rnorm(p)
delta = log(1)
alpha = 1
eta = 1


## Compute Wm ----
tdiffm = matrix(0, m+1, m+1)
for(i in 1:(m+1)){
  for(j in 1:i){
    tdiffm[i,j] = tdiffm[j,i] = abs(knts[i] - knts[j])
  }
}
Sigmam = exp(- tdiffm / rho)
invSigmam = solve(Sigmam)
cholSigmam = chol(Sigmam)
Wm = t(cholSigmam) %*% rnorm(m+1)
lam0m = as.vector( exp(Xm %*% beta + exp(delta) * Wm) )


## Approximation ----
n = length(ts)
indlam0 = sapply(1:n, function(i) which(knts >= ts[i])[1] - 1 - 1)
lam0 = compLam0(ts, maxT, lam0m, indlam0, knts)
intLam0 = compIntLam0(maxT, lam0m)
# par(mfrow = c(1, 1))
# plot(knts, lam0m, type = 'l')
# points(ts, lam0, col = 2)


## Others ---
betaInd = 1:p
deltaInd = p+1
alphaInd = p+2
etaInd = p+3

sigma2 = rep(0.2^2, 3)
adapIter = rep(1, 3)
COVbeta = diag(p)
COVdelta = COVeta = 1

shape_kappa = 2
rate_kappa = 1
lb_eta = 3 / 20
ub_eta = 3 / min(diff(ts))
shape_alpha = rate_alpha = 0.001

niters = 100000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
outers = c(0, seq(1000, niters, by = 1000))

start = 1; postSamples = postBranching = postWm = c(); Accprob = 0; rtime = 0
# load(filename); start = which(outers == nrow(postSamples))

# i = 2; outers = seq(0, 2000, by = 100); adaptInterval = 50

sourceCpp('src/RcppFtns.cpp')

for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = fitLGCPSE(
    outeri, ts, Xm, maxT, knts, tdiffm, beta, delta, rho, 
    alpha, eta, Wm, indlam0, shape_alpha, rate_alpha, lb_eta, 
    ub_eta, sigma2, COVbeta, COVdelta, COVeta, updateCOV, adaptInterval, 
    adaptFactorExponent, adapIter) 
  rtime = rtime + proc.time()[3] - ptm
  
  postSamples = rbind(postSamples, dummy$postSamples)
  postBranching = rbind(postBranching, dummy$postBranching)
  postWm = rbind(postWm, dummy$postWm)
  Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
  
  nSamples = nrow(postSamples)
  beta = postSamples[nSamples, betaInd]
  delta = postSamples[nSamples, deltaInd]
  alpha = postSamples[nSamples, alphaInd]
  eta = postSamples[nSamples, etaInd]
  
  sigma2 = dummy$sigma2
  adapIter = dummy$adapIter
  Wm = dummy$Wm
  COVbeta = dummy$COVbeta
  COVdelta = dummy$COVdelta
  COVeta = dummy$COVeta
  
  save(rtime, postSamples, postBranching, postWm, Accprob, beta, delta, alpha, eta, 
       lb_eta, ub_eta, shape_alpha, rate_alpha,
       sigma2, adapIter, Wm, COVbeta, COVdelta, COVeta,
       updateCOV, adaptInterval, adaptFactorExponent, adapIter,
       file = filename)
}



# c(trpar$beta, trpar$kappa, trpar$eta, trpar$tau, trpar$alpha)
# 
# # load(filename)
# dim(postSamples)
# burn = 1
# # burn = 1000
# colMeans(postSamples[-(1:burn),])
# 
# ts.plot(postSamples[-(1:burn),betaInd[1]]); abline(h = trpar$beta[1], col = 2)
# ts.plot(postSamples[-(1:burn),betaInd[2]]); abline(h = trpar$beta[2], col = 2)
# ts.plot(postSamples[-(1:burn),betaInd[3]]); abline(h = trpar$beta[3], col = 2)
# ts.plot(postSamples[-(1:burn),betaInd[4]]); abline(h = trpar$beta[4], col = 2)
# ts.plot(postSamples[-(1:burn),betaInd[5]]); abline(h = trpar$beta[5], col = 2)
# ts.plot(postSamples[-(1:burn),betaInd[6]]); abline(h = trpar$beta[6], col = 2)
# ts.plot(postSamples[-(1:burn),kappaInd]); abline(h = trpar$kappa, col = 2)
# # ts.plot(postSamples[-(1:burn),rhoInd]); abline(h = trpar$rho, col = 2)
# ts.plot(postSamples[-(1:burn),alphaInd]); abline(h = trpar$alpha, col = 2)
# ts.plot(postSamples[-(1:burn),etaInd]); abline(h = trpar$eta, col = 2)
# 
# Accprob


