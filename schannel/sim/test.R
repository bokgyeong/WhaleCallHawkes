rm(list = ls())
library(Rcpp); library(RcppArmadillo)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

# run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_ID = 3
set.seed(run_ID)



load(paste0('sim/data/sim', datasets[run_ID], '.RData'))
filename = paste0('sim/test/sim', datasets[run_ID], 'NHPPSE.RData')

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
alpha = 1
eta = 1


lam0m = as.vector( exp(Xm %*% beta ) )


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
alphaInd = p+1
etaInd = p+2

sigma2 = rep(0.2^2, 2)
adapIter = rep(1, 2)
COVbeta = diag(p)
COVeta = 1

lb_eta = 3 / 20
ub_eta = 3 / min(diff(ts))
shape_alpha = rate_alpha = 0.001

niters = 20000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
outers = c(0, seq(1000, niters, by = 1000))

start = 1; postSamples = postBranching = c(); Accprob = 0; rtime = 0
# load(filename); start = which(outers == nrow(postSamples))

# i = 2; outers = seq(0, 2000, by = 100); adaptInterval = 50

sourceCpp('src/RcppFtns.cpp')

for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = fitNHPPSE(
    outeri, ts, Xm, maxT, knts, beta, alpha, eta, indlam0, 
    shape_alpha, rate_alpha, lb_eta, ub_eta, sigma2, COVbeta, 
    COVeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter) 
  rtime = rtime + proc.time()[3] - ptm
  
  postSamples = rbind(postSamples, dummy$postSamples)
  postBranching = rbind(postBranching, dummy$postBranching)
  Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
  
  nSamples = nrow(postSamples)
  beta = postSamples[nSamples, betaInd]
  alpha = postSamples[nSamples, alphaInd]
  eta = postSamples[nSamples, etaInd]
  
  sigma2 = dummy$sigma2
  adapIter = dummy$adapIter
  COVbeta = dummy$COVbeta
  COVeta = dummy$COVeta
  
  save(rtime, postSamples, postBranching, Accprob, beta, alpha, eta, 
       lb_eta, ub_eta, shape_alpha, rate_alpha,
       sigma2, adapIter, COVbeta, COVeta,
       updateCOV, adaptInterval, adaptFactorExponent, adapIter,
       sback, knts, rho, Xm,
       file = filename)
}


