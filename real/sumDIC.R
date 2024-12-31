rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
library(Rcpp); library(RcppArmadillo)
hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 2) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 2) }

fold = 'real/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.loglik = paste0(fold, 'loglik/')
path.fig = paste0(fold, 'fig/')
path.sum = paste0(fold, 'sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


datai = 'ccb'
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')

# =============================================================================-
# Compute DIC ----
# =============================================================================-
load(paste0(path.data, datai, '.RData'))

ts = data$ts # unit is minutes
marks = data$marks
soundspeed = 1.5 * 60 # km per min
maxT = ceiling(max(ts))
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3

indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)


# burn = 75000
burn = 50

data.dic = c()

# -----------------------------------------------------------------------------=
## Model (1) NHPP ----
# -----------------------------------------------------------------------------=
fiti = fits[1]
load(paste0(path.fit, datai, '_', fiti, '.RData'))
load(paste0(path.loglik, datai, '_', fiti, '_loglik.RData'))

print(paste0('data ', datai, ', fit ', fiti, ': ', length(postLogLik)))

p = nrow(beta)
K = nrow(distmat)

addburn = burn - (nrow(postBeta0) - length(postLogLik))
if(addburn > 0){
  postLogLik = postLogLik[-(1:addburn)]  
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp(path.cpp)

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, soundspeed, maxT, lam0m, indlam0, knts, rep(0, K), 1, 1)

data.dic = rbind(
  data.dic,
  data.frame(
    fit = fiti, ExpDeviance = mean(-2*postLogLik),
    lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
    EffNumPar = mean(-2*postLogLik) - devianceAtMean,
    DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
  )
)




# -----------------------------------------------------------------------------=
## Model (2) NHPP+GP ----
# -----------------------------------------------------------------------------=
fiti = fits[2]
load(paste0(path.fit, datai, '_', fiti, '.RData'))
load(paste0(path.loglik, datai, '_', fiti, '_loglik.RData'))

print(paste0('data ', datai, ', fit ', fiti, ': ', length(postLogLik)))

addburn = burn - (nrow(postBeta0) - length(postLogLik))
if(addburn > 0){
  postLogLik = postLogLik[-(1:addburn)]  
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postLogdelta = postLogdelta[-(1:burn),]

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp(path.cpp)

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) + exp(mean(postLogdelta[,k])) * colMeans(postWm) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, soundspeed, maxT, lam0m, indlam0, knts, rep(0, K), 1, 1)

data.dic = rbind(
  data.dic,
  data.frame(
    fit = fiti, ExpDeviance = mean(-2*postLogLik),
    lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
    EffNumPar = mean(-2*postLogLik) - devianceAtMean,
    DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
  )
)



# -----------------------------------------------------------------------------=
## Model (3) NHPP+SE ----
# -----------------------------------------------------------------------------=
fiti = fits[3]
load(paste0(path.fit, datai, '_', fiti, '.RData'))
load(paste0(path.loglik, datai, '_', fiti, '_loglik.RData'))

print(paste0('data ', datai, ', fit ', fiti, ': ', length(postLogLik)))

addburn = burn - (nrow(postBeta0) - length(postLogLik))
if(addburn > 0){
  postLogLik = postLogLik[-(1:addburn)]  
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postZeta = postZeta[-(1:burn),]
postEta = postEta[-(1:burn)]
postPhi = postPhi[-(1:burn)]

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp(path.cpp)

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, soundspeed, maxT, lam0m, indlam0, knts, colMeans(postZeta), mean(postEta), mean(postPhi))

data.dic = rbind(
  data.dic,
  data.frame(
    fit = fiti, ExpDeviance = mean(-2*postLogLik),
    lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
    EffNumPar = mean(-2*postLogLik) - devianceAtMean,
    DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
  )
)



# -----------------------------------------------------------------------------=
## Model (4) NHPP+GP+SE ----
# -----------------------------------------------------------------------------=
fiti = fits[4]
load(paste0(path.fit, datai, '_', fiti, '.RData'))
load(paste0(path.loglik, datai, '_', fiti, '_loglik.RData'))

print(paste0('data ', datai, ', fit ', fiti, ': ', length(postLogLik)))

addburn = burn - (nrow(postBeta0) - length(postLogLik))
if(addburn > 0){
  postLogLik = postLogLik[-(1:addburn)]  
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postLogdelta = postLogdelta[-(1:burn),]
postZeta = postZeta[-(1:burn),]
postEta = postEta[-(1:burn)]
postPhi = postPhi[-(1:burn)]

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp(path.cpp)

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) + exp(mean(postLogdelta[,k])) * colMeans(postWm) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, soundspeed, maxT, lam0m, indlam0, knts, colMeans(postZeta), mean(postEta), mean(postPhi))

data.dic = rbind(
  data.dic,
  data.frame(
    fit = fiti, ExpDeviance = mean(-2*postLogLik),
    lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
    EffNumPar = mean(-2*postLogLik) - devianceAtMean,
    DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
  )
)


save(data.dic, file =  paste0(path.sum, datai, 'DIC.RData'))





# =============================================================================-
# Table ----
# =============================================================================-
load(paste0(path.sum, datai, 'DIC.RData'))

newlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')

data.dic$fit = factor(data.dic$fit, levels = unique(data.dic$fit), labels = newlabs)

tab.dic = data.dic %>% 
  mutate(
    ExpDeviance = format(round(ExpDeviance), nsmall = 0),
    EffNumPar = format(round(EffNumPar), nsmall = 0),
    DIC = format(round(DIC), nsmall = 0),
    lb = format(round(lb), nsmall = 0),
    ub = format(round(ub), nsmall = 0)
  ) %>% 
  mutate(HPD = paste0('(', lb, ', ', ub, ')')) %>% 
  select(fit, ExpDeviance, HPD, EffNumPar, DIC)


tab.dic %>% 
  select(-c(EffNumPar)) %>% 
  xtable() %>% 
  print(booktabs = F, include.rownames = F)

