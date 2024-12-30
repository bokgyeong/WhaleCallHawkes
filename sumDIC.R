rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
library(Rcpp); library(RcppArmadillo)
hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 2) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 2) }


path.data = paste0('data/')
path.fit = paste0('fit/')
path.loglik = paste0('loglik/')
path.fig = paste0('fig/')
path.sum = paste0('sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


# datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
datasets = c('NHPP')
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


path.r = paste0('src/RFtns.R')
path.cpp = paste0('src/RcppFtns.cpp')
# =============================================================================-
# Compute DID ----
# =============================================================================-
# burn = 10000 # decided based on trace plots of -2logL
burn = 100

data.dic = c()

for(runID in 1:length(datasets)){

  # ---------------------------------------------------------------------------=
  ## dataset ----
  # ---------------------------------------------------------------------------=
  datai = paste0('sim', datasets[runID])
  load(paste0(path.data, datai, '.RData'))

  ts = data$ts[,1] # unit is minutes
  marks = data$ts[,2]-1
  soundspeed = 1.5 * 60 # km per min
  sback = 20
  knts = unique(c(0, seq(0, maxT, by = sback), maxT))
  m = length(knts) - 1
  # rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
  # rho_w = 60 # effective range is rho_w * 3
  p = ncol(Xm[[1]])
  K = nrow(distmat)

  indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)

  # -----------------------------------------------------------------------------=
  ## Model (1) NHPP ----
  # -----------------------------------------------------------------------------=
  fiti = 'NHPP'
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

  postBetanew = foreach(k = 1:K) %do% {
    sapply(1:p, function(j) postBeta[[j]][,k])
  }

  sourceCpp(path.cpp)

  lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) ) }
  devianceAtMean = -2 * compLogLiki(ts, marks, distmat, soundspeed, maxT, lam0m, indlam0, knts, rep(0, K), 1, 1)

  data.dic = rbind(
    data.dic,
    data.frame(
      data = datai, fit = fiti, ExpDeviance = mean(-2*postLogLik),
      lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
      EffNumPar = mean(-2*postLogLik) - devianceAtMean,
      DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
    )
  )


  # -----------------------------------------------------------------------------=
  ## Model (2) NHPP+GP ----
  # -----------------------------------------------------------------------------=
  fiti = 'LGCP'
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
      data = datai, fit = fiti, ExpDeviance = mean(-2*postLogLik),
      lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
      EffNumPar = mean(-2*postLogLik) - devianceAtMean,
      DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
    )
  )



  # -----------------------------------------------------------------------------=
  ## Model (3) NHPP+SE ----
  # -----------------------------------------------------------------------------=
  fiti = 'NHPPSE'
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
      data = datai, fit = fiti, ExpDeviance = mean(-2*postLogLik),
      lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
      EffNumPar = mean(-2*postLogLik) - devianceAtMean,
      DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
    )
  )



  # -----------------------------------------------------------------------------=
  ## Model (4) NHPP+GP+SE ----
  # -----------------------------------------------------------------------------=
  fiti = 'LGCPSE'
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
      data = datai, fit = fiti, ExpDeviance = mean(-2*postLogLik),
      lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
      EffNumPar = mean(-2*postLogLik) - devianceAtMean,
      DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)
    )
  )

  save(data.dic, file = paste0(path.sum, 'simDIC.RData'))
}




# =============================================================================-
# Table ----
# =============================================================================-
load(paste0(path.sum, 'simDIC.RData'))

# datlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')
datlabs = c('(i) NHPP')
fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')

data.dic$data = factor(data.dic$data, levels = unique(data.dic$data), labels = datlabs)
data.dic$fit = factor(data.dic$fit, levels = unique(data.dic$fit), labels = fitlabs)

data.dic$ExpDeviance = format(round(data.dic$ExpDeviance), nsmall = 0)
data.dic$EffNumPar = format(round(data.dic$EffNumPar), nsmall = 0)
data.dic$DIC = format(round(data.dic$DIC), nsmall = 0)
data.dic$lb = format(round(data.dic$lb), nsmall = 0)
data.dic$ub = format(round(data.dic$ub), nsmall = 0)


data.dic %>%
  mutate(HPD = paste0('(', lb, ', ', ub, ')')) %>%
  select(data, fit, ExpDeviance, HPD, EffNumPar, DIC) %>%
  xtable() %>%
  print(booktabs = T, include.rownames = F)


data.dic %>% 
  select(data, fit, DIC) %>% 
  pivot_wider(names_from = fit, values_from = DIC) %>%
  xtable() %>%
  print(booktabs = T, include.rownames = F)


data.dic %>% 
  select(data, fit, DIC) %>% 
  mutate(DIC = format(round(as.numeric(DIC)/10000, 2), nsmall = 2)) %>% 
  pivot_wider(names_from = fit, values_from = DIC) %>%
  xtable() %>%
  print(booktabs = T, include.rownames = F)






