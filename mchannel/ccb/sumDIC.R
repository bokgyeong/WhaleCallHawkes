rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
# library(ggh4x) # facet_grid2
library(spgs) # chisq.unif.test
library(batchmeans); library(foreach)
library(xtable)
library(Rcpp); library(RcppArmadillo)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 2) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 2) }

fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

burn = 10000 # decided based on trace plots of -2logL


# =============================================================================-
# Compute DIC ----
# =============================================================================-

load(paste0('ccb/data/ccb.RData'))

ts = data$ts # unit is minutes
marks = data$marks
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
  Xm[[i]] = cbind(noise[,i+1],
                  sin(2*pi*(knts + 30)/(8*60)),
                  cos(2*pi*(knts + 30)/(8*60)),
                  sin(2*pi*(knts + 30)/(12*60)),
                  cos(2*pi*(knts + 30)/(12*60)),
                  sin(2*pi*(knts + 30)/(24*60)),
                  cos(2*pi*(knts + 30)/(24*60))) # should be greater than sback*4
}
p = ncol(Xm[[1]])
indlam0 = sapply(1:length(ts), function(i) which(knts >= ts[i])[1] - 1 - 1)


data.dic = c()

# -----------------------------------------------------------------------------=
## Model (1) NHPP ----
# -----------------------------------------------------------------------------=
fitj = 1
load(paste0('ccb/fit/ccbfit', fits[fitj], '.RData'))
load(paste0('ccb/loglik/ccbfit', fits[fitj], 'loglik.RData'))
print(paste0('fit ', fits[fitj], ': ', length(postLogLik)))

if((burn - (nrow(postBeta0) - length(postLogLik))) > 0){
  postLogLik = postLogLik[-(1:(burn - (nrow(postBeta0) - length(postLogLik))))] 
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp('src/RcppFtns.cpp')

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, maxT, lam0m, indlam0, knts, rep(0, K), 1, 1)

data.dic = rbind(data.dic,
                 data.frame(fit = fits[fitj], ExpDeviance = mean(-2*postLogLik),
                            lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
                            EffNumPar = mean(-2*postLogLik) - devianceAtMean,
                            DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))



# -----------------------------------------------------------------------------=
## Model (2) NHPP+GP ----
# -----------------------------------------------------------------------------=
fitj = 2
load(paste0('ccb/fit/ccbfit', fits[fitj], '.RData'))
load(paste0('ccb/loglik/ccbfit', fits[fitj], 'loglik.RData'))
print(paste0('fit ', fits[fitj], ': ', length(postLogLik)))

if((burn - (nrow(postBeta0) - length(postLogLik))) > 0){
  postLogLik = postLogLik[-(1:(burn - (nrow(postBeta0) - length(postLogLik))))] 
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postDelta = postDelta[-(1:burn),]

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp('src/RcppFtns.cpp')

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) + exp(mean(postDelta[,k])) * colMeans(postWm) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, maxT, lam0m, indlam0, knts, rep(0, K), 1, 1)

data.dic = rbind(data.dic,
                 data.frame(fit = fits[fitj], ExpDeviance = mean(-2*postLogLik),
                            lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
                            EffNumPar = mean(-2*postLogLik) - devianceAtMean,
                            DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))



# -----------------------------------------------------------------------------=
## Model (3) NHPP+SE ----
# -----------------------------------------------------------------------------=
fitj = 3
load(paste0('ccb/fit/ccbfit', fits[fitj], '.RData'))
load(paste0('ccb/loglik/ccbfit', fits[fitj], 'loglik.RData'))
print(paste0('data , fit ', fits[fitj], ': ', length(postLogLik)))

if((burn - (nrow(postBeta0) - length(postLogLik))) > 0){
  postLogLik = postLogLik[-(1:(burn - (nrow(postBeta0) - length(postLogLik))))] 
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

sourceCpp('src/RcppFtns.cpp')

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, maxT, lam0m, indlam0, knts, colMeans(postZeta), mean(postEta), mean(postPhi))

data.dic = rbind(data.dic,
                 data.frame(fit = fits[fitj], ExpDeviance = mean(-2*postLogLik),
                            lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
                            EffNumPar = mean(-2*postLogLik) - devianceAtMean,
                            DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))



# -----------------------------------------------------------------------------=
## Model (4) NHPP+GP+SE ----
# -----------------------------------------------------------------------------=
fitj = 4
load(paste0('ccb/fit/ccbfit', fits[fitj], '.RData'))
load(paste0('ccb/loglik/ccbfit', fits[fitj], 'loglik.RData'))
print(paste0('data , fit ', fits[fitj], ': ', length(postLogLik)))

if((burn - (nrow(postBeta0) - length(postLogLik))) > 0){
  postLogLik = postLogLik[-(1:(burn - (nrow(postBeta0) - length(postLogLik))))] 
}

postBeta0 = postBeta0[-(1:burn),]
for(i in 1:p){
  postBeta[[i]] = postBeta[[i]][-(1:burn),]
}
postWm = postWm[-(1:burn),]
postDelta = postDelta[-(1:burn),]
postZeta = postZeta[-(1:burn),]
postEta = postEta[-(1:burn)]
postPhi = postPhi[-(1:burn)]

postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][,k])
}

sourceCpp('src/RcppFtns.cpp')

lam0m = foreach(k = 1:K, .combine = cbind) %do% { exp( mean(postBeta0[,k]) + Xm[[k]] %*% colMeans(postBetanew[[k]]) + exp(mean(postDelta[,k])) * colMeans(postWm) ) }
devianceAtMean = -2 * compLogLiki(ts, marks, distmat, maxT, lam0m, indlam0, knts, colMeans(postZeta), mean(postEta), mean(postPhi))

data.dic = rbind(data.dic,
                 data.frame(fit = fits[fitj], ExpDeviance = mean(-2*postLogLik),
                            lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
                            EffNumPar = mean(-2*postLogLik) - devianceAtMean,
                            DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))


save(data.dic, file = 'ccb/fig/ccbDIC.RData')





# =============================================================================-
# Summary ----
# =============================================================================-
load('ccb/fig/ccbDIC.RData')

newlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+CC', '(iv) NHPP+GP+CC')

# data.dic = data.dic[-2,]

data.dic$fit = factor(data.dic$fit, levels = fits, labels = newlabs)

data.dic$ExpDeviance = format(round(data.dic$ExpDeviance), nsmall = 0)
data.dic$EffNumPar = format(round(data.dic$EffNumPar), nsmall = 0)
data.dic$DIC = format(round(data.dic$DIC), nsmall = 0)
data.dic$lb = format(round(data.dic$lb), nsmall = 0)
data.dic$ub = format(round(data.dic$ub), nsmall = 0)


data.dic = data.dic %>%
  mutate(HPD = paste0('(', lb, ', ', ub, ')')) %>% 
  select(fit, ExpDeviance, HPD, EffNumPar, DIC)


data.dic %>% 
  xtable() %>% 
  print(booktabs = F, include.rownames = F)

