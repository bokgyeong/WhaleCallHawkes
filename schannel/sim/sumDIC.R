rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
# library(ggh4x) # facet_grid2
library(spgs) # chisq..test
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

datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

burn = 10000 # decided based on trace plots of -2logL


# =============================================================================-
# Compute DID ----
# =============================================================================-

# data.dic = c()
# 
# for(i in 1:length(datasets)){
#   
#   # ---------------------------------------------------------------------------=
#   ## dataset ----
#   # ---------------------------------------------------------------------------=
#   load(paste0('sim/data/sim', datasets[i], '.RData'))  
#   ts = data$ts
#   # maxT = ceiling(max(ts))
#   sback = 20 # unit = min
#   knts = unique(c(0, seq(0, maxT, by = sback), maxT))
#   m = length(knts) - 1
#   phi = 60 # effective range is 3 hour
#   p = ncol(Xm)
#   
#   
#   # -----------------------------------------------------------------------------=
#   ## Model (1) NHPP ----
#   # -----------------------------------------------------------------------------=
#   j = 1
#   load(paste0('sim/fit/sim', datasets[i], fits[j], '.RData'))
#   load(paste0('sim/loglik/sim', datasets[i], fits[j], 'loglik.RData'))
#   print(paste0('data ', datasets[i], ', fit ', fits[j], ': ', length(postLogLik)))
#   
#   if((burn - (nrow(postSamples) - length(postLogLik))) > 0){
#     postLogLik = postLogLik[-(1:(burn - (nrow(postSamples) - length(postLogLik))))]
#     postSamples = postSamples[-(1:burn),]  
#   }
#   
#   
#   p = length(beta)
#   betaInd = 1:p
#   
#   sourceCpp('src/RcppFtns.cpp')
#   
#   lam0m = exp( Xm %*% colMeans(postSamples[,betaInd]) );
#   devianceAtMean = -2 * compLogLiki(ts, maxT, lam0m, indlam0, knts, 0, 1)
#   
#   data.dic = rbind(data.dic,
#                    data.frame(data = datasets[i], fit = fits[j], ExpDeviance = mean(-2*postLogLik),
#                               lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
#                               EffNumPar = mean(-2*postLogLik) - devianceAtMean,
#                               DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))
#   
#   
#   
#   # -----------------------------------------------------------------------------=
#   ## Model (2) NHPP+GP ----
#   # -----------------------------------------------------------------------------=
#   j = 2
#   load(paste0('sim/fit/sim', datasets[i], fits[j], '.RData'))
#   load(paste0('sim/loglik/sim', datasets[i], fits[j], 'loglik.RData'))
#   print(paste0('data ', datasets[i], ', fit ', fits[j], ': ', length(postLogLik)))
#   
#   if((burn - (nrow(postSamples) - length(postLogLik))) > 0){
#     postLogLik = postLogLik[-(1:(burn - (nrow(postSamples) - length(postLogLik))))]
#     postSamples = postSamples[-(1:burn),]
#     postWm = postWm[-(1:burn),]
#   }
#   
#   p = length(beta)
#   
#   betaInd = 1:p
#   deltaInd = p+1
#   
#   sourceCpp('src/RcppFtns.cpp')
#   
#   lam0m = exp( Xm %*% colMeans(postSamples[,betaInd]) + exp(mean(postSamples[,deltaInd])) * colMeans(postWm) );
#   devianceAtMean = -2 * compLogLiki(ts, maxT, lam0m, indlam0, knts, 0, 1)
#   
#   data.dic = rbind(data.dic,
#                    data.frame(data = datasets[i], fit = fits[j], ExpDeviance = mean(-2*postLogLik),
#                               lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
#                               EffNumPar = mean(-2*postLogLik) - devianceAtMean,
#                               DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))
#   
#   
#   
#   # -----------------------------------------------------------------------------=
#   ## Model (3) NHPP+SE ----
#   # -----------------------------------------------------------------------------=
#   j = 3
#   load(paste0('sim/fit/sim', datasets[i], fits[j], '.RData'))
#   load(paste0('sim/loglik/sim', datasets[i], fits[j], 'loglik.RData'))
#   print(paste0('data ', datasets[i], ', fit ', fits[j], ': ', length(postLogLik)))
# 
#   if((burn - (nrow(postSamples) - length(postLogLik))) > 0){
#     postLogLik = postLogLik[-(1:(burn - (nrow(postSamples) - length(postLogLik))))]
#     postSamples = postSamples[-(1:burn),]
#   }
# 
#   p = length(beta)
# 
#   betaInd = 1:p
#   alphaInd = p+1
#   etaInd = p+2
# 
#   sourceCpp('src/RcppFtns.cpp')
# 
#   lam0m = exp( Xm %*% colMeans(postSamples[,betaInd]) );
#   devianceAtMean = -2 * compLogLiki(ts, maxT, lam0m, indlam0, knts, mean(postSamples[,alphaInd]), mean(postSamples[,etaInd]))
# 
#   data.dic = rbind(data.dic,
#                    data.frame(data = datasets[i], fit = fits[j], ExpDeviance = mean(-2*postLogLik),
#                               lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
#                               EffNumPar = mean(-2*postLogLik) - devianceAtMean,
#                               DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))
# 
# 
# 
#   # -----------------------------------------------------------------------------=
#   ## Model (4) NHPP+GP+SE ----
#   # -----------------------------------------------------------------------------=
#   j = 4
#   load(paste0('sim/fit/sim', datasets[i], fits[j], '.RData'))
#   load(paste0('sim/loglik/sim', datasets[i], fits[j], 'loglik.RData'))
#   print(paste0('data ', datasets[i], ', fit ', fits[j], ': ', length(postLogLik)))
# 
#   if((burn - (nrow(postSamples) - length(postLogLik))) > 0){
#     postLogLik = postLogLik[-(1:(burn - (nrow(postSamples) - length(postLogLik))))]
#     postSamples = postSamples[-(1:burn),]
#     postWm = postWm[-(1:burn),]
#   }
# 
#   p = length(beta)
# 
#   betaInd = 1:p
#   deltaInd = p+1
#   alphaInd = p+2
#   etaInd = p+3
# 
#   sourceCpp('src/RcppFtns.cpp')
# 
#   lam0m = exp( Xm %*% colMeans(postSamples[,betaInd]) + exp(mean(postSamples[,deltaInd])) * colMeans(postWm) );
#   devianceAtMean = -2 * compLogLiki(ts, maxT, lam0m, indlam0, knts, mean(postSamples[,alphaInd]), mean(postSamples[,etaInd]))
# 
#   data.dic = rbind(data.dic,
#                    data.frame(data = datasets[i], fit = fits[j], ExpDeviance = mean(-2*postLogLik),
#                               lb = hpd1(-2*postLogLik), ub = hpd2(-2*postLogLik),
#                               EffNumPar = mean(-2*postLogLik) - devianceAtMean,
#                               DIC = mean(-2*postLogLik) - devianceAtMean + mean(-2*postLogLik)))
#   
#   
#   save(data.dic, file = 'sim/fig/simDIC.RData')
# }










# =============================================================================-
# Summary ----
# =============================================================================-
load('sim/fig/simDIC.RData')

datlabs = fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+SE', '(iv) NHPP+GP+SE')
# fitlabs = c('(i) NHPP', '(ii) NHPP+GP')
# fits = fits[1:2]

data.dic$data = factor(data.dic$data, levels = datasets, labels = datlabs)
data.dic$fit = factor(data.dic$fit, levels = fits, labels = fitlabs)

data.dic$ExpDeviance = format(round(data.dic$ExpDeviance), nsmall = 0)
data.dic$EffNumPar = format(round(data.dic$EffNumPar), nsmall = 0)
data.dic$DIC = format(round(data.dic$DIC), nsmall = 0)
data.dic$lb = format(round(data.dic$lb), nsmall = 0)
data.dic$ub = format(round(data.dic$ub), nsmall = 0)


data.dic = data.dic %>%
  mutate(HPD = paste0('(', lb, ', ', ub, ')')) %>%
  select(data, fit, ExpDeviance, HPD, EffNumPar, DIC)


data.dic %>%
  xtable() %>%
  print(booktabs = F, include.rownames = F)


data.dic %>% 
  select(data, fit, DIC) %>% 
  pivot_wider(names_from = fit, values_from = DIC) %>%
  xtable() %>%
  print(booktabs = T, include.rownames = F)




