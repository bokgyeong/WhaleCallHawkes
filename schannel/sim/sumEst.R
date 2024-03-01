rm(list = ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# sum = function(x){ paste0(format(round(bm(x)$est, 2), nsmall = 2), ' (', 
#                           format(round(HPDinterval(as.mcmc(x))[1], 2), nsmall = 2), ', ',
#                           format(round(HPDinterval(as.mcmc(x))[2], 2), nsmall = 2), ')') }
bmmean = function(x) { format(round(bm(x)$est, 1), nsmall = 1) }
hpd = function(x){ paste0('(', 
                          format(round(HPDinterval(as.mcmc(x))[1], 1), nsmall = 1), ', ',
                          format(round(HPDinterval(as.mcmc(x))[2], 1), nsmall = 1), ')') }
hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 1) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 1) }


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


datai = datasets[1]
load(paste0('simUnif2/data/sim', datai, '.RData'))

p = ncol(Xm)
burn = 10000

row.names = c('Intercept', 'Noise',
              paste0('8h ', c('sine', 'cosine')), 
              paste0('12h ', c('sine', 'cosine')), 
              paste0('24h ', c('sine', 'cosine')), 'Variance')
row.names.SE = c('Intercept', 'Noise',
                 paste0('8h ', c('sine', 'cosine')), 
                 paste0('12h ', c('sine', 'cosine')), 
                 paste0('24h ', c('sine', 'cosine')), 'Variance', 'SE intensity', 'SE Decay')

par = c(paste0('$\\beta_{', 0:(p-1), '}$'), '$\\kappa$')
par.SE = c(paste0('$\\beta_{', 0:(p-1), '}$'), '$\\kappa$', '$\\alpha^{*}$', '$\\eta$')



# -----------------------------------------------------------------------------=
# NHPP ----
# -----------------------------------------------------------------------------=
load(paste0('simUnif2/fit/sim', datai, 'NHPP.RData'))

betaInd = 1:p
posterior = postSamples[-(1:burn), c(betaInd)]

sumNHPP = data.frame(
  par = par,
  truth = format(c(trpar$beta, trpar$kappa), nsmall = 1),
  dummytruth = c(trpar$beta, trpar$kappa),
  mean = c(apply(posterior, 2, bmmean), NA),
  dummyhpd = c(apply(posterior, 2, hpd), NA),
  hpd1 = c(apply(posterior, 2, hpd1), NA),
  hpd2 = c(apply(posterior, 2, hpd2), NA),
  # row.names = c(paste0('beta', 0:(p-1)), 'alpha', 'eta')) %>% 
  row.names = row.names) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  # select(c(par, truth, mean, hpd, sig))
  select(c(par, truth, mean, hpd))




# -----------------------------------------------------------------------------=
# LGCP ----
# -----------------------------------------------------------------------------=
load(paste0('simUnif2/fit/sim', datai, 'LGCP.RData'))

betaInd = 1:p
kappaInd = p+1
phiInd = p+2

posterior = postSamples[-(1:burn), c(betaInd, kappaInd)]


sumLGCP = data.frame(
  dummytruth = c(trpar$beta, ifelse(is.na(trpar$kappa), 0, trpar$kappa)),
  mean = c(apply(posterior, 2, bmmean)),
  dummyhpd = c(apply(posterior, 2, hpd)),
  hpd1 = c(apply(posterior, 2, hpd1)),
  hpd2 = c(apply(posterior, 2, hpd2))) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  # select(c(mean, hpd, sig))
  select(c(mean, hpd))





# -----------------------------------------------------------------------------=
# NHPP + SE ----
# -----------------------------------------------------------------------------=
load(paste0('simUnif2/fit/sim', datai, 'NHPPSE.RData'))

betaInd = 1:p
etaInd = p+1
alphaInd = p+2

postSamples[,alphaInd] = postSamples[,alphaInd] / postSamples[,etaInd]
posterior = postSamples[-c(1:burn), c(betaInd, alphaInd, etaInd)]

sumNHPPSE = data.frame(
  par = par.SE,
  truth = format(round(c(trpar$beta, trpar$kappa, trpar$alpha/trpar$eta, trpar$eta), 1), nsmall = 1),
  dummytruth = c(trpar$beta, ifelse(is.na(trpar$kappa), 0, trpar$kappa), ifelse(is.na(trpar$alpha/trpar$eta), 0, trpar$alpha/trpar$eta), ifelse(is.na(trpar$eta), 0, trpar$eta)),
  mean = c(apply(posterior, 2, bmmean)[1:p], NA, apply(posterior, 2, bmmean)[(p+1):(p+2)]),
  dummyhpd = c(apply(posterior, 2, hpd)[1:p], NA, apply(posterior, 2, hpd)[(p+1):(p+2)]),
  hpd1 = c(apply(posterior, 2, hpd1)[1:p], NA, apply(posterior, 2, hpd1)[(p+1):(p+2)]),
  hpd2 = c(apply(posterior, 2, hpd2)[1:p], NA, apply(posterior, 2, hpd2)[(p+1):(p+2)]),
  row.names = row.names.SE) %>%
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  # select(c(par, truth, mean, hpd, sig))
  select(c(par, truth, mean, hpd))



# -----------------------------------------------------------------------------=
# LGCP + SE ----
# -----------------------------------------------------------------------------=
load(paste0('simUnif2/fit/sim', datai, 'LGCPSE.RData'))

betaInd = 1:p
kappaInd = p+1
phiInd = p+2
etaInd = p+3
alphaInd = p+4

postSamples[,alphaInd] = postSamples[,alphaInd] / postSamples[,etaInd]
posterior = postSamples[-c(1:burn), c(betaInd, kappaInd, alphaInd, etaInd)]


sumLGCPSE = data.frame(
  dummytruth = c(trpar$beta, ifelse(is.na(trpar$kappa), 0, trpar$kappa), ifelse(is.na(trpar$alpha/trpar$eta), 0, trpar$alpha/trpar$eta), ifelse(is.na(trpar$eta), 0, trpar$eta)),
  mean = c(apply(posterior, 2, bmmean)),
  dummyhpd = c(apply(posterior, 2, hpd)),
  hpd1 = c(apply(posterior, 2, hpd1)),
  hpd2 = c(apply(posterior, 2, hpd2))) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  # select(c(mean, hpd, sig))
  select(c(mean, hpd))



# -----------------------------------------------------------------------------=
# Table
# -----------------------------------------------------------------------------=
cbind(sumNHPP, sumLGCP) %>% 
  xtable() %>% 
  print(booktabs = T, sanitize.text.function = function(x) {x})

cbind(sumNHPPSE, sumLGCPSE) %>% 
  xtable() %>% 
  print(booktabs = T, sanitize.text.function = function(x) {x})
  



