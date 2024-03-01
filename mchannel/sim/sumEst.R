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

bmmean = function(x) { format(round(bm(x)$est, 1), nsmall = 1) }
hpd = function(x){ paste0('(', 
                          format(round(HPDinterval(as.mcmc(x))[1], 1), nsmall = 1), ', ',
                          format(round(HPDinterval(as.mcmc(x))[2], 1), nsmall = 1), ')') }
hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 1) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 1) }


datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


datai = datasets[2]
load(paste0('sim2/data/sim', datai, '.RData'))

if(datai %in% datasets[c(1, 3)]){
  trpar$tau = c(trpar$tau[1:(p+1)], 0)
}

p = ncol(Xm[[1]])
k = 10


burn = 10000

cov.names = c('Intercept', 'Noise',
              paste0('8h ', c('sine', 'cosine')), 
              paste0('12h ', c('sine', 'cosine')), 
              paste0('24h ', c('sine', 'cosine')))

row.names = c(cov.names, paste0(cov.names, ' variance'),
              'GP coefficient variance', 'GP coefficient')
row.names.se = c(row.names, 'SE intensity', 'Decay in time', 'Decay in space')

par = c(paste0('$\\beta_{', 0:p, ',', k, '}$'), paste0('$\\tau_{', 0:p, '}$'), '\\tau_{\\delta}', paste0('\\delta_{', k, '}'))
par.se = c(par, paste0('$\\alpha_{', k, '}$'), '$\\eta$', '$\\phi$')



# -----------------------------------------------------------------------------=
# (i) NHPP ----
# -----------------------------------------------------------------------------=
load(paste0('sim2/fit/sim', datai, 'fit', fits[1], '.RData'))


posterior = cbind(postBeta0[-(1:burn),k], 
                  foreach(i = 1:p, .combine = 'cbind') %do% { postBeta[[i]][-(1:burn),k] },
                  postTau[-(1:burn),])


sumNHPP = data.frame(
  row.names = row.names,
  par = par,
  truth = format(round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k]), 1), nsmall = 1),
  dummytruth = round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k]), 1),
  mean = c(apply(posterior, 2, bmmean), NA, NA),
  dummyhpd = c(apply(posterior, 2, hpd), NA, NA),
  hpd1 = c(apply(posterior, 2, hpd1), NA, NA),
  hpd2 = c(apply(posterior, 2, hpd2), NA, NA)) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  select(c(par, truth, mean, hpd, sig))




# -----------------------------------------------------------------------------=
# (ii) NHPP + GP ----
# -----------------------------------------------------------------------------=
load(paste0('sim2/fit/sim', datai, 'fit', fits[2], '.RData'))


posterior = cbind(postBeta0[-(1:burn),k], 
                  foreach(i = 1:p, .combine = 'cbind') %do% { postBeta[[i]][-(1:burn),k] },
                  postTau[-(1:burn),], postDelta[-(1:burn),k])


sumNHPPGP = data.frame(
  par = par,
  truth = format(round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k]), 1), nsmall = 1),
  dummytruth = round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k]), 1),
  mean = c(apply(posterior, 2, bmmean)),
  dummyhpd = c(apply(posterior, 2, hpd)),
  hpd1 = c(apply(posterior, 2, hpd1)),
  hpd2 = c(apply(posterior, 2, hpd2))) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  select(c(mean, hpd, sig))




# -----------------------------------------------------------------------------=
# (iii) NHPP + SE ----
# -----------------------------------------------------------------------------=
load(paste0('sim2/fit/sim', datai, 'fit', fits[3], '.RData'))

posterior = cbind(postBeta0[-(1:burn),k], 
                  foreach(i = 1:p, .combine = 'cbind') %do% { postBeta[[i]][-(1:burn),k] },
                  postTau[-(1:burn),], postZeta[-(1:burn),k], postEta[-(1:burn)], postPhi[-(1:burn)])


sumNHPPSE = data.frame(
  row.names = row.names.se,
  par = par.se,
  truth = format(round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k], trpar$zeta[k], trpar$eta, trpar$phi), 1), nsmall = 1),
  dummytruth = round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k], trpar$zeta[k], ifelse(is.na(trpar$eta), 0, trpar$eta), ifelse(is.na(trpar$phi), 0, trpar$phi)), 1),
  mean = c(apply(posterior[,1:(1+p+1+p)], 2, bmmean), NA, NA, apply(posterior[,-(1:(1+p+1+p))], 2, bmmean)),
  dummyhpd = c(apply(posterior[,1:(1+p+1+p)], 2, hpd), NA, NA, apply(posterior[,-(1:(1+p+1+p))], 2, hpd)),
  hpd1 = c(apply(posterior[,1:(1+p+1+p)], 2, hpd1), NA, NA, apply(posterior[,-(1:(1+p+1+p))], 2, hpd1)),
  hpd2 = c(apply(posterior[,1:(1+p+1+p)], 2, hpd2), NA, NA, apply(posterior[,-(1:(1+p+1+p))], 2, hpd2))) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  select(c(par, truth, mean, hpd, sig))





# -----------------------------------------------------------------------------=
# (iv) NHPP + GP + SE ----
# -----------------------------------------------------------------------------=
load(paste0('sim2/fit/sim', datai, 'fit', fits[4], '.RData'))

posterior = cbind(postBeta0[-(1:burn),k], 
                  foreach(i = 1:p, .combine = 'cbind') %do% { postBeta[[i]][-(1:burn),k] },
                  postTau[-(1:burn),], postDelta[-(1:burn),k], postZeta[-(1:burn),k], postEta[-(1:burn)], postPhi[-(1:burn)])


sumNHPPGPSE = data.frame(
  par = par.se,
  truth = format(round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k], trpar$zeta[k], trpar$eta, trpar$phi), 1), nsmall = 1),
  dummytruth = round(c(trpar$beta0[k], trpar$beta[,k], trpar$tau, trpar$delta[k], trpar$zeta[k], ifelse(is.na(trpar$eta), 0, trpar$eta), ifelse(is.na(trpar$phi), 0, trpar$phi)), 1),
  mean = c(apply(posterior, 2, bmmean)),
  dummyhpd = c(apply(posterior, 2, hpd)),
  hpd1 = c(apply(posterior, 2, hpd1)),
  hpd2 = c(apply(posterior, 2, hpd2))) %>% 
  mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>% 
  mutate(hpd = ifelse( (hpd1 <= dummytruth) & (dummytruth <= hpd2), paste0('\\textbf{', dummyhpd, '}'), dummyhpd)) %>%
  select(c(mean, hpd, sig))




# -----------------------------------------------------------------------------=
# Table
# -----------------------------------------------------------------------------=
cbind(sumNHPP, sumNHPPGP)
cbind(sumNHPPSE, sumNHPPGPSE)


cbind(sumNHPP, sumNHPPGP) %>% 
  xtable() %>% 
  print(booktabs = T, sanitize.text.function = function(x) {x})

cbind(sumNHPPSE, sumNHPPGPSE) %>% 
  xtable() %>% 
  print(booktabs = T, sanitize.text.function = function(x) {x})





