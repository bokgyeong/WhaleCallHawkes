rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

datai = 2
fiti = 2

load(paste0('sim/data/sim', datasets[datai], '.RData'))
load(paste0('sim/fit/sim', datasets[datai], fits[fiti],'.RData'))

p = ncol(Xm)



# beta ----
betaInd = 1:p
cname = paste0('beta_', 0:(p-1))
trp = data.frame(Parameter = cname, Value = trpar$beta)
dummyBeta = postSamples[,betaInd]
colnames(dummyBeta) = cname
tc.beta = as.data.frame(dummyBeta) %>% 
  add_column(Iteration = 1:nrow(dummyBeta)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  # filter(Iteration > 10000) %>%
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.beta


# delta ---
deltaInd = p+1
ts.plot(postSamples[,deltaInd])
abline(h = trpar$delta, col = 2)


# W ----
# indW = 1:9
indW = 10:18
dummyW = postWm[,indW]
cname = paste0('W_', 1:ncol(dummyW))
trp = data.frame(Parameter = cname, Value = trpar$Wm[indW])
colnames(dummyW) = cname
tc.Wm = as.data.frame(dummyW) %>% 
  add_column(Iteration = 1:nrow(dummyW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.Wm


# alpha and eta ----
alphaInd = p+2
etaInd = p+3
dummySE = postSamples[,c(alphaInd, etaInd)]
cname = c('alpha', 'eta')
trp = data.frame(Parameter = cname, Value = c(trpar$alpha, trpar$eta))
colnames(dummySE) = cname
tc.SE = as.data.frame(dummySE) %>% 
  add_column(Iteration = 1:nrow(dummySE)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  # filter(Iteration > 10000) %>%
  # filter(Iteration <= 10000) %>%
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.SE



