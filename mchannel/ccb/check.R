rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 2

load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/fit/ccbfit', datasets[runID],'.RData'))

p = nrow(beta)
K = nrow(distmat)



# beta0 ----
cname = paste0('beta_0', 1:ncol(postBeta0))
colnames(postBeta0) = cname
tc.beta0 = as.data.frame(postBeta0) %>% 
  add_column(Iteration = 1:nrow(postBeta0)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  facet_wrap(~Parameter, scales = 'free')
tc.beta0


# beta ----
u = 2
dummyBeta = postBeta[[u]]
cname = paste0('beta_{', 1:K,  ',', u, '}')
colnames(dummyBeta) = cname
tc.beta = as.data.frame(dummyBeta) %>% 
  add_column(Iteration = 1:nrow(dummyBeta)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  # filter(Iteration > 5000) %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  facet_wrap(~Parameter, scales = 'free')
tc.beta


# delta ----
cname = paste0('delta_', 1:ncol(postDelta))
colnames(postDelta) = cname
tc.delta = as.data.frame(postDelta) %>% 
  add_column(Iteration = 1:nrow(postDelta)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  facet_wrap(~Parameter, scales = 'free')
tc.delta


# W ----
dummyW = postWm[,1:9]
cname = paste0('W_', 1:ncol(dummyW))
colnames(dummyW) = cname
tc.Wm = as.data.frame(dummyW) %>% 
  add_column(Iteration = 1:nrow(dummyW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  facet_wrap(~Parameter, scales = 'free')
tc.Wm


# beta0 + delta*W ----
marks = data$ts[,2]
dummyBDW = sapply(1:ncol(dummyW), function(i) 
  postBeta0[,marks[i]] + postDelta[,marks[i]] * postWm[,i] )
cname = paste0('beta_{0,m_', 1:ncol(dummyW), '}+delta_{m_', 1:ncol(dummyW), '}W_', 1:ncol(dummyW))
colnames(dummyBDW) = cname
tc.BDW = as.data.frame(dummyBDW) %>% 
  add_column(Iteration = 1:nrow(dummyBDW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  facet_wrap(~Parameter, scales = 'free')
tc.BDW



# beta0 + delta*W ----
marks = data$ts[,2]
dummyDW = sapply(1:ncol(dummyW), function(i) postDelta[,marks[i]] * postWm[,i] )
cname = paste0('delta_{m_', 1:ncol(dummyW), '}W_', 1:ncol(dummyW))
colnames(dummyDW) = cname
tc.DW = as.data.frame(dummyDW) %>% 
  add_column(Iteration = 1:nrow(dummyBDW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  facet_wrap(~Parameter, scales = 'free')
tc.DW





# alpha, eta, and phi ----
dummySE = cbind(postZeta, postEta, postPhi)
cname = c(paste0('alpha_', 1:K), 'eta', 'phi')
colnames(dummySE) = cname
tc.SE = as.data.frame(dummySE) %>% 
  add_column(Iteration = 1:nrow(dummySE)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  filter(Iteration > 10000) %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  # facet_wrap(~Parameter, scales = 'free')
  facet_wrap(~Parameter)
tc.SE


quantile(postZeta[-(1:10000),3], probs = c(0.025, 0.975))

library(coda)
sapply(1:10, function(i) HPDinterval(as.mcmc(postZeta[-(1:10000), i])))
sort(table(data$marks+1))




# table(data$ts[,2])
