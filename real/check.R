rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)


fold = 'sim'
fold.data = 'data'
fold.fit = 'fit'
fold.fig = 'fig'

path.fit = paste0('/work/bk232/upcallHawkes/mchannel_rev/', fold, '/', fold.fit, '/')
path.fig = paste0(fold, '/', fold.fig, '/')
ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)


# datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
# datai = paste0('sim', datasets[3])
# fiti = datasets[3]

datasets = c('LGCPSE')
datai = paste0('sim', datasets[1])
fiti = 'LGCPSE'


load(paste0(fold, '/', fold.data, '/', datai, '.RData'))
load(paste0(path.fit, datai, '_', fiti, '.RData'))

dim(postBeta0)

p = nrow(beta)
K = nrow(distmat)


# Tilde ----
cname = c(paste0('betaTilde_', 0:p), 'log(deltaTilde)')
trp = data.frame(Parameter = cname, Value = c(mean(trpar$beta0), rowMeans(trpar$beta), log(mean(trpar$delta))))
colnames(postTilde) = cname
tc.tilde = as.data.frame(postTilde) %>% 
  add_column(Iteration = 1:nrow(postTilde)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.tilde




# beta0 ----
cname = paste0('beta_0', 1:ncol(postBeta0))
trp = data.frame(Parameter = cname, Value = trpar$beta0)
colnames(postBeta0) = cname
tc.beta0 = as.data.frame(postBeta0) %>% 
  add_column(Iteration = 1:nrow(postBeta0)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.beta0



# delta ----
cname = paste0('log(delta_', 1:ncol(postLogdelta), ')')
trp = data.frame(Parameter = cname, Value = log(trpar$delta))
colnames(postLogdelta) = cname
tc.delta = as.data.frame(postLogdelta) %>% 
  add_column(Iteration = 1:nrow(postLogdelta)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.delta


# W ----
dummyW = postWm[,1:9]
cname = paste0('W_', 1:ncol(dummyW))
trp = data.frame(Parameter = cname, Value = trpar$Wm[1:9])
colnames(dummyW) = cname
tc.Wm = as.data.frame(dummyW) %>% 
  add_column(Iteration = 1:nrow(dummyW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.Wm


# beta0 + delta*W ----
marks = data$ts[,2]
dummyBDW = sapply(1:ncol(dummyW), function(i) 
  postBeta0[,marks[i]] + exp(postLogdelta[,marks[i]]) * postWm[,i] )
cname = paste0('beta_{0,m_', 1:ncol(dummyW), '}+delta_{m_', 1:ncol(dummyW), '}W_', 1:ncol(dummyW))

trp = data.frame(Parameter = cname, 
                 Value = sapply(1:ncol(dummyW), function(i) trpar$beta0[marks[i]] + trpar$delta[marks[i]]*trpar$Wm[i]))
colnames(dummyBDW) = cname
tc.BDW = as.data.frame(dummyBDW) %>% 
  add_column(Iteration = 1:nrow(dummyBDW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.BDW



# delta*W ----
marks = data$ts[,2]
dummyDW = sapply(1:ncol(dummyW), function(i) exp(postLogdelta[,marks[i]]) * postWm[,i] )
cname = paste0('delta_{m_', 1:ncol(dummyW), '}W_', 1:ncol(dummyW))

trp = data.frame(Parameter = cname, Value = sapply(1:ncol(dummyW), function(i) trpar$delta[marks[i]]*trpar$Wm[i]))
colnames(dummyDW) = cname
tc.DW = as.data.frame(dummyDW) %>% 
  add_column(Iteration = 1:nrow(dummyBDW)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.DW





# alpha, eta, and phi ----
dummySE = cbind(postZeta, postEta, postPhi)
cname = c(paste0('alpha_', 1:K), 'eta', 'phi')
trp = data.frame(Parameter = cname, Value = c(trpar$zeta, trpar$eta, trpar$phi))
colnames(dummySE) = cname
tc.SE = as.data.frame(dummySE) %>% 
  add_column(Iteration = 1:nrow(dummySE)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  # filter(Iteration > 10000) %>%
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.SE

sum(data$branching == 0)
sum(postBranching[nrow(postBranching),] == 0)

# sourceCpp('src/RcppFtns.cpp')
# ind = which(postZeta[,1] > 5)
# for(i in ind){
#   print(compSumIntHl(1-1, data$ts[,1], data$ts[,2], maxT, distmat, postEta[i], postPhi[i]))
# }

# rtctAlphaSumIntHik(data$ts[1,1], data$ts[,1], data$ts[,2], distmat, postZeta, postEta, postPhi)



# beta ----
u = 2
dummyBeta = postBeta[[u]]
cname = paste0('beta_{', 1:K,  ',', u, '}')
trp = data.frame(Parameter = cname, Value = trpar$beta[u,])
colnames(dummyBeta) = cname
tc.beta = as.data.frame(dummyBeta) %>% 
  add_column(Iteration = 1:nrow(dummyBeta)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value') %>% 
  # filter(Iteration > 5000) %>% 
  ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_hline(aes(yintercept = Value), trp, color = 'red', linetype = 'dashed') +
  facet_wrap(~Parameter, scales = 'free')
tc.beta


# table(data$ts[,2])
