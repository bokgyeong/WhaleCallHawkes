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
hpd1 = function(x){ HPDinterval(as.mcmc(x))[1] }
hpd2 = function(x){ HPDinterval(as.mcmc(x))[2] }



load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/fit/ccbfitLGCPSE.RData'))

p = dim(beta)[1]
K = nrow(distmat)


burn = 10000 # decide based on the trace plots of -2logL

# =============================================================================-
# For each HP ----
# =============================================================================-

cov.names = c('Intercept', 'Noise',
              paste0('8h ', c('sine', 'cosine')), 
              paste0('12h ', c('sine', 'cosine')), 
              paste0('24h ', c('sine', 'cosine')))

row.names = c(cov.names, paste0(cov.names, ' variance'),
              'GP coefficient variance', 'GP coefficient')
row.names.se = c(row.names, 'SE intensity', 'Decay in time', 'Decay in space')



par = c(paste0('$\\beta_{', 0:p, ',k}$'), paste0('$\\tau_{', 0:p, '}$'), '\\tau_{\\delta}', paste0('\\delta_{k}'))
par.se = c(par, paste0('$\\alpha_{k}$'), '$\\eta$', '$\\phi$')






sumLGCPSE = list()
for(k in 1:K){
  posterior = cbind(postBeta0[-(1:burn),k], 
                    foreach(i = 1:p, .combine = 'cbind') %do% { postBeta[[i]][-(1:burn),k] },
                    postTau[-(1:burn),], postDelta[-(1:burn),k], postZeta[-(1:burn),k], postEta[-(1:burn)], postPhi[-(1:burn)])
  
  
  
  dummy = data.frame(
    row = row.names.se,
    par = par.se,
    Mean = c(apply(posterior, 2, bmmean)),
    CI = c(apply(posterior, 2, hpd)),
    hpd1 = c(apply(posterior, 2, hpd1)),
    hpd2 = c(apply(posterior, 2, hpd2))) %>%
    mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', ''))
  
  if(k %in% c(1, 6)){
    sumLGCPSE[[k]] = dummy %>% 
      select(c(row, par, Mean, CI, sig))    
    
  } else {
    sumLGCPSE[[k]] = dummy %>% 
      select(c(Mean, CI, sig))    
  }
}




## Tables ----

cbind(sumLGCPSE[[1]], sumLGCPSE[[2]], sumLGCPSE[[3]], sumLGCPSE[[4]], sumLGCPSE[[5]]) %>% 
  xtable() %>% 
  print(booktabs = T, include.rownames = F, sanitize.text.function = function(x) {x})

cbind(sumLGCPSE[[6]][-c(9, 11, 12),], sumLGCPSE[[7]][-c(9, 11, 12),], 
      sumLGCPSE[[8]][-c(9, 11, 12),], sumLGCPSE[[9]][-c(9, 11, 12),],
      sumLGCPSE[[10]][-c(9, 11, 12),]) %>% 
  xtable() %>% 
  print(booktabs = T, include.rownames = F, sanitize.text.function = function(x) {x})



r = 10
bmmean_r = function(x) { format(round(bm(x)$est, r), nsmall = r) }
hpd_r = function(x){ paste0('(', 
                          format(round(HPDinterval(as.mcmc(x))[1], r), nsmall = r), ', ',
                          format(round(HPDinterval(as.mcmc(x))[2], r), nsmall = r), ')') }
data.frame(
  row = c(paste0('MARU', 1:K, ' Excitement'), 'Temporal decay', 'Spatial decay'),
  # par = c(paste0('$\\alpha_{', 1:K, '}/\\eta$'), '$\\eta$', '$\\phi$'),
  # Mean = c(apply(cbind(postZeta[-(1:burn),]/postEta[-(1:burn)], postEta[-(1:burn)], postPhi[-(1:burn)]), 2, bmmean_r)),
  # CI = c(apply(cbind(postZeta[-(1:burn),]/postEta[-(1:burn)], postEta[-(1:burn)], postPhi[-(1:burn)]), 2, hpd_r))) %>% 
  par = c(paste0('$\\alpha_{', 1:K, '}'), '$\\eta$', '$\\phi$'),
  # Mean = c(apply(cbind(postZeta[-(1:burn),], postEta[-(1:burn)], postPhi[-(1:burn)]), 2, bmmean_r)),
  Mean = c(apply(cbind(postZeta[-(1:burn),], postEta[-(1:burn)], postPhi[-(1:burn)]), 2, median)),
  CI = c(apply(cbind(postZeta[-(1:burn),], postEta[-(1:burn)], postPhi[-(1:burn)]), 2, hpd_r))) %>%
  xtable() %>% 
  print(booktabs = T, include.rownames = F, sanitize.text.function = function(x) {x})




max(distmat)
min(distmat[upper.tri(distmat, diag = F)])

distmat = format(round(distmat, 1), nsmall = 1)
distmat[lower.tri(distmat, diag = T)] = NA
distmat = as.data.frame(distmat)
colnames(distmat) = rownames(distmat) = paste0(1:10)
distmat %>% 
  xtable() %>% 
  print(booktabs = T)




## Trace plots for alphastars ----
dat.post = data.frame(cbind(sapply(1:K, function(k) postZeta[-(1:burn),k]/postEta[-(1:burn)]), 
                            postEta[-(1:burn)], postPhi[-(1:burn)]))
aeplabs = c(paste0('alpha^{*}_', 1:10), 'eta', 'phi')
colnames(dat.post) = aeplabs
dat.post = dat.post %>% 
  mutate('Iteration' = (burn+1):nrow(postBeta0)) %>% 
  pivot_longer(!Iteration, names_to = 'Parameter', values_to = 'Value')
dat.post$Parameter = factor(dat.post$Parameter, levels = aeplabs)

### common axis ---
groups = list(c(1:10))
dummy = c()
for(i in 1:length(groups)){
  
  range_ = c(0, dat.post %>% 
               filter(Parameter %in% aeplabs[groups[[i]]]) %>% 
               select(Value) %>% 
               max())
  
  for(j in groups[[i]]){
    dummy = rbind(dummy, data.frame(Parameter = aeplabs[j], Value = range_, Iteration = dat.post$Iteration)) 
  }
}
dummy$Parameter = factor(dummy$Parameter, levels = aeplabs)

plot.post = dat.post %>% ggplot() +
  geom_line(aes(x = Iteration, y = Value)) +
  geom_blank(aes(x = Iteration, y = Value), dummy)+
  facet_wrap(~Parameter, scale = 'free', nrow = 3)
plot.post

ggsave(plot = plot.post, width = 8, height = 4, filename = 'ccbUnif/fig/ccbunifAlphaStarModel6.eps')


table(data$marks+1)


# =============================================================================-
# Across all HPs ----
# =============================================================================-
# load('ccbUnif/fit/ccbfitNHPPbkwSE.RData')
# 
# p = dim(beta)[1]
# K = nrow(distmat)
# 
# row.names.se = c(paste0('HP', 1:K, ' intercept'),
#                  foreach(i = 1:K, .combine = c) %do% { paste0('HP', i, c(' noise', ' wind speed', ' wind gust',
#                                                                          paste0(' 6h ', c('sine', 'cosine')), 
#                                                                          paste0(' 2h ', c('sine', 'cosine')), 
#                                                                          paste0(' 24h ', c('sine', 'cosine')))) },
#                  'Intercept variance',
#                  paste0(c('Noise', 'Wind speed', 'Wind gust',
#                           paste0('6h ', c('sine', 'cosine')), 
#                           paste0('2h ', c('sine', 'cosine')), 
#                           paste0('24h ', c('sine', 'cosine'))), ' variance'),
#                  'GP variance',
#                  paste0('HP', 1:K, ' SE intensity'), 
#                  'SE decay in time', 'SE decay in space')
# 
# par.se = c(sapply(0:p, function(i) paste0('$\\beta_{,', i, 1:K, '}$')),
#            '$\\tau_{00}$',
#            paste0('$\\tau_{', 0:(p-1), '}$'),
#            '$\\kappa$', 
#            paste0('$\\alpha^{*}_{', 1:K, '}$'), '$\\eta$', '$\\phi$')
# 
# 
# burn = 1000
# 
# posterior = cbind(postBeta0[-(1:burn),], 
#                   foreach(i = 1:p, .combine = 'cbind') %do% { postBeta[[i]][-(1:burn),] },
#                   postTau[-(1:burn),], postKappa[-(1:burn)],
#                   postZeta[-(1:burn),]/postEta[-(1:burn)], 
#                   postEta[-(1:burn)], postPhi[-(1:burn)])
# 
# 
# sumNHPPbkwSE = data.frame(
#   row = row.names.se,
#   par = par.se,
#   Mean = c(apply(posterior, 2, bmmean)),
#   CI = c(apply(posterior, 2, hpd)),
#   hpd1 = c(apply(posterior, 2, hpd1)),
#   hpd2 = c(apply(posterior, 2, hpd2))) %>%
#   mutate(sig = ifelse( (hpd1 > 0) | (0 > hpd2), '*', '')) %>%
#   select(c(row, par, Mean, CI, sig))
# 
# 
# 
# 
# 
# ### Table ----
# cbind(sumNHPPbkwSE[1:20,], sumNHPPbkwSE[-(1:20),]) %>% 
#   xtable() %>% 
#   print(booktabs = T, sanitize.text.function = function(x) {x}, include.rownames = FALSE)





