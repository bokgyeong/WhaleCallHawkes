rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
# library(ggh4x) # facet_grid2
library(spgs) # chisq.unif.test
library(batchmeans); library(foreach)
library(xtable)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


datasets = fits = c('NHPPbk', 'NHPPbkw', 'NHPPbwk', 'NHPPbkwk', 
                    'NHPPbkSE', 'NHPPbkwSE', 'NHPPbwkSE', 'NHPPbkwkSE')

# -----------------------------------------------------------------------------=
# Box plot of betas ----
# -----------------------------------------------------------------------------=
burn = 5000

data.beta = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim2/fit/sim', datasets[i], 'fit', fits[j], '.RData'))
    
    if(fits[j] %in% c('NHPPbwk', 'NHPPbwkSE')){
      p = ncol(postBeta)
      K = ncol(postBeta0)
      for(i in 1:p){
        for(k in 1:K){
          data.beta = rbind(data.beta, 
                            data.frame(data = datasets[i], fit = fits[j], Name = paste0('beta', i),
                                       HP = factor(k), value = postBeta[-(1:burn),i])) 
        }
      }
    } else {
      p = length(postBeta)
      K = ncol(postBeta0)
      for(i in 1:p){
        for(k in 1:K){
          data.beta = rbind(data.beta, 
                            data.frame(data = datasets[i], fit = fits[j], Name = paste0('beta', i),
                                       HP = factor(k), value = postBeta[[i]][-(1:burn),k])) 
        }
      }
    } 
  }
}

data.beta$data = factor(data.beta$data, levels = datasets, 
                       labels = c('(1) X*Bk', '(2) X*Bk + W', '(3) X*B + Wk', '(4) X*Bk + Wk',
                                  '(5) X*Bk + SE', '(6) X*Bk + W + SE', '(7) X*B + Wk + SE', '(8) X*Bk + Wk + SE'))
data.beta$fit = factor(data.beta$fit, levels = fits, 
                       labels = c('(1) X*Bk', '(2) X*Bk + W', '(3) X*B + Wk', '(4) X*Bk + Wk',
                                  '(5) X*Bk + SE', '(6) X*Bk + W + SE', '(7) X*B + Wk + SE', '(8) X*Bk + Wk + SE'))

save(data.beta, file = 'sim2/fig/sim2Beta.RData')


plot.beta.NHPPbk = data.beta %>% 
  filter(fit == 'X*Bk') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbk

ggsave(plot = plot.beta.NHPPbk, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbk_beta.eps')

plot.beta.NHPPbkw = data.beta %>% 
  filter(fit == 'X*Bk + W') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkw

ggsave(plot = plot.beta.NHPPbkw, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbkw_beta.eps')

plot.beta.NHPPbwk = data.beta %>% 
  filter(fit == 'X*B + Wk') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbwk

ggsave(plot = plot.beta.NHPPbwk, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbwk_beta.eps')

plot.beta.NHPPbkwk = data.beta %>% 
  filter(fit == 'X*Bk + Wk') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkwk

ggsave(plot = plot.beta.NHPPbkwk, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbkwk_beta.eps')


plot.beta.NHPPbkSE = data.beta %>% 
  filter(fit == 'X*Bk + SE') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkSE

ggsave(plot = plot.beta.NHPPbkSE, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbkSE_beta.eps')

plot.beta.NHPPbkwSE = data.beta %>% 
  filter(fit == 'X*Bk + W + SE') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkwSE

ggsave(plot = plot.beta.NHPPbkwSE, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbkwSE_beta.eps')

plot.beta.NHPPbwkSE = data.beta %>% 
  filter(fit == 'X*B + Wk + SE') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbwkSE

ggsave(plot = plot.beta.NHPPbwkSE, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbwkSE_beta.eps')

plot.beta.NHPPbkwkSE = data.beta %>% 
  filter(fit == 'X*Bk + Wk + SE') %>% 
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkwkSE

ggsave(plot = plot.beta.NHPPbkwkSE, width = 6.2, height = 3.6, filename = 'sim2/fig/ccbNHPPbkwkSE_beta.eps')

