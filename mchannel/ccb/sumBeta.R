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


fits = c('NHPPbk', 'NHPPbkw', 'NHPPbwk', 'NHPPbkwk',
         'NHPPbkSE', 'NHPPbkwSE', 'NHPPbwkSE', 'NHPPbkwkSE')


# -----------------------------------------------------------------------------=
# Box plot of betas ----
# -----------------------------------------------------------------------------=
# burn = 1000
# 
# data.beta = c()
# for(j in 1:length(fits)){
#   load(paste0('ccbRSRwind/fit/ccbfit', fits[j], '.RData'))   
#   
#   if(fits[j] %in% c('NHPPbwk', 'NHPPbwkSE')){
#     p = ncol(postBeta)
#     K = ncol(postBeta0)
#     for(i in 1:p){
#       for(k in 1:K){
#         data.beta = rbind(data.beta, 
#                           data.frame(fit = fits[j], Name = paste0('beta', i),
#                                      HP = factor(k), value = postBeta[-(1:burn),i])) 
#       }
#     }
#   } else {
#     p = length(postBeta)
#     K = ncol(postBeta0)
#     for(i in 1:p){
#       for(k in 1:K){
#         data.beta = rbind(data.beta, 
#                           data.frame(fit = fits[j], Name = paste0('beta', i),
#                                      HP = factor(k), value = postBeta[[i]][-(1:burn),k])) 
#       }
#     }
#   } 
# }
# save(data.beta, file = 'ccbRSRwind/fig/ccbBeta.RData')




# -----------------------------------------------------------------------------=
# Plot ----
# -----------------------------------------------------------------------------=

load('ccbRSRwind/fig/ccbBeta.RData')

newlabs = c('(1) Xk*Bk', '(2) Xk*Bk+W', '(3) Xk*B+Wk', '(4) Xk*Bk+Wk',
            '(5) Xk*Bk+SE', '(6) Xk*Bk+W+SE', '(7) Xk*B+Wk+SE', '(8) Xk*Bk+Wk+SE')
data.beta$fit = factor(data.beta$fit, levels = fits, labels = newlabs)


# plot.beta.NHPPbk = data.beta %>% 
#   filter(fit == '(1) Xk*Bk') %>% 
#   ggplot(aes(x = HP, y = value)) +
#   geom_boxplot(outlier.size = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   facet_wrap(~Name, scales = 'free') +
#   labs(x = 'Hydrophone', y = '')
# plot.beta.NHPPbk
# 
# ggsave(plot = plot.beta.NHPPbk, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel1.eps')
# 
# plot.beta.NHPPbkw = data.beta %>% 
#   filter(fit == 'Xk*Bk+W') %>% 
#   ggplot(aes(x = HP, y = value)) +
#   geom_boxplot(outlier.size = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   facet_wrap(~Name, scales = 'free') +
#   labs(x = 'Hydrophone', y = '')
# plot.beta.NHPPbkw
# 
# ggsave(plot = plot.beta.NHPPbkw, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel2.eps')
# 
# plot.beta.NHPPbwk = data.beta %>% 
#   filter(fit == 'Xk*B+Wk') %>% 
#   ggplot(aes(x = HP, y = value)) +
#   geom_boxplot(outlier.size = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   facet_wrap(~Name, scales = 'free') +
#   labs(x = 'Hydrophone', y = '')
# plot.beta.NHPPbwk
# 
# ggsave(plot = plot.beta.NHPPbwk, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel3.eps')
# 
# plot.beta.NHPPbkwk = data.beta %>% 
#   filter(fit == 'Xk*Bk+Wk') %>% 
#   ggplot(aes(x = HP, y = value)) +
#   geom_boxplot(outlier.size = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   facet_wrap(~Name, scales = 'free') +
#   labs(x = 'Hydrophone', y = '')
# plot.beta.NHPPbkwk
# 
# ggsave(plot = plot.beta.NHPPbkwk, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel4.eps')
# 
# 
# plot.beta.NHPPbkSE = data.beta %>% 
#   filter(fit == 'Xk*Bk+SE') %>% 
#   ggplot(aes(x = HP, y = value)) +
#   geom_boxplot(outlier.size = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   facet_wrap(~Name, scales = 'free') +
#   labs(x = 'Hydrophone', y = '')
# plot.beta.NHPPbkSE
# 
# ggsave(plot = plot.beta.NHPPbkSE, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel5.eps')

plot.beta.NHPPbkwSE = data.beta %>%
  filter(fit == '(6) Xk*Bk+W+SE') %>%
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, scales = 'free') +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkwSE

ggsave(plot = plot.beta.NHPPbkwSE, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel6.eps')

# plot.beta.NHPPbwkSE = data.beta %>% 
#   filter(fit == 'Xk*B+Wk+SE') %>% 
#   ggplot(aes(x = HP, y = value)) +
#   geom_boxplot(outlier.size = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   facet_wrap(~Name, scales = 'free') +
#   labs(x = 'Hydrophone', y = '')
# plot.beta.NHPPbwkSE
# 
# ggsave(plot = plot.beta.NHPPbwkSE, width = 6.2, height = 3.6, filename = 'ccbRSRwind/fig/ccbBoxBetaModel7.eps')


data.beta$Name = factor(data.beta$Name, levels = paste0('beta', 1:7), labels = paste0('beta', 0:6))

plot.beta.NHPPbkwkSE = data.beta %>%
  filter(fit == '(8) Xk*Bk+Wk+SE') %>%
  ggplot(aes(x = HP, y = value)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Name, nrow = 2) +
  labs(x = 'Hydrophone', y = '')
plot.beta.NHPPbkwkSE

ggsave(plot = plot.beta.NHPPbkwkSE, width = 9, height = 4, filename = 'ccbRSRwind/fig/ccbBoxBetaModel8.eps')

