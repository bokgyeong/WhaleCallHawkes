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


datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
fits = c('NHPP', 'LGCP')

burn = 10000

data.loglik = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim/loglik/sim', datasets[i], fits[j], 'loglik.RData'))
    
    data.loglik = rbind(data.loglik, 
                        data.frame(data = datasets[i], fit = fits[j], Iteration = burn + 1:length(postLogLik), negTwoLoglik = -2*postLogLik))
  }
}


datalabs = fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+SE', '(iv) NHPP+GP+SE')
fitlabs = c('(i) NHPP', '(ii) NHPP+GP')

data.loglik$data = factor(data.loglik$data, levels = datasets, labels = datalabs)
data.loglik$fit = factor(data.loglik$fit, levels = fits, labels = fitlabs)


# -----------------------------------------------------------------------------=
# Posterior distribution of -2 * loglik ----
# -----------------------------------------------------------------------------=

burn = 10000

trace.negTwoLogLik = data.loglik %>% 
  filter(Iteration > burn) %>%
  ggplot(aes(x = Iteration, y = negTwoLoglik)) +
  geom_line() +
  facet_wrap(~data + fit, scales = 'free', nrow = 4) +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
  labs(x = 'Iteration', y = '-2logL')
trace.negTwoLogLik

ggsave(plot = trace.negTwoLogLik, width = 10, height = 6.5,
       filename = 'sim/fig/usimNegTwoLogLikTrace.eps')


hist.negTwoLogLik = data.loglik %>% 
  filter(Iteration > burn) %>%
  ggplot(aes(x = negTwoLoglik)) +
  geom_histogram(aes(y = ..density..), color="black", fill="white") +
  facet_wrap(~data + fit, scales = 'free', nrow = 4) +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
  labs(x = '-2logL', y = 'Density')
hist.negTwoLogLik

ggsave(plot = hist.negTwoLogLik, width = 10, height = 6.5,
       filename = 'sim/fig/usimNegTwoLogLikHist.eps')



data.loglik %>% 
  group_by(data, fit) %>% 
  summarise(Mean = paste0(format(round(bm(negTwoLoglik)[[1]], 1), nsmall = 1), 
                          ' (', format(round(bm(negTwoLoglik)[[2]], 1), nsmall = 1), ')')) %>% 
  pivot_wider(names_from = fit, values_from = Mean) %>% 
  as.data.frame() %>% 
  xtable() %>% 
  print(booktabs = F, include.colnames = F)
  


# -----------------------------------------------------------------------------=
# DIC ----
# -----------------------------------------------------------------------------=
data.dic$ExpDeviance = format(round(data.dic$ExpDeviance), nsmall = 0)
data.dic$EffNumPar = format(round(data.dic$EffNumPar), nsmall = 0)
data.dic$DIC = format(round(data.dic$DIC), nsmall = 0)

data.dic %>% 
  xtable() %>% 
  print(booktabs = T, include.rownames = F)



# -----------------------------------------------------------------------------=
# BIC ----
# -----------------------------------------------------------------------------=
# data.bic$ExpDeviance = format(round(data.bic$ExpDeviance), nsmall = 0)
# data.bic$NumPars = format(round(data.bic$NumPars), nsmall = 0)
# data.bic$BIC = format(round(data.bic$BIC), nsmall = 0)
# 
# data.bic %>% 
#   xtable() %>% 
#   print(booktabs = F, include.rownames = F)
# 
# 
# 
# cbind(data.dic, data.bic %>% select(-c(fit, ExpDeviance))) %>%
#   select(fit, ExpDeviance, EffNumPar, NumPars, DIC, BIC) %>% 
#   xtable() %>%
#   print(booktabs = F, include.rownames = F)
