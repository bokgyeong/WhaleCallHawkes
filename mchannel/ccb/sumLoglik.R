rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
bmmean = function(x) { format(round(bm(x)$est, 0), nsmall = 0) }
hpd = function(x){ paste0('(', 
                          format(round(HPDinterval(as.mcmc(x))[1], 1), nsmall = 1), ', ',
                          format(round(HPDinterval(as.mcmc(x))[2], 1), nsmall = 1), ')') }



load(paste0('ccb/data/ccb.RData'))
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
# fits = c('NHPP', 'LGCP', 'NHPPSE')

burn = 10000

data.loglik = c()
for(j in 1:length(fits)){
  load(paste0('ccb/loglik/ccbfit', fits[j], 'loglik.RData'))
  print(paste0(fits[j], ': ', length(postLogLik)))
  
  data.loglik = rbind(data.loglik, 
                      data.frame(fit = fits[j], Iteration = burn + 1:length(postLogLik), negTwoLoglik = -2*postLogLik))
}


# newlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+SE', '(iv) NHPP+GP+SE')
newlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+SE')

data.loglik$fit = factor(data.loglik$fit, levels = fits, labels = newlabs)



# -----------------------------------------------------------------------------=
# Posterior distribution of -2 * loglik ----
# -----------------------------------------------------------------------------=

burn = 10000

## Trace plot ----
trace.negTwoLogLik = data.loglik %>% 
  # filter(Iteration > burn) %>%
  ggplot(aes(x = Iteration, y = negTwoLoglik)) +
  geom_line() +
  facet_wrap(~fit, scales = 'free', nrow = 1) +
  # facet_wrap(~fit, nrow = 1) +
  labs(x = 'Iteration', y = '-2logL')
trace.negTwoLogLik


ggsave(plot = trace.negTwoLogLik, width = 10, height = 3.2, filename = 'ccb/fig/ccbNegTwoLogLikTrace.eps')



## Histogram ----

histall.negTwoLogLik = data.loglik %>% 
  filter(Iteration > burn) %>%
  ggplot(aes(x = negTwoLoglik)) +
  geom_histogram(aes(y = ..density..), bins = 30, color="black", fill="white") +
  facet_wrap(~fit, scales = 'free', nrow = 1) +
  labs(x = '-2logL', y = 'Density')
histall.negTwoLogLik

ggsave(plot = histall.negTwoLogLik, width = 10, height = 3.2, filename = 'ccb/fig/ccbNegTwoLogLikHistAll.eps')



groups = list(c(2, 6:8))
dummy = c()
for(i in 1:length(groups)){
  
  range_neg2logL = data.loglik %>% 
    filter(Iteration > burn) %>%
    filter(fit %in% newlabs[groups[[i]]]) %>% 
    select(negTwoLoglik) %>% 
    range()
  
  for(j in groups[[i]]){
    dummy = rbind(dummy, data.frame(fit = newlabs[j], negTwoLoglik = range_neg2logL)) 
  }
}

hist.negTwoLogLik = data.loglik %>% 
  filter(fit %in% newlabs[c(2, 6:8)]) %>% 
  filter(Iteration > burn) %>%
  ggplot(aes(x = negTwoLoglik)) +
  geom_histogram(aes(y = ..density..), bins = 40, color="black", fill="white") +
  facet_wrap(~fit, scales = 'free', nrow = 1) +
  geom_blank(aes(x = negTwoLoglik), data = dummy) +
  labs(x = '-2logL', y = 'Density')
hist.negTwoLogLik

ggsave(plot = hist.negTwoLogLik, width = 10, height = 2, filename = 'ccb/fig/ccbNegTwoLogLikHist.eps')




  



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
# 
# cbind(data.dic, data.bic %>% select(-c(fit, ExpDeviance))) %>%
#   select(fit, ExpDeviance, EffNumPar, NumPars, DIC, BIC) %>% 
#   xtable() %>%
#   print(booktabs = F, include.rownames = F)


