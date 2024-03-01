rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(ggh4x) # facet_grid2
library(spgs) # chisq.unif.test
library(batchmeans); library(foreach)
library(xtable)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


# -----------------------------------------------------------------------------=
# Histogram ----
# -----------------------------------------------------------------------------=

data.num = c()
for(i in 1:length(datasets)){
  load(paste0('sim/data/sim', datasets[i], '1.RData'))
  data.num = rbind(data.num, data.frame(data = datasets[i], ts = data$ts / 60))
}

data.num$data = factor(data.num$data, levels = datasets, labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))

plot.num = data.num %>% 
  ggplot(aes(x = ts)) +
  # geom_histogram(binwidth = 3, color = "black", fill = "white") +
  geom_histogram(binwidth = 2, color = "black", fill = "white") +
  facet_wrap(~data, scale = 'free_x') +
  scale_x_continuous(breaks = seq(0, 5 * 24, by = 12)) +
  labs(x = 'Time (hour)', y = 'Number of calls')
plot.num

ggsave(plot = plot.num, width = 5.5, height = 3.5, filename = 'sim/fig/usimNum.eps')



# -----------------------------------------------------------------------------=
# fitted XB vs true XB ----
# -----------------------------------------------------------------------------=

p = ncol(Xm)
betaInd = 1:p

data.xb = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim/data/sim', datasets[i], '1.RData'))
    load(paste0('sim/fit/sim', datasets[i], '1', fits[j], '.RData')) 
    trXB = Xm %*% trpar$beta 
    XB = postSamples[,betaInd] %*% t(Xm)
    XBci = t(sapply(1:ncol(XB), function(ii) HPDinterval(as.mcmc(XB[,ii]))[1:2]))
    data.xb = rbind(data.xb, data.frame(data = datasets[i], fit = fits[j], truth = trXB,
                                        ts = knts, mean = colMeans(XB), lb = XBci[,1], ub = XBci[,2]))
  }
}

data.xb$data = factor(data.xb$data, levels = datasets, labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.xb$fit = factor(data.xb$fit, levels = fits, labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))

plot.xb = data.xb %>% 
  ggplot(aes(x = ts)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) + 
  geom_line(aes(y = mean)) +
  geom_line(aes(y = truth), color = 'red', linetype = 2) +
  facet_wrap(~ data + fit, scales = 'free') +
  labs(x = expression(t), y = expression(X(t)*beta)) +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
plot.xb

ggsave(plot = plot.xb, width = 5.6, height = 5.5, device = cairo_ps,
       filename = 'sim/fig/usimXB.eps')




# -----------------------------------------------------------------------------=
# XB vs XBGP ----
# -----------------------------------------------------------------------------=

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

data.xbgp = c()
for(i in 1:length(datasets)){
  load(paste0('sim/data/sim', datasets[i], '1.RData'))
  
  trXB = Xm %*% trpar$beta 
  data.xbgp = rbind(data.xbgp, data.frame(data = datasets[i], Name = 'exp(XB)', ts = knts, lam0 = exp(trXB)))
  data.xbgp = rbind(data.xbgp, data.frame(data = datasets[i], Name = 'XB', ts = knts, lam0 = trXB))
  
  if(datasets[i] %in% c('LGCP', 'LGCPSE')){
    trXBGP = Xm %*% trpar$beta + trpar$Wm  
    data.xbgp = rbind(data.xbgp, data.frame(data = datasets[i], Name = 'exp(XB + W)', ts = knts, lam0 = exp(trXBGP)))
    data.xbgp = rbind(data.xbgp, data.frame(data = datasets[i], Name = 'XB + W', ts = knts, lam0 = trXBGP))
  }
}

data.xbgp$data = factor(data.xbgp$data, levels = datasets, labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.xbgp$Name = factor(data.xbgp$Name, levels = c('exp(XB)', 'exp(XB + W)', 'XB', 'XB + W'))

plot.xbgp = data.xbgp %>% 
  filter(Name %in% c('XB', 'XB + W')) %>% 
  ggplot(aes(x = ts)) +
  geom_line(aes(y = lam0, color = Name, linetype = Name)) +
  facet_wrap(~ data) +
  labs(x = 'Time (min)', y = '', color = '', linetype = '') +
  theme(legend.position = 'bottom')
plot.xbgp

ggsave(plot = plot.xbgp, width = 5.5, height = 3.5,
       filename = 'sim/fig/usimXBGP.eps')


plot.expxbgp = data.xbgp %>% 
  filter(Name %in% c('exp(XB)', 'exp(XB + W)')) %>% 
  ggplot(aes(x = ts)) +
  geom_line(aes(y = lam0, color = Name, linetype = Name)) +
  facet_wrap(~ data) +
  labs(x = 'Time (min)', y = '', color = '', linetype = '') +
  theme(legend.position = 'bottom')
plot.expxbgp


ggsave(plot = plot.expxbgp, width = 5.5, height = 3.5,
       filename = 'sim/fig/usimExpXBGP.eps')




# -----------------------------------------------------------------------------=
# mu(t) / lambda(t) ----
# -----------------------------------------------------------------------------=

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

data.mutrig = c()
for(i in 1:length(datasets)){
  load(paste0('sim/data/sim', datasets[i], '1.RData'))
  tsout = seq(0, maxT, length.out = 2000)[-1]
  indlam0 = sapply(1:length(tsout), function(i) which(knts >= tsout[i])[1] - 1 - 1)
  
  sourceCpp('src/RcppFtns.cpp')
  
  if(datasets[i] %in% c('LGCP', 'LGCPSE')){
    lam0m = as.vector(exp(Xm %*% trpar$beta + trpar$Wm))
  } else {
    lam0m = as.vector(exp(Xm %*% trpar$beta))
  }
  mut = compLam0(tsout, maxT, lam0m, indlam0, knts) 
    
  if(datasets[i] %in% c('LGCPSE', 'NHPPSE')){
    trigt = compSumAHiter(data$ts, tsout, trpar$alpha, trpar$eta)
  } else {
    trigt = 0
  }
    
  data.mutrig = rbind(data.mutrig, data.frame(data = datasets[i], ts = tsout, mu = mut, trig = trigt))
}

data.mutrig = data.mutrig %>% mutate(mubylam = mu / (mu + trig), loglambymu = log(mu + trig) - log(mu))
data.mutrig$data = factor(data.mutrig$data, levels = datasets, labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))

plot.mutrig = data.mutrig %>% 
  ggplot(aes(x = ts, y = mubylam)) +
  geom_line() +
  facet_wrap(~data) +
  labs(x = 'Time (min)', y = expression(mu(t)/lambda(t)))
plot.mutrig

ggsave(plot = plot.mutrig, width = 5.5, height = 3.5,
       filename = 'sim/fig/usimMuByLam.eps')


plot.loglambymu = data.mutrig %>% 
  ggplot(aes(x = ts, y = loglambymu)) +
  geom_line() +
  facet_wrap(~data) +
  labs(x = 'Time (min)', y = expression(log(lambda(t)/mu(t))))
plot.loglambymu

ggsave(plot = plot.loglambymu, width = 5.5, height = 3.5,
       filename = 'sim/fig/usimMuBLogLamByMu.eps')
