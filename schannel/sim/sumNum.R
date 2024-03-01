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

hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 2) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 2) }


datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


# -----------------------------------------------------------------------------=
# Compute
# -----------------------------------------------------------------------------=

burn = 10000

data.num = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim/data/sim', datasets[i], '.RData'))
    load(paste0('sim/num/sim', datasets[i], fits[j], 'num.RData'))

    ts = data$ts
    branching = data$branching

    data.num = rbind(data.num,
                     data.frame(
                       data = datasets[i], fit = fits[j],
                       Iteration = burn + 1:nrow(postNum),
                       numBack = postNum[,1], numSE = postNum[,2], numTotal = rowSums(postNum),
                       trnumBack = sum(branching == 0), trnumSE = sum(branching != 0), trnumTotal = length(ts)
                     ))
    
  }
}

save(data.num, file = 'sim/num/simNum.RData')




# -----------------------------------------------------------------------------=
# Posterior distribution of the number of calls ----
# -----------------------------------------------------------------------------=

load('sim/num/simNum.RData')

datalabs = fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+CC', '(iv) NHPP+GP+CC')

data.num$data = factor(data.num$data, levels = datasets, labels = datalabs)
data.num$fit = factor(data.num$fit, levels = fits, labels = fitlabs)


data.num.to = data.num %>% 
  group_by(data, fit, Iteration) %>% 
  summarise(tonumBack = sum(numBack), tonumSE = sum(numSE), tonumTotal = sum(numTotal),
            trtonumBack = sum(trnumBack), trtonumSE = sum(trnumSE), trtonumTotal = sum(trnumTotal))



## Models (iii) and (iv) ----

hist.total34 = data.num.to %>% 
  filter(fit == fitlabs[3:4]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumTotal, y=..density.., color = fit, fill = fit), alpha = 0.5) +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Total number of calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  # scale_color_grey() +
  # scale_fill_grey() +
  theme(legend.position = 'none')

hist.back34 = data.num.to %>% 
  filter(fit == fitlabs[3:4]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumBack, y=..density.., color = fit, fill = fit), alpha = 0.5) +
  geom_vline(aes(xintercept = trtonumBack), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of contact calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme(legend.position = 'none')

hist.SE34 = data.num.to %>% 
  filter(fit == fitlabs[3:4]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumSE, y=..density.., color = fit, fill = fit), alpha = 0.5) +
  geom_vline(aes(xintercept = trtonumSE), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of countercalls', y = 'Density', color = 'Fitting model', fill = 'Fitting model') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme(legend.position = 'bottom')


hist.all34 = ggarrange(hist.total34, hist.back34, hist.SE34, ncol = 1)
ggsave(plot = hist.all34, width = 9, height = 5.5, device = cairo_ps,
       filename = 'sim/fig/usimNum34.eps')




## Models (1) and (2) ----

hist.total12 = data.num.to %>% 
  filter(fit == fitlabs[1:2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumTotal, y=..density.., color = fit, fill = fit), alpha = 0.5) +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Total number of calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  # scale_color_grey() +
  # scale_fill_grey() +
  theme(legend.position = 'none')

hist.back12 = data.num.to %>% 
  filter(fit == fitlabs[1:2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumBack, y=..density.., color = fit, fill = fit), alpha = 0.5) +
  geom_vline(aes(xintercept = trtonumBack), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of contact calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme(legend.position = 'none')

hist.SE12 = data.num.to %>% 
  filter(fit == fitlabs[1:2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumSE, y=..density.., color = fit, fill = fit), alpha = 0.5) +
  geom_vline(aes(xintercept = trtonumSE), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of countercalls', y = 'Density', color = 'Fitting model', fill = 'Fitting model') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme(legend.position = 'bottom')


hist.all12 = ggarrange(hist.total12, hist.back12, hist.SE12, ncol = 1)
ggsave(plot = hist.all12, width = 9, height = 5.5, device = cairo_ps,
       filename = 'sim/fig/usimNum12.eps')




## (ii) NHPP + GP ----
hist.total2 = data.num.to %>% 
  filter(fit == fitlabs[2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumTotal, y=..density..), color="black", fill="white") +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Total number of calls', y = 'Density')

ggsave(plot = hist.total2, width = 9, height = 2,
       filename = 'sim/fig/msimNum2.eps')



## (i) NHPP ----
hist.total1 = data.num.to %>% 
  filter(fit == fitlabs[1]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumTotal, y=..density..), color="black", fill="white") +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Total number of calls', y = 'Density')

ggsave(plot = hist.total1, width = 9, height = 2,
       filename = 'sim/fig/msimNum1.eps')









## Total ----

hist.total = data.num.to %>% 
  ggplot() +
  geom_histogram(aes(x = tonumTotal, y=..density..), color="black", fill="white") +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data, fit), scales = 'free', nrow = 4) +
  labs(x = 'Total number of calls') +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
hist.total

ggsave(plot = hist.total, width = 5.6, height = 5.5,
       filename = 'sim/fig/msimNumTotal.eps')



## Background ----

hist.back = data.num.to %>% 
  filter(fit %in% fitlabs[3:4]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumBack, y=..density..), color="black", fill="white") +
  geom_vline(aes(xintercept = trtonumBack), color = 'red', linetype = 'dashed') +
  # facet_wrap(vars(data, fit), scales = 'free', ncol = 4) +
  facet_wrap(vars(data, fit), scales = 'free', nrow = 2, dir = 'v') +
  labs(x = 'Expected number of contact calls', y = 'Density') +
  # theme(strip.text = element_text(size = 7)) +
  # theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
  theme(strip.text = element_text(size = 9)) +
  theme(strip.text.x = element_text(margin = margin(.04, 0, .04, 0, "cm")))
hist.back

# ggsave(plot = hist.back, width = 5.6, height = 5.5,
#        filename = 'sim/fig/msimNumBack.eps')

ggsave(plot = hist.back, width = 9.5, height = 4,
       filename = 'sim/fig/msimNumBackM3M4.pdf')



## SE ----

hist.SE = data.num.to %>% 
  ggplot() +
  geom_histogram(aes(x = tonumSE, y=..density..), color="black", fill="white") +
  geom_vline(aes(xintercept = trtonumSE), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data, fit), scales = 'free', nrow = 4) +
  labs(x = 'Total number of countercalls') +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
hist.SE

ggsave(plot = hist.SE, width = 5.6, height = 5.5,
       filename = 'sim/fig/msimNumSE.eps')





