rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable); library(cowplot)
library(sf); library(tigris)

fold = 'sim/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.num = paste0(fold, 'num/')
path.fig = paste0(fold, 'fig/')
path.sum = paste0(fold, 'sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


# datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
datasets = c('NHPP')
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')
# =============================================================================-
# Compute expected numbers ----
# =============================================================================-

K = 10

data.num = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    datai = paste0('sim', datasets[i])
    fiti = fits[j]

    load(paste0(path.data, datai, '.RData'))
    load(paste0(path.num, datai, '_', fiti, '_num.RData'))

    ts = data$ts[,1]
    marks = data$ts[,2]
    branching = data$branching

    trnumTotal = table(marks)
    trnumBack = table(marks[branching == 0])
    if(datasets[i] %in% datasets[1:2] ){
      trnumSE = rep(0, K)
    } else {
      trnumSE = table(marks[branching != 0])
    }


    for(k in 1:K){
      numBackatk = postNum[,k]

      if( fits[j] %in% fits[1:2] ){
        numSEatk = rep(0, nrow(postNum))
      } else {
        numSEatk = rowSums(sapply(1:K, function(m) postNum[,K*m + k]))
      }
      numTotalatk = numBackatk + numSEatk

      data.num = rbind(
        data.num,
        data.frame(
          data = datasets[i], fit = fits[j],
          Iteration = 1:nrow(postNum),
          MARU = as.factor(k),
          numBack = numBackatk, numSE = numSEatk, numTotal = numTotalatk,
          trnumBack = trnumBack[k], trnumSE = trnumSE[k], trnumTotal = trnumTotal[k]
        )
      )
    }
  }
}

save(data.num, file = paste0(path.sum, 'simNum.RData'))




# =============================================================================-
# Histogram ----
# =============================================================================-

load(paste0(path.sum, 'simNum.RData'))

# datlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')
datlabs = c('(i) NHPP')
fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')

data.num$data = factor(data.num$data, levels = unique(data.num$data), labels = datlabs)
data.num$fit = factor(data.num$fit, levels = unique(data.num$fit), labels = fitlabs)


data.num.to = data.num %>% 
  group_by(data, fit, Iteration) %>% 
  summarise(tonumBack = sum(numBack), tonumSE = sum(numSE), tonumTotal = sum(numTotal),
            trtonumBack = sum(trnumBack), trtonumSE = sum(trnumSE), trtonumTotal = sum(trnumTotal))


## Models (i) and (ii) ----

hist.total12 = data.num.to %>% 
  filter(fit == fitlabs[1:2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumTotal, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Total number of calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(0.1, 0, 0.1, 0, unit = 'cm')),
  )

hist.back12 = data.num.to %>% 
  filter(fit == fitlabs[1:2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumBack, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_vline(aes(xintercept = trtonumBack), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of contact calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(0.1, 0, 0.1, 0, unit = 'cm')),
  )

hist.SE12 = data.num.to %>% 
  filter(fit == fitlabs[1:2]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumSE, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_vline(aes(xintercept = trtonumSE), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of countercalls', y = 'Density', color = 'Fitting model', fill = 'Fitting model') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(0.1, 0, 0.1, 0, unit = 'cm')),
  )


hist.all12 = ggarrange(hist.total12, hist.back12, hist.SE12, ncol = 1)
# ggsave(plot = hist.all12, width = 6, height = 4.8,
#        filename = paste0(path.fig, 'msimNum12.pdf'))





## Models (iii) and (iv) ----

hist.total34 = data.num.to %>% 
  filter(fit == fitlabs[3:4]) %>% 
  ggplot() +
  # geom_histogram(aes(x = tonumTotal, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_histogram(aes(x = tonumTotal, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_vline(aes(xintercept = trtonumTotal), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Total number of calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(0.1, 0, 0.1, 0, unit = 'cm')),
  )

hist.back34 = data.num.to %>% 
  filter(fit == fitlabs[3:4]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumBack, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_vline(aes(xintercept = trtonumBack), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of contact calls', y = 'Density') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(0.1, 0, 0.1, 0, unit = 'cm')),
  )

hist.SE34 = data.num.to %>% 
  filter(fit == fitlabs[3:4]) %>% 
  ggplot() +
  geom_histogram(aes(x = tonumSE, y=..density.., color = fit, fill = fit), alpha = 0.5, position = 'identity') +
  geom_vline(aes(xintercept = trtonumSE), color = 'red', linetype = 'dashed') +
  facet_wrap(vars(data), scales = 'free', nrow = 1) +
  labs(x = 'Number of countercalls', y = 'Density', color = 'Fitting model', fill = 'Fitting model') +
  scale_color_manual(values=c("gray50", "gray10")) +
  scale_fill_manual(values=c("gray70", "gray30")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(0.1, 0, 0.1, 0, unit = 'cm')),
  )


hist.all34 = ggarrange(hist.total34, hist.back34, hist.SE34, ncol = 1)
# ggsave(plot = hist.all34, width = 6, height = 4.8,
#        filename = paste0(path.fig, 'msimNum34.pdf'))




