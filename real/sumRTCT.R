rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
get_meanCI = function(x){
  mean = mean(x)
  CI = as.vector(HPDinterval(as.mcmc(x)))
  return(c(mean, CI))
}

fold = 'real/'
path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.rtct = paste0(fold, 'rtct/')
path.fig = paste0(fold, 'fig/')
path.sum = paste0(fold, 'sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


datai = 'ccb'
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')


# =============================================================================-
# Compute RTCT ----
# =============================================================================-
load(paste0(path.data, datai, '.RData'))

ts = data$ts
marks = data$marks + 1
K = nrow(distmat)


data.d = c()
# load(paste0(path.sum, datai, 'rtct.RData'))

for(j in 1:length(fits)){
  fiti = fits[j]
  load(paste0(path.rtct, datai, '_', fiti, '_rtct.RData'))

  print(paste0(fiti, ': n = ', nrow(postCompen[[1]])))

  # order statistics
  orderstat = foreach(iter = 1:nrow(postCompen[[1]]), .combine = rbind) %do% {
    sort(postCompen[[1]][iter,])
  }

  # posterior mean and CI
  rtctres = foreach(i = 1:ncol(orderstat), .combine = rbind) %do% {
    get_meanCI(orderstat[, i])
  }

  data.d = rbind(
    data.d,
    data.frame(
      fit = fiti, ts = ts, hp = factor(marks),
      d = rtctres[,1], lb = rtctres[,2], ub = rtctres[,3]
    )
  )

  save(data.d, file = paste0(path.sum, datai, 'rtct.RData'))
}





# =============================================================================-
# Q-Q plot ----
# =============================================================================-
load(paste0(path.sum, datai, 'rtct.RData'))


newlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')

data.d$fit = factor(data.d$fit, levels = fits, labels = newlabs)


# data.comp.hp$fit = factor(data.comp.hp$fit, levels = fits, labels = newlabs)
# data.d.hp$fit = factor(data.d.hp$fit, levels = fits, labels = newlabs)



## Q-Q plot with uncertainty band ----
data.qq = data.d %>%
  group_by(fit) %>%
  mutate(Sample = d, Theoretical = log(length(d)) - log(length(d) - (1:length(d)-0.5)))

range.qq = data.qq %>%
  group_by(fit) %>% 
  reframe(x = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))),
          y = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))))


## range of x-axis and y-axis ----
# new.data.qq = data.qq %>% 
#   group_by(fit) %>% 
#   mutate(newub = ifelse(ub > max(Sample), max(max(Sample), max(Theoretical)), ub))


## common axis ---
groups = list(c(2, 3, 4))
dummy = c()
for(i in 1:length(groups)){
  
  range_ = c(0, data.qq %>%
               filter(fit %in% newlabs[groups[[i]]]) %>%
               ungroup() %>% 
               select(Sample, Theoretical) %>%
               max())
  
  for(j in groups[[i]]){
    dummy = rbind(dummy, data.frame(fit = newlabs[j], Sample = range_, Theoretical = range(range.qq$x)))
  }
}

new.data.qq = data.qq %>%
  group_by(fit) %>%
  mutate(
    newub = ifelse(
      (fit %in% newlabs[groups[[1]]]) & (ub > max(dummy$Sample)), 
      max(max(dummy$Sample), max(Theoretical)), ub
    )
  )



plot.qqband = new.data.qq %>%
  ggplot(aes(x = Theoretical)) +
  geom_ribbon(aes(ymin = lb, ymax = newub), fill = "grey75", alpha = 0.8) +
  geom_point(aes(y = Sample), size = 0.3) +
  geom_line(aes(x = x, y = y), range.qq) +
  geom_blank(aes(x = Theoretical, y = Sample), dummy) +
  facet_wrap(~fit, scales = 'free', nrow = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 7)
  )
  
plot.qqband

# ggsave(plot = plot.qqband, width = 6.5, height = 1.9,
#        filename = paste0(path.fig, 'ccbQQband.pdf'))




# =============================================================================-
# mean squared difference between empirical and expected ----
# =============================================================================-

data.qq %>%
  group_by(fit) %>% 
  summarise(msd = round(mean((Sample - Theoretical)^2), 3)) %>% 
  t() %>% 
  as.data.frame() %>% 
  xtable() %>% 
  print(booktabs = T, include.colnames = F) 




