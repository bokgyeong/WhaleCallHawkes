rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
get_meanCI = function(x){
  mean = mean(x)
  CI = as.vector(HPDinterval(as.mcmc(x)))
  return(c(mean, CI))
}

path.data = paste0('data/')
path.fit = paste0('fit/')
path.rtct = paste0('rtct/')
path.fig = paste0('fig/')
path.sum = paste0('sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


# datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
datasets = c('NHPP')
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


path.r = paste0('src/RFtns.R')
path.cpp = paste0('src/RcppFtns.cpp')
# =============================================================================-
# Compute RTCT ----
# =============================================================================-

data.d = c()

for(ii in 1:length(datasets)){
  for(j in 1:length(fits)){
    datai = datasets[ii]
    fiti = fits[j]

    load(paste0(path.data, 'sim', datai, '.RData'))
    load(paste0(path.rtct, 'sim', datai, '_', fiti, '_rtct.RData'))

    print(paste0(datai, ' & ', fiti, ': n = ', nrow(postCompen[[1]])))

    ts = data$ts[,1]
    marks = data$ts[,2] + 1
    K = nrow(distmat)

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
        data = datai, fit = fiti, ts = ts, hp = factor(marks),
        d = rtctres[,1], lb = rtctres[,2], ub = rtctres[,3]
      )
    )

    save(data.d, file = paste0(path.sum, 'simrtct.RData'))
  }
}





# =============================================================================-
# Q-Q plot ----
# =============================================================================-

load(paste0(path.sum, 'simrtct.RData'))

# datlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')
datlabs = c('(i) NHPP')
fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+E', '(iv) NHPP+GP+E')

data.d$data = factor(data.d$data, levels = unique(data.d$data), labels = datlabs)
data.d$fit = factor(data.d$fit, levels = unique(data.d$fit), labels = fitlabs)



## Q-Q plot with uncertainty band ----
data.qq = data.d %>%
  group_by(data, fit) %>%
  mutate(Sample = d, Theoretical = log(length(d)) - log(length(d) - (1:length(d)-0.5)))


range.qq = data.qq %>%
  group_by(data,fit) %>% 
  reframe(x = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))),
          y = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))))
  # reframe(x = c(0, min(c(max(Sample), max(Theoretical)))),
  #         y = c(0, min(c(max(Sample), max(Theoretical)))))
  # reframe(x = c(0, max(c(Sample, Theoretical))), 
  #         y = c(0, max(c(Sample, Theoretical))))



## range of x-axis and y-axis ----
new.data.qq = data.qq %>% 
  group_by(data, fit) %>% 
  mutate(newub = ifelse(ub > max(Sample), max(max(Sample), max(Theoretical)), ub))



plot.qqband = new.data.qq %>%
  ggplot(aes(x = Theoretical)) +
  geom_ribbon(aes(ymin = lb, ymax = newub), fill = "grey70", alpha = 0.8) +
  geom_point(aes(y = Sample), size = 0.3) +
  # geom_abline(intercept = 0, slope = 1) +
  geom_line(aes(x = x, y = y), range.qq) +
  facet_wrap(~data + fit, scales = 'free') +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    # strip.background = element_blank(),
    strip.background = element_rect(color = NA),
    strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")),
  )
plot.qqband


# ggsave(plot = plot.qqband, width = 5.6, height = 5.5,
#        filename = paste0(path.fig, '/msimQQband.pdf'))



# =============================================================================-
# mean squared difference between empirical and expected ----
# =============================================================================-

data.qq %>%
  group_by(data, fit) %>% 
  summarise(msd = format(round(mean((Sample - Theoretical)^2), 4), nsmall = 4)) %>%
  # summarise(msd = mean((Sample - Theoretical)^2)) %>% 
  pivot_wider(names_from = fit, values_from = msd) %>% 
  as.data.frame() %>% 
  xtable() %>% 
  print(booktabs = T, include.rownames = F) 



