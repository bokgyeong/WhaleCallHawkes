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


datasets = fits = c('NHPPbk', 'NHPPbkw', 'NHPPbwk', 
                    'NHPPbkSE', 'NHPPbkwSE', 'NHPPbwkSE')

data.comp = c()
data.d = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim2/rtct/sim', datasets[i], 'fit', fits[j], 'RTCT.RData'))   
    data.comp = rbind(data.comp, data.frame(data = datasets[i], fit = fits[j], lb = postCompen[,1], Lam = postCompen[,2], ub = postCompen[,3]))
    data.d = rbind(data.d, data.frame(data = datasets[i], fit = fits[j], lb = postCompen[,4], d = postCompen[,5], ub = postCompen[,6]))
  }
}

data.comp$data = factor(data.comp$data, levels = datasets, 
                        labels = c('(1) X*Bk', '(2) X*Bk + W', '(3) X*B + Wk', 
                                   '(4) X*Bk + SE', '(5) X*Bk + W + SE', '(6) X*B + Wk + SE'))
data.comp$fit = factor(data.comp$fit, levels = fits, 
                       labels = c('(1) X*Bk', '(2) X*Bk + W', '(3) X*B + Wk', 
                                  '(4) X*Bk + SE', '(5) X*Bk + W + SE', '(6) X*B + Wk + SE'))
data.d$data = factor(data.d$data, levels = datasets, 
                     labels = c('(1) X*Bk', '(2) X*Bk + W', '(3) X*B + Wk', 
                                '(4) X*Bk + SE', '(5) X*Bk + W + SE', '(6) X*B + Wk + SE'))
data.d$fit = factor(data.d$fit, levels = fits, 
                    labels = c('(1) X*Bk', '(2) X*Bk + W', '(3) X*B + Wk', 
                               '(4) X*Bk + SE', '(5) X*Bk + W + SE', '(6) X*B + Wk + SE'))


# -----------------------------------------------------------------------------=
# Q-Q plot ----
# -----------------------------------------------------------------------------=

## Q-Q plot with uncertainty band ----
data.qq = data.d %>%
  group_by(data,fit) %>%
  reframe(order = order(d), Sample = d[order], lb = lb[order], ub = ub[order],
          Theoretical = log(length(d)) - log(length(d) - (1:length(d)-0.5)))


range.qq = data.qq %>%
  group_by(data,fit) %>% 
  reframe(x = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))),
          y = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))))
  # reframe(x = c(0, min(c(max(Sample), max(Theoretical)))),
  #         y = c(0, min(c(max(Sample), max(Theoretical)))))
  # reframe(x = c(0, max(c(Sample, Theoretical))), 
  #         y = c(0, max(c(Sample, Theoretical))))


plot.qqband = data.qq %>%
  ggplot(aes(x = Theoretical)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
  geom_point(aes(y = Sample), size = 0.8) +
  # geom_abline(intercept = 0, slope = 1) +
  geom_line(aes(x = x, y = y), range.qq) +
  facet_wrap(~data + fit, scales = 'free') +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
plot.qqband


ggsave(plot = plot.qqband, width = 12, height = 12, device = cairo_ps,
       filename = 'sim2/fig/sim2QQband6.eps')



## Q-Q plot ----
# plot.qq = data.d %>% 
  # ggplot(aes(sample = d)) +
  # stat_qq(distribution = qexp, size = 0.8) +
  # stat_qq_line(distribution = qexp) +
  # geom_abline(intercept = 0, slope = 1) +
plot.qq = data.qq %>% 
  ggplot(aes(x = Theoretical)) +
  geom_point(aes(y = Sample), size = 0.8) +
  geom_line(aes(x = x, y = y), range.qq) +
  facet_wrap(~ data + fit, scales = 'free') +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
plot.qq

ggsave(plot = plot.qq, width = 12, height = 12, device = cairo_ps,
       filename = 'sim2/fig/sim2QQ6.eps')




## mean squared difference between empirical and expected ----
data.qq %>%
  group_by(data, fit) %>% 
  summarise(msd = format(round(mean((Sample - Theoretical)^2), 3), nsmall = 3)) %>%
  pivot_wider(
    names_from = fit,
    values_from = msd
  ) %>%
  xtable() %>% 
  print(booktabs = T, include.rownames = F)
