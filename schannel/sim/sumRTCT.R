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


datasets = fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')


data.comp = c()
data.d = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim/rtct/sim', datasets[i], fits[j], 'rtct.RData'))
    data.comp = rbind(data.comp, data.frame(data = datasets[i], fit = fits[j], ts = ts[1:nrow(postCompen)], 
                                            lb = postCompen[,1], t = postCompen[,2], ub = postCompen[,3]))
    data.d = rbind(data.d, data.frame(data = datasets[i], fit = fits[j], 
                                      lb = postCompen[,4], d = postCompen[,5], ub = postCompen[,6]))
  }
}

datalabs = fitlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+SE', '(iv) NHPP+GP+SE')


data.comp$fit = factor(data.comp$fit, levels = fits, labels = fitlabs)
data.comp$data = factor(data.comp$data, levels = datasets, labels = datalabs)
data.d$fit = factor(data.d$fit, levels = fits, labels = fitlabs)
data.d$data = factor(data.d$data, levels = datasets, labels = datalabs)




# -----------------------------------------------------------------------------=
# Q-Q plot ----
# -----------------------------------------------------------------------------=

data.qq = data.d %>% 
  group_by(data, fit) %>% 
  reframe(order = order(d), Sample = d[order], lb = lb[order], ub = ub[order],
            Theoretical = log(length(d)) - log(length(d) - (1:length(d) - 0.5)))


## reference line ----
range.qq = data.qq %>%
  group_by(data, fit) %>% 
  reframe(x = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))),
          y = c(0, ifelse( max(Sample) > max(Theoretical), min(c(max(Sample), max(Theoretical))), max(c(max(Sample), max(Theoretical))))))



## range of x-axis and y-axis ----
new.data.qq = data.qq %>% 
  group_by(data, fit) %>% 
  mutate(newub = ifelse(ub > max(Sample), max(max(Sample), max(Theoretical)), ub))

  
## Q-Q plot with uncertainty band ----
plot.qqband = new.data.qq %>%
  ggplot(aes(x = Theoretical)) +
  geom_ribbon(aes(ymin = lb, ymax = newub), fill = "grey70", alpha = 0.8) +
  geom_point(aes(y = Sample), size = 0.3) +
  geom_line(aes(x = x, y = y), range.qq) +
  facet_wrap(~data + fit, scales = 'free', nrow = 4) +
  # theme_bw() +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm")))
plot.qqband

ggsave(plot = plot.qqband, width = 5.6, height = 5.5, device = cairo_ps,
       filename = paste0('sim/fig/usimQQband.eps'))



## mean squared difference between empirical and expected ----
data.qq %>%
  group_by(data, fit) %>% 
  summarise(msd = format(round(mean((Sample - Theoretical)^2), 3), nsmall = 3)) %>%
  # summarise(msd = mean((Sample - Theoretical)^2)) %>% 
  pivot_wider(names_from = fit, values_from = msd) %>% 
  as.data.frame() %>% 
  xtable() %>% 
  print(booktabs = T, include.rownames = F) 



# -----------------------------------------------------------------------------=
# Compensator plot ----
# -----------------------------------------------------------------------------=
# plot.comp = data.comp %>% 
#   ggplot(aes(x = ts)) +
#   geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) + 
#   geom_line(aes(y = t)) +
#   geom_abline(intercept = 0, slope = 1, color = 'red') +
#   facet_wrap(~data + fit, scales = 'free') +
#   theme(strip.text = element_text(size = 7)) +
#   theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
#   labs(x = expression(t), y = expression(t^'*'))
# plot.comp
# 
# ggsave(plot = plot.comp, width = 5.6, height = 5.5, 
#        filename = 'sim/fig/usimComp.eps')



# -----------------------------------------------------------------------------=
# Chi-square test ----
# -----------------------------------------------------------------------------=
# data.pval$median = format(round(data.pval$median, 2), nsmall = 2)
# data.pval$lb = format(round(data.pval$lb, 2), nsmall = 2)
# data.pval$ub = format(round(data.pval$ub, 2), nsmall = 2)
# 
# data.pval = data.pval %>%
#   mutate(pval = ifelse(median >= 0.05, paste0('\\textbf{', median, '} (', lb, ', ', ub, ')'), 
#                        paste0(median, ' (', lb, ', ', ub, ')')))
# 
# 
# 
# test.rtct = t(data.frame(NHPP = data.pval %>% filter(data == 'NHPP') %>% select(pval) %>% unlist(),
#                          LGCP =data.pval %>% filter(data == 'NHPP + GP') %>% select(pval) %>% unlist(),
#                          NHPPSE = data.pval %>% filter(data == 'NHPP + SE') %>% select(pval) %>% unlist(),
#                          LGCPSE = data.pval %>% filter(data == 'NHPP + GP + SE') %>% select(pval) %>% unlist()))
# 
# colnames(test.rtct) = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE')
# rownames(test.rtct) = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE')
# 
# test.rtct %>% 
#   xtable() %>% 
#   print(booktabs = T, sanitize.text.function = function(x) {x})
