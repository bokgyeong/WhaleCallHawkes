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


load(paste0('ccb/data/ccb.RData'))
fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

marks = data$marks + 1
K = nrow(distmat)

data.comp = c()
data.d = c()
for(j in 1:length(fits)){
  load(paste0('ccb/rtct/ccbfit', fits[j], 'rtct.RData'))
  print(paste0(fits[j], ': n = ', nrow(postCompen)))
  
  data.comp = rbind(data.comp, data.frame(fit = fits[j], ts = ts[1:nrow(postCompen)], hp = factor(marks[1:nrow(postCompen)]),
                                          lb = postCompen[,1], Lam = postCompen[,2], ub = postCompen[,3]))
  data.d = rbind(data.d, data.frame(fit = fits[j], ts = ts[1:nrow(postCompen)], hp = factor(marks[1:nrow(postCompen)]), 
                                    lb = postCompen[,4], d = postCompen[,5], ub = postCompen[,6]))
  # for(k in 1:K){
  #   data.comp.hp = rbind(data.comp.hp, data.frame(fit = fits[j], HP = factor(k), lb = postCompen[,6+(k-1)*3+1], Lam = postCompen[,6+(k-1)*3+2], ub = postCompen[,6+(k-1)*3+3]))
  #   data.d.hp = rbind(data.d.hp, data.frame(fit = fits[j], HP = factor(k), lb = postCompen[,6+10*3+(k-1)*3+1], d = postCompen[,6+10*3+(k-1)*3+2], ub = postCompen[,6+10*3+(k-1)*3+3]))
  # }
}


newlabs = c('(i) NHPP', '(ii) NHPP+GP', '(iii) NHPP+CC', '(iv) NHPP+GP+CC')

# data.d = data.d %>% mutate(u = 1 - exp(-d))
data.comp$fit = factor(data.comp$fit, levels = fits, labels = newlabs)
data.d$fit = factor(data.d$fit, levels = fits, labels = newlabs)


# data.comp.hp$fit = factor(data.comp.hp$fit, levels = fits, labels = newlabs)
# data.d.hp$fit = factor(data.d.hp$fit, levels = fits, labels = newlabs)



# -----------------------------------------------------------------------------=
# Q-Q plot ----
# -----------------------------------------------------------------------------=

## Q-Q plot with uncertainty band ----
data.qq = data.d %>%
  group_by(fit) %>%
  reframe(order = order(d), Sample = d[order], lb = lb[order], ub = ub[order], hp = hp[order],
          Theoretical = log(length(d)) - log(length(d) - (1:length(d)-0.5)))

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
               select(Sample, Theoretical) %>%
               max())
  
  for(j in groups[[i]]){
    dummy = rbind(dummy, data.frame(fit = newlabs[j], Sample = range_, Theoretical = range(range.qq$x)))
  }
}

new.data.qq = data.qq %>%
  group_by(fit) %>%
  mutate(newub = ifelse((fit %in% newlabs[groups[[1]]]) & (ub > max(dummy$Sample)), 
                        max(max(dummy$Sample), max(Theoretical)), ub))



plot.qqband = new.data.qq %>%
  ggplot(aes(x = Theoretical)) +
  geom_ribbon(aes(ymin = lb, ymax = newub), fill = "grey70", alpha = 0.8) +
  geom_point(aes(y = Sample), size = 0.5) +
  geom_line(aes(x = x, y = y), range.qq) +
  geom_blank(aes(x = Theoretical, y = Sample), dummy) +
  facet_wrap(~fit, scales = 'free', nrow = 1)
  
plot.qqband

ggsave(plot = plot.qqband, device = cairo_ps,
       width = 6.5, height = 1.9,
       filename = 'ccb/fig/ccbQQband.eps')





## Q-Q plot ----
# plot.qq = data.qq %>%
#   ggplot(aes(x = Theoretical)) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   geom_line(aes(x = x, y = y), range.qq) +
#   facet_wrap(~fit, scales = 'free', nrow = 2)
# plot.qq
# 
# ggsave(plot = plot.qq, device = cairo_ps,
#        width = 6.5, height = 3.5,
#        filename = 'ccb/fig/ccbwindQQ.eps')




## Q-Q plot by HP ----
# plot.qqband.hp = data.qq %>%
#   # filter(fit == '(7) Xk*B + Wk + SE') %>% 
#   filter(fit == '(6) Xk*B + Wk + SE') %>% 
#   ggplot(aes(x = Theoretical)) +
#   geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~hp, scales = 'free', nrow = 2)
# plot.qqband.hp
# 
# ggsave(plot = plot.qqband.hp, width = 7, height = 3.5, device = cairo_ps,
#        filename = 'ccb/fig/ccbwindQQbandHP.eps')


## mean squared difference between empirical and expected ----
data.qq %>%
  group_by(fit) %>% 
  summarise(msd = round(mean((Sample - Theoretical)^2), 2)) %>% 
  t() %>% 
  as.data.frame() %>% 
  xtable() %>% 
  print(booktabs = T, include.colnames = F) 








## Q-Q plot with colors ----
# data.qq.col = data.d
# 
# indday = seq(0, 9*60*24, by = 24*60)
# data.qq.col$indday = sapply(1:nrow(data.qq.col), function(i) factor(which(indday >= data.qq.col$ts[i])[1] - 1))
# 
# hr = 6
# indh = seq(0, 9*60*24, by = hr*60)
# data.qq.col$indh = sapply(1:nrow(data.qq.col), function(i) which(indh >= data.qq.col$ts[i])[1] - 1)
# data.qq.col$indh = data.qq.col$indh %% (24/hr)
# data.qq.col$indh[data.qq.col$indh == 0] = (24/hr)
# data.qq.col$indh = factor(data.qq.col$indh)
# 
# hr2 = 12
# indh2 = seq(0, 9*60*24, by = hr2*60)
# data.qq.col$indh2 = sapply(1:nrow(data.qq.col), function(i) factor(which(indh2 >= data.qq.col$ts[i])[1] - 1))
# 
# 
# data.qq.col = data.qq.col %>%
#   group_by(fit) %>%
#   reframe(order = order(d), 
#           hp = hp[order], indday = indday[order], indh = indh[order], indh2 = indh2[order],
#           Sample = d[order], lb = lb[order], ub = ub[order],
#           Theoretical = log(length(d)) - log(length(d) - (1:length(d)-0.5)),
#           slope = diff(quantile(Sample, probs = c(0.75, 0.25))) / diff(quantile(Theoretical, probs = c(0.75, 0.25))),
#           yint = (sum(quantile(Sample, probs = c(0.75, 0.25))) - slope * sum(quantile(Theoretical, probs = c(0.75, 0.25))))/2)
# 
# 
# plot.qq.hp = data.qq.col %>%
#   filter(fit == 'NHPP + Wk + SE') %>% 
#   ggplot(aes(x = Theoretical)) +
#   # ggplot(aes(x = Theoretical, label = hp, color = hp)) +
#   # geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   # geom_text(aes(y = Sample), hjust = -0.2, vjust = -0.2, size = 3) + 
#   geom_abline(aes(intercept = yint, slope = slope)) +
#   # scale_color_brewer(palette = 'Paired') +
#   # facet_wrap(~fit, scales = 'free') +
#   # labs(color = 'HP')
#   facet_wrap(~hp, scales = 'free')
# plot.qq.hp
# 
# ggsave(plot = plot.qq.hp, width = 6, height = 5, 
#        # filename = 'ccb/fig/ccbwindQQhp.eps')
#        filename = 'ccb/fig/ccbwindQQhp2.eps')
# 
# 
# plot.qq.day = data.qq.col %>%
#   filter(fit == 'NHPP + Wk + SE') %>% 
#   ggplot(aes(x = Theoretical)) +
#   # ggplot(aes(x = Theoretical, color = indday, label = indday)) +
#   # geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   # geom_text(aes(y = Sample), hjust = -0.2, vjust = -0.2, size = 3) + 
#   geom_abline(aes(intercept = yint, slope = slope)) +
#   # scale_color_brewer(palette = 'Paired') +
#   # facet_wrap(~fit, scales = 'free') +
#   # labs(color = 'Day')
#   facet_wrap(~indday, scales = 'free')
# plot.qq.day
# 
# ggsave(plot = plot.qq.day, width = 6, height = 5, 
#        # filename = 'ccb/fig/ccbwindQQday.eps')
#        filename = 'ccb/fig/ccbwindQQday2.eps')
# 
# 
# 
# plot.qq.h = data.qq.col %>%
#   filter(fit == 'NHPP + Wk + SE') %>% 
#   ggplot(aes(x = Theoretical)) +
#   # ggplot(aes(x = Theoretical, color = indh, label = indh)) +
#   # geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   # geom_text(aes(y = Sample), hjust = -0.2, vjust = -0.2, size = 3) + 
#   geom_abline(aes(intercept = yint, slope = slope)) +
#   # scale_color_brewer(palette = 'Paired') +
#   # facet_wrap(~fit, scales = 'free') +
#   # labs(color = paste0(hr, ' hour'))
#   facet_wrap(~indh, scales = 'free')
# plot.qq.h
# 
# ggsave(plot = plot.qq.h, width = 6, height = 5, 
#        # filename = paste0('ccb/fig/ccbwindQQ', hr, 'hr.eps'))
#        filename = paste0('ccb/fig/ccbwindQQ', hr, 'hr2.eps'))
# 
# 
# 
# plot.qq.h2 = data.qq.col %>%
#   filter(fit == 'NHPP + Wk + SE') %>% 
#   ggplot(aes(x = Theoretical)) +
#   # ggplot(aes(x = Theoretical, color = indh2, label = indh2)) +
#   # geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   # geom_text(aes(y = Sample), hjust = -0.2, vjust = -0.2, size = 3) + 
#   geom_abline(aes(intercept = yint, slope = slope)) +
#   # scale_color_brewer(palette = 'Paired') +
#   # facet_wrap(~fit, scales = 'free') +
#   # labs(color = paste0(hr2, ' hour'))
#   facet_wrap(~indh2, scales = 'free')
# plot.qq.h2




# -----------------------------------------------------------------------------=
# Q-Q plot by HP ----
# -----------------------------------------------------------------------------=
# 
# ## Q-Q plot with uncertainty band ----
# data.hp.qq = data.d.hp %>%
#   group_by(fit, HP) %>%
#   reframe(order = order(d), Sample = d[order], lb = lb[order], ub = ub[order],
#           Theoretical = log(length(d)) - log(length(d) - (1:length(d)-0.5)))
# 
# 
# plot.qq.hp.band = data.hp.qq %>%
#   filter(fit == '(7) Xk*B + Wk + SE') %>% 
#   ggplot(aes(x = Theoretical)) +
#   geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) +
#   geom_point(aes(y = Sample), size = 0.8) +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~HP, scales = 'free', nrow = 2)
# plot.qq.hp.band
# 
# 
# # ggsave(plot = plot.qq.hp.band, width = 7, height = 3.5, device = cairo_ps,
# #        filename = 'ccb/fig/ccbwindQQbyHPband.eps')
# 
# 
# 
# ## Q-Q plot ----
# plot.qq.hp = data.d.hp %>%
#   filter(fit == '(7) Xk*B + Wk + SE') %>% 
#   ggplot(aes(sample = d)) +
#   stat_qq(distribution = qexp, size = 0.8) +
#   # stat_qq_line(distribution = qexp) +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~ HP, scales = 'free', nrow = 2)
# plot.qq.hp
# 
# # ggsave(plot = plot.qq.hp, width = 7, height = 3.5, device = cairo_ps,
# #        filename = 'ccb/fig/ccbwindQQbyHP.eps')


