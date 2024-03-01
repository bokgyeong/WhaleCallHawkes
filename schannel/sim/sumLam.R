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

data.mu = c()
data.trig = c()
data.lam = c()
data.mubylam = c()
for(i in 1:length(datasets)){
  for(j in 1:length(fits)){
    load(paste0('sim/lam/sim', datasets[i], '1', fits[j], 'lam.RData'))
    
    if(fits[j] %in% c('LGCPSE', 'NHPPSE')){
      data.mu = rbind(data.mu, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,1], median = postLam[,2], ub = postLam[,3]))
      data.trig = rbind(data.trig, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,4], median = postLam[,5], ub = postLam[,6]))
      data.lam = rbind(data.lam, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,7], median = postLam[,8], ub = postLam[,9]))
      data.mubylam = rbind(data.mubylam, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,10], median = postLam[,11], ub = postLam[,12]))
      
    } else {
      data.mu = rbind(data.mu, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,1], median = postLam[,2], ub = postLam[,3]))
      data.trig = rbind(data.trig, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = 0, median = 0, ub = 0))
      data.lam = rbind(data.lam, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,1], median = postLam[,2], ub = postLam[,3]))
      data.mubylam = rbind(data.mubylam, data.frame(data = datasets[i], fit = fits[j], ts = tsout[1:nrow(postLam)], lb = 1, median = 1, ub = 1))
    }
  }
}


data.mu$fit = factor(data.mu$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.mu$data = factor(data.mu$data, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.trig$fit = factor(data.trig$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.trig$data = factor(data.trig$data, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.lam$fit = factor(data.lam$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.lam$data = factor(data.lam$data, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.mubylam$fit = factor(data.mubylam$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.mubylam$data = factor(data.mubylam$data, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))



plot.mubylam = data.mubylam %>% 
  ggplot(aes(x = ts)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) + 
  geom_line(aes(y = median)) +
  facet_wrap(~data + fit) +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
  labs(x = 'Time (min)', y = expression(mu(t)/lambda(t)))
plot.mubylam  

ggsave(plot = plot.mubylam, width = 5.6, height = 5.5, device = cairo_ps,
       filename = 'sim/fig/usimFittedMuByLam.eps')



