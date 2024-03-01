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


fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

data.mu = c()
data.trig = c()
data.lam = c()
data.mubylam = c()
for(j in 1:length(fits)){
  load(paste0('nopp/lam/nopp', fits[j], 'lam.RData'))
  
  if(fits[j] %in% c('LGCPSE', 'NHPPSE')){
    data.mu = rbind(data.mu, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,1], median = postLam[,2], ub = postLam[,3]))
    data.trig = rbind(data.trig, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,4], median = postLam[,5], ub = postLam[,6]))
    data.lam = rbind(data.lam, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,7], median = postLam[,8], ub = postLam[,9]))
    data.mubylam = rbind(data.mubylam, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,10], median = postLam[,11], ub = postLam[,12]))
    
  } else {
    data.mu = rbind(data.mu, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,1], median = postLam[,2], ub = postLam[,3]))
    data.trig = rbind(data.trig, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = 0, median = 0, ub = 0))
    data.lam = rbind(data.lam, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = postLam[,1], median = postLam[,2], ub = postLam[,3]))
    data.mubylam = rbind(data.mubylam, data.frame(fit = fits[j], ts = tsout[1:nrow(postLam)], lb = 1, median = 1, ub = 1))
  }
}



data.mu$fit = factor(data.mu$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.trig$fit = factor(data.trig$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.lam$fit = factor(data.lam$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))
data.mubylam$fit = factor(data.mubylam$fit, levels = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE'), labels = c('NHPP', 'NHPP + GP', 'NHPP + SE', 'NHPP + GP + SE'))



plot.mubylam = data.mubylam %>% 
  ggplot(aes(x = ts)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "grey70", alpha = 0.8) + 
  geom_line(aes(y = median)) +
  facet_wrap(~fit) +
  labs(x = 'Time (min)', y = expression(mu(t)/lambda(t)))
plot.mubylam  

ggsave(plot = plot.mubylam, width = 5.5, height = 3.5, device = cairo_ps,
       filename = 'nopp/fig/noppFittedMuByLam.eps')



