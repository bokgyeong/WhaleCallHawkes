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


load(paste0('ccb/data/ccb.RData'))

fits = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

ts = data$ts # unit is minutes
marks = data$marks
maxT = ceiling(max(ts))
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
rho_w = 60 # effective range is rho_w * 3

## covariates ----
noise = data.frame(ts = knts) %>%
  left_join(noise)
noise[,-1] = scale(noise[,-1])

K = nrow(distmat)
Xm = list()
for(i in 1:K){
  Xm[[i]] = cbind(noise[,i+1],
                  sin(2*pi*(knts + 30)/(8*60)),
                  cos(2*pi*(knts + 30)/(8*60)),
                  sin(2*pi*(knts + 30)/(12*60)),
                  cos(2*pi*(knts + 30)/(12*60)),
                  sin(2*pi*(knts + 30)/(24*60)),
                  cos(2*pi*(knts + 30)/(24*60))) # should be greater than sback*4
}
p = ncol(Xm[[1]])

burn = 10000


fit = fits[4]




# -----------------------------------------------------------------------------=
# Intenstiy functions ----
# -----------------------------------------------------------------------------=

load(paste0('ccb/lam/ccbfit', fit, 'lam.RData'))


dat.lam = c()
for(k in 1:K){
  
  dat.lam = rbind(dat.lam, 
                  data.frame(Name = 'Backgound', MARU = paste0('MARU ', k), Time = tsnew,
                             lb = postBack[,3*(k-1)+1], mean = postBack[,3*(k-1)+1+1], ub = postBack[,3*(k-1)+1+1+1]))
  
  dat.lam = rbind(dat.lam, 
                  data.frame(Name = 'Countercall', MARU = paste0('MARU ', k), Time = tsnew,
                             lb = postSE[,3*(k-1)+1], mean = postSE[,3*(k-1)+1+1], ub = postSE[,3*(k-1)+1+1+1]))
  
  dat.lam = rbind(dat.lam, 
                  data.frame(Name = 'Total call', MARU = paste0('MARU ', k), Time = tsnew,
                             lb = postLam[,3*(k-1)+1], mean = postLam[,3*(k-1)+1+1], ub = postLam[,3*(k-1)+1+1+1]))
  
  
}


dat.lam$MARU = factor(dat.lam$MARU, levels = paste0('MARU ', 1:K))

dat.lam = dat.lam %>% 
  mutate(UTC = as.POSIXct(Time * 60, origin = '2010-04-02 00:30:00.000', tz = "UTC", format = '%Y-%m-%d %H:%M:%S'))


shade = data.frame(dusk = c(as.POSIXlt("2010-04-02 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), 
                            seq.POSIXt(as.POSIXlt("2010-04-02 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-10 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')), 
                   dawn = c(as.POSIXlt("2010-04-02 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                            seq.POSIXt(as.POSIXlt("2010-04-03 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-11 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')),
                   top = Inf,
                   bottom = -Inf)


plot.back = dat.lam %>% 
  filter(Name == 'Backgound') %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  # geom_ribbon(aes(x = UTC, ymin = lb, ymax = ub), fill = "lightsteelblue", alpha = 0.6) +
  geom_line(aes(x = UTC, y = mean)) +
  facet_wrap(~MARU, nrow = 2) +
  labs(x = 'Date', y = 'Background intensity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))
plot.back

ggsave(plot = plot.back, width = 10, height = 3.5, device = cairo_ps,
       filename = 'ccb/fig/ccbLamBack.eps')




plot.both = dat.lam %>% 
  filter(Name != 'Total call') %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_line(aes(x = UTC, y = mean, color = Name, alpha = Name)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_alpha_manual(values = c(1, 0.4)) +
  facet_wrap(~MARU, nrow = 2) +
  labs(x = 'Date', y = 'Intensity', color = '', alpha = '') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5)) +
  theme(legend.position = 'bottom')
plot.both

ggsave(plot = plot.both, width = 10, height = 3.9, device = cairo_ps,
       filename = 'ccb/fig/ccbLamBackSE.eps')

  
  
  
plot.total = dat.lam %>% 
  filter(Name == 'Total call') %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  # geom_ribbon(aes(x = UTC, ymin = lb, ymax = ub), fill = "lightsteelblue", alpha = 0.6) +
  geom_line(aes(x = UTC, y = mean)) +
  facet_wrap(~MARU, nrow = 2) +
  labs(x = 'Date', y = 'Total call intensity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))
plot.total


ggsave(plot = plot.total, width = 10, height = 3.5, device = cairo_ps,
       filename = 'ccb/fig/ccbLamTotal.eps')



  
  
# -----------------------------------------------------------------------------=
# Rug plot 
# -----------------------------------------------------------------------------=

dat.raw = c()
for(k in 1:K){
  dat.raw = rbind(dat.raw,
                  data.frame(
                    MARU = factor(paste0('MARU ', k)),
                    Time = data$ts[(data$marks + 1) == k]
                  ))
}

head(dat.raw)

dat.raw = dat.raw %>% 
  mutate(UTC = as.POSIXct(Time * 60, origin = '2010-04-02 00:30:00.000', tz = "UTC", format = '%Y-%m-%d %H:%M:%S'))



plot.total.rug = dat.lam %>% 
  filter(Name == 'Total call') %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  # geom_ribbon(aes(x = UTC, ymin = lb, ymax = ub), fill = "lightsteelblue", alpha = 0.6) +
  geom_line(aes(x = UTC, y = mean)) +
  geom_rug(aes(x = UTC), dat.raw, length = unit(0.1, "npc")) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0), add = c(0.1, 0.02))) +
  facet_wrap(~MARU, ncol = 2) +
  labs(x = 'Date', y = 'Total call intensity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))

plot.total.rug  

ggsave(plot = plot.total.rug, width = 8, height = 8, device = cairo_ps,
       filename = 'ccb/fig/ccbLamTotalRug.eps')

