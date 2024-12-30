rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)

path.data = paste0('data/')
path.lam = paste0('lam/')
path.fig = paste0('fig/')
path.sum = paste0('sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 1
datai = paste0('sim', datasets[runID])
fiti = 'LGCPSE'


path.r = paste0('src/RFtns.R')
path.cpp = paste0('src/RcppFtns.cpp')

# =============================================================================-
# Load and summarize results ----
# =============================================================================-
load(paste0(path.data, datai, '.RData'))
load(paste0(path.lam, datai, '_', fiti, '_lam.RData'))


ts = data$ts[,1] # unit is minutes
marks = data$ts[,2]-1
soundspeed = 1.5 * 60 # km per min
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
# rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3

K = nrow(distmat)


dat.lam = c()
for(k in 1:K){
  
  dat.lam = rbind(
    dat.lam, 
    data.frame(
      Name = 'Backgound', MARU = paste0('MARU ', k), Time = tsnew,
      lb = postBack[,3*(k-1)+1], mean = postBack[,3*(k-1)+1+1], ub = postBack[,3*(k-1)+1+1+1]
    )
  )
  
  dat.lam = rbind(
    dat.lam, 
    data.frame(
      Name = 'Excitement', MARU = paste0('MARU ', k), Time = tsnew,
      lb = postSE[,3*(k-1)+1], mean = postSE[,3*(k-1)+1+1], ub = postSE[,3*(k-1)+1+1+1]
    )
  )
  
  dat.lam = rbind(
    dat.lam, 
    data.frame(
      Name = 'Total call', MARU = paste0('MARU ', k), Time = tsnew,
      lb = postLam[,3*(k-1)+1], mean = postLam[,3*(k-1)+1+1], ub = postLam[,3*(k-1)+1+1+1]
    )
  )
  
}

dat.lam$MARU = factor(dat.lam$MARU, levels = paste0('MARU ', 1:K))

dat.lam = dat.lam %>% 
  mutate(UTC = as.POSIXct(Time * 60, origin = '2010-04-02 00:30:00.000', tz = "UTC", format = '%Y-%m-%d %H:%M:%S'))


# =============================================================================-
# Figures ----
# =============================================================================-

library(suncalc)
dummy = getSunlightTimes(date = seq.Date(as.Date('2010-04-02'), as.Date('2010-04-06'), by = 1),
                         lat = 42, lon = -70.4, keep = c("sunrise", "sunset"), tz = "EST")

shade = data.frame(
  dusk = as.POSIXlt(c('2010-04-02 00:30:00', as.character(dummy$sunset)), tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
  dawn = as.POSIXlt(c(as.character(dummy$sunrise), '2010-04-07 00:00:00'), tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
  top = Inf,
  bottom = -Inf
)


## background and excitement process intensities ----
plot.both = dat.lam %>% 
  filter(Name != 'Total call') %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_line(aes(x = UTC, y = mean, color = Name, alpha = Name)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_alpha_manual(values = c(1, 0.4)) +
  facet_wrap(~MARU, ncol = 2) +
  labs(x = 'Date', y = 'Intensity', color = '', alpha = '') +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5)) +
  theme(
    axis.text= element_text(size = 7)
  ) +
  theme(legend.position = 'bottom')
plot.both

# ggsave(plot = plot.both, width = 7, height = 9, 
#        filename = paste0(path.fig, 'ccbLamBackSE.pdf'))


  
## total intensities ----
dat.raw = c()
for(k in 1:K){
  dat.raw = rbind(
    dat.raw,
    data.frame(
      MARU = factor(paste0('MARU ', k)),
      Time = ts[(marks + 1) == k]
    )
  )
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
  # theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))
  theme(
    axis.text= element_text(size = 7)
  )
plot.total.rug  

# ggsave(plot = plot.total.rug, width = 7, height = 9,
#        filename = paste0(path.fig, 'ccbLamTotalRug.pdf'))

