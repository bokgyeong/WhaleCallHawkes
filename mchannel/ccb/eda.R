rm(list = ls())
library(Rcpp); library(RcppArmadillo)
library(tidyverse); library(readr)
library(foreach)
library(sf); library(tigris)
options(tigris_use_cache = TRUE)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')


load(paste0('ccb/data/ccb.RData'))
K = nrow(distmat)

# -----------------------------------------------------------------------------=
# KDE ----
# -----------------------------------------------------------------------------=

kde = foreach(k = 1:K, .combine = cbind) %do% {
  ts = data$ts[ (data$marks + 1) == k ]
  density(ts, from = 0, to = ceiling(max(data$ts)))$y
}

kde = as.data.frame(kde)
colnames(kde) = factor(paste0('MARU', 1:10), levels = paste0('MARU', 1:10))

plot.kde = kde %>% 
  add_column(x = density(ts, from = 0, to = ceiling(max(data$ts)))$x) %>% 
  pivot_longer(!x, names_to = 'MARU', values_to = 'Value') %>% 
  mutate(MARU_ = factor(MARU, levels = paste0('MARU', 1:10))) %>% 
  ggplot() +
  geom_line(aes(x = x, y = Value)) +
  facet_wrap(~MARU_) +
  labs(x = 'Minute', y = 'Density', title = 'Kernal density estimation')

ggsave(plot.kde, width = 7, height = 5, file = 'ccb/fig_meeting/kde.pdf')



# -----------------------------------------------------------------------------=
# Cross-correlation ----
# -----------------------------------------------------------------------------=

# cormat = cor(kde)
# 
# ggcorrplot(cormat, type = "lower",
#            outline.col = "white",
#            ggtheme = ggplot2::theme_gray,
#            colors = c("#6D9EC1", "white", "#E46726"))


counties_sf <- counties(cb = TRUE)
ccb = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')
coord = read_csv("ccb/data/CCB_2010.csv")

coord = coord %>% dplyr::select(Latitude, Longitude) %>% add_column(HP = as.factor(1:10))

map.hp = coord %>% ggplot() + 
  geom_sf(data = ccb) + 
  geom_label(aes(x = Longitude, y = Latitude, label = HP)) +
  # geom_point(aes(x = Longitude, y = Latitude)) +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1))
map.hp


# 
# f.dat = coord %>% 
#   select(Site, Latitude, Longitude) %>% 
#   right_join(dat, by = join_by(Site == marks)) %>% 
#   arrange(ind, Site)


cormat = cor(kde)

# dat.cor = data.frame()
# ind = 1
# for(l in 1:(K-1)){
#   for(k in (l+1):K){
#     # if(cormat[l,k] > 0.3){
#       dat.cor = rbind(dat.cor,
#                       data.frame(pair = as.factor(ind), corr = cormat[l,k],
#                                  lonl = coord[l,2], latl = coord[l,1], 
#                                  lonk = coord[k,2], latk = coord[k,1]))
#       ind = ind + 1 
#     # }
#   }
# }


dat.cor = data.frame()
ind = 1
for(l in 1:K){
  for(k in 1:K){
    if(l != k){
    dat.cor = rbind(dat.cor,
                    data.frame(MARU = paste0('MARU ', l), MARUs = as.factor(k), 
                               corr = cormat[l,k],
                               lonl = coord[l,2], latl = coord[l,1], 
                               lonk = coord[k,2], latk = coord[k,1]))
    ind = ind + 1 
    }
  }
}

head(dat.cor)
dat.cor$MARU = factor(dat.cor$MARU, levels = paste0('MARU ', 1:10))

# library(RColorBrewer)
# library(scales)

plot.cor = coord %>% ggplot() + 
  geom_sf(data = ccb) + 
  geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude.1, yend = Latitude.1, linewidth = corr, color = corr), dat.cor) +
  geom_label(aes(x = Longitude, y = Latitude, label = HP), size = 2.8) +
  coord_sf(xlim = c(-70.5, -70.17), ylim = c(41.86, 42.04)) +
  scale_x_continuous(breaks = seq(-70.5, -70.17, by = 0.1)) +
  scale_linewidth(range = c(0.5, 3)) +
  labs(linewidth = 'Correlation coefficient', color = 'Correlation coefficient') +
  scale_colour_gradient2(guide = "legend") +
  facet_wrap(~MARU, nrow = 2) +
  theme(legend.position = 'bottom')


ggsave(plot.cor, width = 9.4, height = 4, 
       file = 'ccb/fig_meeting/cor.pdf')


library(scales)
plot.cor10 = coord %>% 
  ggplot() + 
  geom_sf(data = ccb) + 
  geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude.1, yend = Latitude.1, linewidth = corr, color = corr), 
               dat.cor %>% filter(MARU == 'MARU 10')) +
  geom_label(aes(x = Longitude, y = Latitude, label = HP), size = 2.8) +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1)) +
  scale_x_continuous(breaks = seq(-70.7, -69.98, by = 0.15)) +
  scale_linewidth(range = c(0.5, 3)) +
  labs(linewidth = 'Pairwise correlation', color = 'Pairwise correlation') +
  theme_bw() +
  # scale_colour_gradient2(high = muted('red'), low = muted('blue'), mid = 0, guide = "legend") +
  scale_colour_gradient2(high = 'red', low = 'blue', mid = 0, guide = "legend") +
  theme(legend.position = 'right')


ggsave(plot.cor10, width = 5.6, height = 3, 
       # file = 'ccb/fig_meeting/cor.pdf')
       device = cairo_ps, file = 'ccb/fig/cor10.eps')




plot.cor = coord %>% ggplot() + 
  geom_sf(data = ccb) + 
  geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude.1, yend = Latitude.1, linewidth = corr, color = corr), dat.cor) +
  geom_label(aes(x = Longitude, y = Latitude, label = HP), size = 2.8) +
  coord_sf(xlim = c(-70.5, -70.17), ylim = c(41.86, 42.04)) +
  scale_x_continuous(breaks = seq(-70.5, -70.17, by = 0.1)) +
  scale_linewidth(range = c(0.5, 3)) +
  labs(linewidth = 'Correlation coefficient', color = 'Correlation coefficient') +
  scale_colour_gradient2(guide = "legend") +
  facet_wrap(~MARU, nrow = 2) +
  theme(legend.position = 'bottom')


ggsave(plot.cor, width = 9.4, height = 4, 
       file = 'ccb/fig_meeting/cor.pdf')
# file = 'ccb/fig/cor.eps')





plot(cormat[upper.tri(cormat, diag = F)], distmat[upper.tri(distmat, diag = F)])
plot(distmat[upper.tri(distmat, diag = F)], cormat[upper.tri(cormat, diag = F)])



ind = which(distmat[upper.tri(distmat, diag = F)] < 9)
ind = which(distmat[upper.tri(distmat, diag = F)] < 9)

plot(distmat[upper.tri(distmat, diag = F)][ind], cormat[upper.tri(cormat, diag = F)][ind])


cormat[upper.tri(cormat, diag = F)]




