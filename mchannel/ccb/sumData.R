rm(list=ls())
library(sf); library(tigris)
library(tidyverse); library(readr)
library(Rcpp); library(RcppArmadillo)

options(tigris_use_cache = TRUE)

source('src/RFtns.R')
sourceCpp('src/RcppFtnsUnif.cpp')



counties_sf <- counties(cb = TRUE)
ccb = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')
coord = read_csv("ccb/data/CCB_2010.csv")

coord = coord %>% select(Latitude, Longitude) %>% add_column(HP = as.factor(1:10))

map.hp = coord %>% ggplot() + 
  geom_sf(data = ccb) + 
  geom_label(aes(x = Longitude, y = Latitude, label = HP)) +
  # geom_point(aes(x = Longitude, y = Latitude)) +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1))
map.hp


ggsave(plot = map.hp, width = 5, height = 3.7, filename = 'ccbUnif/fig/ccbHP.eps')




# =============================================================================-
# Number of calls per HP per hour
# =============================================================================-
load('ccbUnif/data/ccb.RData')
coord = read_csv("ccb/data/CCB_2010.csv")

data$marks = as.numeric(data$marks) + 1

indseq = seq(0+30, 7*24*60+30, by = 3*60)
data$ind = as.factor(sapply(1:nrow(data), function(i) which(data$ts[i] < indseq)[1]-1))

dat = data %>% group_by(ind, marks) %>% count()

f.dat = coord %>% 
  select(Site, Latitude, Longitude) %>% 
  right_join(dat, by = join_by(Site == marks)) %>% 
  arrange(ind, Site)


p.NumPerTime = f.dat %>% 
  filter(ind %in% as.factor(1:17)) %>%
  # filter(ind %in% as.factor(17:32)) %>%
  # filter(ind %in% as.factor(33:48)) %>%
  ggplot() + 
  geom_sf(data = ccb) + 
  # geom_point(aes(x = Longitude, y = Latitude)) +
  geom_point(aes(x = Longitude, y = Latitude, size = n), shape = 1) +
  facet_wrap(vars(ind)) +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p.NumPerTime

ggsave(plot = p.NumPerTime, width = 6.6, height = 5.5, filename = 'ccbUnif/fig/ccbNumPer3h.eps')





# =============================================================================-
# Event times ----
# =============================================================================-
load('ccbUnif/data/ccb.RData')

data$marks = as.numeric(data$marks) + 1

## Number of calls by each HP ----
table(data$marks) %>%
  t() %>% 
  xtable() %>% 
  print(booktabs = T, include.rownames = F)


## Call times of each HP ----
dat = data.frame(HP = data$marks, time = data$UTC)

shade = data.frame(dusk = c(as.POSIXlt("2010-04-02 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), 
                            seq.POSIXt(as.POSIXlt("2010-04-02 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-10 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')), 
                   dawn = c(as.POSIXlt("2010-04-02 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                            seq.POSIXt(as.POSIXlt("2010-04-03 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-11 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')),
                   top = Inf,
                   bottom = -Inf)

callTimes = dat %>% ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_point(aes(x = time, y = HP), size = 0.5) +
  scale_y_continuous(breaks = 1:10) +
  labs(x = 'Time', y = "Hydrophone")
callTimes

ggsave(plot = callTimes, device = cairo_ps, width = 6, height = 3, 
       file = 'ccbUnif/fig/ccbCallTimes.eps')




# =============================================================================-
# ambient noise ----
# =============================================================================-
noise = readRDS('ccbUnif/data/2010_noise_rms_wide.rds')
head(noise)

noise_ = noise
colnames(noise_) = c('UTC', 'day', paste0(1:10))
noise_ = noise_ %>% pivot_longer(!c(UTC, day), names_to = "HP", values_to = "noise")
noise_$HP = factor(noise_$HP, levels = paste0(1:10))

# noise_$UTC = format(noise_$UTC, format = '%Y-%m-%d %H:%M:%S')
# class(noise_$UTC)
# class(shade$dusk)

shade = data.frame(dusk = c(as.POSIXlt("2010-04-02 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), 
                            seq.POSIXt(as.POSIXlt("2010-04-02 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-10 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')), 
                   dawn = c(as.POSIXlt("2010-04-02 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                            seq.POSIXt(as.POSIXlt("2010-04-03 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-11 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')),
                   top = Inf,
                   bottom = -Inf)

plot.noise = noise_ %>%
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_line(aes(x = UTC, y = noise)) +
  facet_wrap(~ HP, nrow = 2) +
  labs(x = 'Time', y = 'Noise') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))
plot.noise

ggsave(plot = plot.noise, width = 9, height = 3, device = cairo_ps,
       filename = 'ccbUnif/fig/ccbNoise.eps')


