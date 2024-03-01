rm(list=ls())
library(sf); library(tigris)
library(tidyverse)
options(tigris_use_cache = TRUE)

counties_sf <- counties(cb = TRUE)
ccb = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')
coord = read_csv("ccb/data/CCB_2010.csv")

coord = coord %>% select(Latitude, Longitude) %>% add_column(HP = as.factor(1:10))

map.hp = coord %>% ggplot() + 
  geom_sf(data = ccb) + 
  # geom_label(aes(x = Longitude, y = Latitude, label = HP)) +
  geom_point(aes(x = Longitude, y = Latitude)) +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1))




# =============================================================================-
# Number of calls per HP per hour
# =============================================================================-
wcall = readRDS('ccb/data/2010_unique_edges.rds')
coord = read_csv("ccb/data/CCB_2010.csv")

head(wcall)
head(coord)
coord = coord %>% mutate(maru_id = factor(Site - 1))

wcall[order(wcall$UTC),]

data = wcall %>% 
  select(UTC, maru_id) %>% 
  inner_join(coord) %>% 
  select(maru_id, UTC, Latitude, Longitude)

op <- options(digits.secs=3)
data$UTC = strptime(data$UTC, "%Y-%m-%d %H:%M:%OS", tz = "UTC")

data = data %>% 
  mutate(Month = format(UTC, '%m'), Day = format(UTC, '%d')) %>%
  filter(Month == '04', Day %in% paste0(0, 1:7)) %>%
  arrange(UTC)

std = as.numeric(strptime('2010-04-02 00:30:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC"))

data = data %>% 
  mutate(ts = ( as.numeric(UTC) - std ) / 60, marks = maru_id) %>% # unit = 1 min
  select(ts, marks, Latitude, Longitude) %>% 
  arrange(ts)

data$marks = as.factor(as.numeric(data$marks)-1)

indseq = seq(0+30, 7*24*60+30, by = 3*60)
data$ind = as.factor(sapply(1:nrow(data), function(i) which(data$ts[i] < indseq)[1]-1))

dat = data %>% group_by(ind, marks) %>% count()

f.dat = coord %>% 
  select(maru_id, Latitude, Longitude) %>% 
  right_join(dat, by = join_by(maru_id == marks)) %>% 
  arrange(ind, maru_id)


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

ggsave(plot = p.NumPerTime, width = 6.6, height = 5.5, filename = 'ccb/fig/ccbNumPer3h.eps')

