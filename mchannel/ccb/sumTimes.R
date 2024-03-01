rm(list=ls())
library(sf); library(tigris)
library(tidyverse)
# options(tigris_use_cache = TRUE)

load('ccb/data/ccb.RData')

data.times = data.frame(time = data$ts, hp = data$marks)

plot.times = data.times %>% 
  ggplot(aes(x = time, y = hp, color = as.factor(hp+1))) +
  geom_point(size = 0.8) +
  labs(x = 'Time', y = '') +
  guides(color = guide_legend(title = "Hydrophone"))

ggsave(plot = plot.times, width = 6, height = 3.2, 
       filename = 'ccb/fig/ccbTimes.eps')
