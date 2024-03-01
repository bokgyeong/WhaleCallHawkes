rm(list = ls())
library(Rcpp); library(RcppArmadillo)
library(tidyverse); library(readr)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')

# Note that the times of the call data and noise data are labeled UTC but they are actually EST

# =============================================================================-
# Call times ----
# =============================================================================-
wcall = rbind(
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090328_upcall_detection_log.txt") %>% mutate(ts = (BeginTimes + EndTimes) / 2 / 60),
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090329_upcall_detection_log.txt") %>% mutate(ts = 1*24*60 + (BeginTimes + EndTimes) / 2 / 60),
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090330_upcall_detection_log.txt") %>% mutate(ts = 2*24*60 + (BeginTimes + EndTimes) / 2 / 60),
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090331_upcall_detection_log.txt") %>% mutate(ts = 3*24*60 + (BeginTimes + EndTimes) / 2 / 60),
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090401_upcall_detection_log.txt") %>% mutate(ts = 4*24*60 + (BeginTimes + EndTimes) / 2 / 60),
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090402_upcall_detection_log.txt") %>% mutate(ts = 5*24*60 + (BeginTimes + EndTimes) / 2 / 60),
  read_table("nopp/data/NEFSC_SBNMS_200903_2_20090403_upcall_detection_log.txt") %>% mutate(ts = 6*24*60 + (BeginTimes + EndTimes) / 2 / 60)
)

data = wcall %>% 
  filter(is.na(Notes)) %>% 
  select(ts, Channel) %>% 
  arrange(ts)

plot(data$ts)
length(data$ts)
length(unique(data$ts))

data$ts = data$ts - 1



# =============================================================================-
# Noise variable ----
# =============================================================================-
noise = readRDS('nopp/data/2009_noise.rds')

std = as.numeric(strptime('2009-03-28 00:01:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC"))

noise = noise %>% 
  mutate(ts = ( as.numeric(UTC) - std) / 60, noise = RMS) %>% # unit = 1 min
  select(UTC, ts, noise) %>% 
  arrange(ts)



## Check if there are missing data in noise ----
ts = data$ts
maxT = ceiling(max(ts))

dummy = data.frame(ts = unique(c(0, seq(0, maxT, by = 1), maxT))) %>%
  left_join(noise)

which(is.na(dummy$noise)) # 900  9823
any(is.na(dummy %>% slice(901:9822) %>% select(noise) %>% unlist()))
dummy[c(901, 9822),]
# 2009-03-28 15:01:00 100.3619
# 2009-04-03 08:45:00  98.6896




# =============================================================================-
# Final dataset ----
# =============================================================================-

std.final = (as.numeric(strptime('2009-03-28 15:01:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC")) - std) / 60

data = data %>% 
  filter(ts >= (as.numeric(strptime('2009-03-28 15:01:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC")) - std) / 60) %>% 
  filter(ts <= (as.numeric(strptime('2009-04-03 08:45:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC")) - std) / 60) %>% 
  mutate(ts = ts - std.final)

plot(data$ts)
length(data$ts)
length(unique(data$ts))


noise = noise %>% 
  filter(ts >= (as.numeric(strptime('2009-03-28 15:01:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC")) - std) / 60) %>% 
  filter(ts <= (as.numeric(strptime('2009-04-03 08:45:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC")) - std) / 60) %>% 
  mutate(ts = ts - std.final)

dim(noise)
dim(unique(noise))
noise = unique(noise)


save(data, noise, file = 'nopp/data/nopp.RData')



# =============================================================================-
# Noise plot ----
# =============================================================================-
load('nopp/data/nopp.RData')
head(noise)

shade = data.frame(dusk = seq.POSIXt(as.POSIXlt("2009-03-28 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 6), 
                   dawn = seq.POSIXt(as.POSIXlt("2009-03-29 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 6),
                   top = Inf,
                   bottom = -Inf)


plot.noise = noise %>%
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_line(aes(x = UTC, y = noise)) +
  labs(x = 'Time', y = 'Noise') +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))

ggsave(plot = plot.noise, width = 5, height = 2.8, device = cairo_ps,
       filename = 'nopp/fig/nopp2009noise.eps')
