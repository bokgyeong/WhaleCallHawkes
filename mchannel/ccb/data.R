rm(list = ls())
library(Rcpp); library(RcppArmadillo)
library(tidyverse); library(readr)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')


# Note that the times of the call data and noise data are labeled UTC but they are actually EST


# -----------------------------------------------------------------------------=
# Up-call times ----
# -----------------------------------------------------------------------------=
wcall = readRDS('ccb/data/2010_unique_edges.rds')
coord = read_csv("ccb/data/CCB_2010.csv")

head(wcall)
head(coord)
unique(wcall$maru_id)


K = nrow(coord)
distmat = matrix(0, K, K)
library(geosphere)
for(k in 1:K){
  for(l in 1:K){
    if(k < l){
      distmat[k,l] = distmat[l,k] = distm(coord[k,c(7,6)], coord[l,c(7,6)]) / 1000
    }
  }
}


head(wcall)
head(coord)
coord = coord %>% mutate(maru_id = factor(Site-1))

head(wcall[order(wcall$UTC),])

data = wcall %>% 
  select(UTC, maru_id) %>% 
  inner_join(coord) %>% 
  select(maru_id, UTC, Latitude, Longitude)

op = options(digits.secs=3)
# data$UTC = strptime(data$UTC, "%Y-%m-%d %H:%M:%S", tz = "UTC")
data$UTC = strptime(data$UTC, "%Y-%m-%d %H:%M:%OS", tz = "UTC")

data = data %>% 
  # mutate(Month = format(UTC, '%m'), Day = format(UTC, '%d')) %>% 
  # filter(Month == '04', Day %in% paste0(0, 1:7)) %>% 
  arrange(UTC)

std = as.numeric(strptime('2010-04-02 00:30:00.000', "%Y-%m-%d %H:%M:%OS", tz = "UTC"))

data = data %>% 
  mutate(ts = ( as.numeric(UTC) - std ) / 60, marks = maru_id) %>% # unit = 1 min
  select(ts, marks, Latitude, Longitude, UTC) %>% 
  arrange(ts)

head(data)
dim(unique(data))
dim(data)

data = unique(data)
data$marks = as.numeric(data$marks)-1




# -----------------------------------------------------------------------------=
# Ambient noise ----
# -----------------------------------------------------------------------------=
noise = readRDS('ccb/data/2010_noise_rms_wide.rds')
head(noise)

head(noise$UTC)
head(as.numeric(noise$UTC))


UTC = noise %>%
  mutate(ts = ( as.numeric(UTC) - std) / 60) %>% # unit = 1 min
  arrange(ts) %>% 
  select(UTC, ts)

noise = noise %>%
  mutate(ts = ( as.numeric(UTC) - std) / 60) %>% # unit = 1 min
  select(-c(UTC, day)) %>% 
  arrange(ts)

head(noise)



# Save data ----
save(data, noise, UTC, distmat, file = paste0('ccb/data/ccb.RData'))






