rm(list = ls())
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(tidyverse)
library(spgs) # chisq.unif.test

datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')
run_ID = 4

load(paste0('sim2/data/sim', datasets[run_ID], '.RData'))
load(paste0('sim2/fit/sim', datasets[run_ID], 'NHPPSE.RData'))
filename = paste0('sim2/rtct/sim', datasets[run_ID], 'NHPPSErtct.RData')
load(filename)

n = length(data$ts) 
save(n, k, devianceAtMean, postLogLik, file = filename)
