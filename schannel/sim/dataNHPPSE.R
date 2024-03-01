rm(list = ls())
library(Rcpp); library(RcppArmadillo); library(tidyverse)

source('src/RFtns.R')
sourceCpp('src/RcppFtns.cpp')


# design matrix ----
load('nopp/data/nopp.RData'); rm(data)
filename = 'sim/data/simNHPPSE.RData'

maxT = 5 * 24 * 60 # unit = min
sback = 20 # unit = min
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
rho = 60 # effective range is 3 hour

noise = data.frame(ts = knts) %>%
  left_join(noise)
noiseVar = as.vector(scale(noise$noise))

Xm = cbind(1,
           noiseVar,
           sin(2*pi*(knts + 15*60 + 1)/(8*60)),
           cos(2*pi*(knts + 15*60 + 1)/(8*60)),
           sin(2*pi*(knts + 15*60 + 1)/(12*60)),
           cos(2*pi*(knts + 15*60 + 1)/(12*60)),
           sin(2*pi*(knts + 15*60 + 1)/(24*60)),
           cos(2*pi*(knts + 15*60 + 1)/(24*60))) # should be greater than sback*4 
p = ncol(Xm)
# par(mfrow = c(2,3))
# for(i in 4:p){
#   plot(knts, Xm[,i], type = 'l')  
# }


# True parameters ----
# beta = c(-1.6, -0.5, 0.1, 0.3, 0.4, 0.0, -0.1, 0.1)
beta = c(-2.3, -1, 0.6, 0.8, 0.9, 0.5, -0.6, 0.6)
delta = NA
alpha = 0.34
eta = 0.51
trpar = list(beta = beta, delta = delta, alpha = alpha, eta = eta)

lam0m = as.vector( exp(Xm %*% beta) )


# Simulate data ----
data = simulateHawkes(alpha, eta, maxT, lam0m, knts)
ts = data$ts
branching = data$branching
n = length(ts)
sum(branching == 0)


# Approximation ----
indlam0 = sapply(1:n, function(i) which(knts >= ts[i])[1] - 1 - 1)
lam0 = compLam0(ts, maxT, lam0m, indlam0, knts)
intLam0 = compIntLam0(maxT, lam0m)
plot(knts, lam0m, type = 'l')
points(ts, lam0, col = 2)


# Save data ----
save(maxT, knts, m, rho, trpar, Xm, lam0m, data, indlam0, 
     file = filename)



