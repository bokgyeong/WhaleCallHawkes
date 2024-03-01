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
# XB and XB + W ----
# -----------------------------------------------------------------------------=
# load(paste0('ccb/fit/ccbfit', fit, '.RData'))
# 
# postBetanew = foreach(k = 1:K) %do% {
#   sapply(1:p, function(j) postBeta[[j]][,k])
# }
# 
# UTC = data.frame(ts = knts) %>%
#   left_join(UTC)
# 
# harmInd = c(2:7)
# 
# data.xb = c()
# for(k in 1:K){
#   XB = postBeta0[-(1:burn),k] + postBetanew[[k]][-(1:burn),] %*% t(Xm[[k]])
#   XBharm = postBetanew[[k]][-(1:burn),harmInd] %*% t(Xm[[k]][,harmInd])
#   XBW = XB + exp(postDelta[-(1:burn),k]) * postWm[-(1:burn),]
# 
#   XBci = t(sapply(1:ncol(XB), function(ii) HPDinterval(as.mcmc(XB[,ii]))[1:2]))
#   XBharmci = t(sapply(1:ncol(XB), function(ii) HPDinterval(as.mcmc(XBharm[,ii]))[1:2]))
#   XBWci = t(sapply(1:ncol(XBW), function(ii) HPDinterval(as.mcmc(XBW[,ii]))[1:2]))
# 
#   data.xb = rbind(data.xb, data.frame(fit = fit, Name = 'Xk*Bk', HP = k, ts = knts, UTC = UTC$UTC, value = colMeans(XB), lb = XBci[,1], ub = XBci[,2]))
#   data.xb = rbind(data.xb, data.frame(fit = fit, Name = 'Xk*Bk with harmonics only', HP = k, ts = knts, UTC = UTC$UTC, value = colMeans(XBharm), lb = XBharmci[,1], ub = XBharmci[,2]))
#   data.xb = rbind(data.xb, data.frame(fit = fit, Name = 'Xk*Bk + W', HP = k, ts = knts, UTC = UTC$UTC, value = colMeans(XBW), lb = XBWci[,1], ub = XBWci[,2]))
# }
# save(data.xb, file = 'ccb/fig/xbM4.RData')


