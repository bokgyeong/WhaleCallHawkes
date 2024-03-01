rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
# library(ggh4x) # facet_grid2
library(spgs) # chisq.unif.test
library(batchmeans); library(foreach)
library(xtable)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

hpd1 = function(x){ round(HPDinterval(as.mcmc(x))[1], 2) }
hpd2 = function(x){ round(HPDinterval(as.mcmc(x))[2], 2) }



# =============================================================================-
# K by K matrix ----
# =============================================================================-

burn = 10000
K = 10

load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/num/ccbfitLGCPSEnum.RData'))

ts = data$ts # unit is minutes
marks = data$marks

data.back = c()
data.SE = c()
for(k in 1:K){
  numBackatk = postNum[,k]
  numSEatk = sapply(1:K, function(m) postNum[,K*m + k])

  data.back = rbind(data.back,
                    data.frame(
                      Received = as.factor(k),
                      meanBack = mean(numBackatk), lbBack = quantile(numBackatk, 0.025), ubBack = quantile(numBackatk, 0.975)
                    ))
  for(l in 1:K){
    data.SE = rbind(data.SE,
                    data.frame(
                      Received = as.factor(k), Excite = as.factor(l),
                      meanSE = mean(numSEatk[,l]), lbSE = quantile(numSEatk[,l], 0.025), ubSE = quantile(numSEatk[,l], 0.975)
                    ))
  }
}



tab.SE = data.SE %>% 
  mutate(intmeanSE = format(round(meanSE, 0), nsmall = 0)) %>% 
  select(Received, Excite, intmeanSE) %>% 
  pivot_wider(names_from = Received, values_from = intmeanSE) %>% 
  select(-Excite) %>% 
  xtable() %>%
  print(booktabs = T)



