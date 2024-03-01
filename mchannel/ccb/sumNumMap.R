rm(list = ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
library(Rcpp); library(RcppArmadillo)
library(sf); library(tigris)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



# =============================================================================-
# Countercalls by excited by each MARU ----
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


options(tigris_use_cache = TRUE)
counties_sf <- counties(cb = TRUE)
ccb = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')

coord = read_csv("ccb/data/CCB_2010.csv")
coord = coord %>% dplyr::select(Latitude, Longitude) %>% add_column(MARU = paste0(1:10))


# num.all = table(marks+1)
# dat.num.se = matrix(0, K, K)
# for(l in 1:K){
#   for(k in 1:K){
#     dummy = exp(- phihat * distmat[l,k]) * astarhat[l] # by one call
#     dat.num.se[l,k] = dummy * as.numeric(num.all[l])
#   }
# }
# 
# colnames(dat.num.se) = paste0('HP', 1:10)
# 
# dat.num.se = dat.num.se %>% 
#   as.data.frame() %>% 
#   add_column(Parent = paste0('HP', 1:10)) %>% 
#   pivot_longer(cols = !Parent,
#                names_to = 'Offspring',
#                values_to = 'numCalls')


data.SE$Excite = factor(data.SE$Excite, levels = paste0(1:10), labels = paste0('MARU ', 1:10))


dat.num.se = coord %>% 
  select(MARU, Latitude, Longitude) %>% 
  right_join(data.SE, by = join_by(MARU == Received))

round(quantile(dat.num.se$meanSE, c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1)))

p.num.se = dat.num.se %>% 
  ggplot() + 
  geom_sf(data = ccb) + 
  geom_text(aes(x = Longitude, y = Latitude, label = ifelse(round(meanSE) == 0, NA, MARU), 
                 size = ifelse(round(meanSE) == 0, NA, round(meanSE)))) +
  # geom_point(aes(x = Longitude, y = Latitude, size = meanSE), shape = 1) +
  # geom_point(aes(x = Longitude, y = Latitude, size = ifelse(round(meanSE) == 0, NA, round(meanSE))), shape = 1) +
  facet_wrap(vars(Excite), nrow = 2) +
  labs(size = 'Expected number countercalls excited') +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1)) +
  # scale_size_continuous(breaks = c(1, 5, 10, 30, 50, 100, 200)) +
  scale_size_continuous(breaks = c(1, 10, 50, 100, 200)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom')
p.num.se

ggsave(plot = p.num.se, width = 8, height = 3.5, filename = 'ccb/fig/ccbNumSE.eps')




# =============================================================================-
# Within- and cross-monitor countercalls ----
# =============================================================================-

burn = 10000
K = 10

load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/num/ccbfitLGCPSEnum.RData'))

ts = data$ts # unit is minutes
marks = data$marks

head(postNum)


data.num = c()
for(k in 1:K){
  
  numBackatk = postNum[,k]
  numEatk = sapply(1:K, function(m) postNum[,K*m + k])
  numEfromk = postNum[,(K*k+1):(K*k + K)]
  
  df = data.frame(
    Received = as.factor(k),
    meanBack = mean(numBackatk), lbBack = quantile(numBackatk, 0.025), ubBack = quantile(numBackatk, 0.975),
    meanE = mean(rowSums(numEatk)), lbE = quantile(rowSums(numEatk), 0.025), ubE = quantile(rowSums(numEatk), 0.975),
    meanPBE = mean(numBackatk / rowSums(numEatk)), lbPBE = quantile(numBackatk/rowSums(numEatk), 0.025), ubPBE = quantile(numBackatk/rowSums(numEatk), 0.975),
    meanTo = mean(numBackatk + rowSums(numEatk)), lbTo = quantile(numBackatk + rowSums(numEatk), 0.025), ubTo = quantile(numBackatk + rowSums(numEatk), 0.975),
    meanWE = mean(numEatk[,k]), lbWE = quantile(numEatk[,k], 0.025), ubWE = quantile(numEatk[,k], 0.975),
    meanCE = mean(rowSums(numEatk[,-k])), lbCE = quantile(rowSums(numEatk[,-k]), 0.025), ubCE = quantile(rowSums(numEatk[,-k]), 0.975)
    # meanWE = mean(numEfromk[,k]), lbWE = quantile(numEfromk[,k], 0.025), ubWE = quantile(numEfromk[,k], 0.975),
    # meanCE = mean(rowSums(numEfromk[,-k])), lbCE = quantile(rowSums(numEfromk[,-k]), 0.025), ubCE = quantile(rowSums(numEfromk[,-k]), 0.975)
  )
  
  data.num = rbind(data.num, df)
}





## summary statistics ----

bm(rowSums(postNum))
HPDinterval(as.mcmc(rowSums(postNum)))
length(ts)

data.num %>% 
  dplyr::select(Received, meanTo, lbTo, ubTo) %>% 
  add_column(Obs = as.vector(table(marks+1)))


data.num %>% 
  dplyr::select(Received, meanBack, meanE)



## plot ----
options(tigris_use_cache = TRUE)
counties_sf <- counties(cb = TRUE)
ccb = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')

coord = read_csv("ccb/data/CCB_2010.csv")
coord = coord %>% dplyr::select(Latitude, Longitude) %>% add_column(MARU = paste0(1:10))


data.wce = data.num %>% 
  dplyr::select(Received, meanWE, meanCE) %>% 
  pivot_longer(!Received, names_to = 'Name', values_to = 'Value')

data.wce$Name = factor(data.wce$Name, levels = c('meanWE', 'meanCE'), labels = c('(a) Within-monitor countercall', '(b) Cross-monitor countercall'))

data.wce = coord %>% 
  dplyr::select(MARU, Latitude, Longitude) %>% 
  right_join(data.wce, by = join_by(MARU == Received))


plot.wce = data.wce %>% 
  ggplot() + 
  geom_sf(data = ccb) + 
  geom_point(aes(x = Longitude, y = Latitude, size = Value), shape = 1, color = 'red') +
  geom_text(aes(x = Longitude, y = Latitude, label = MARU), size = 3.5, color = 'gray50') +
  # scale_size_continuous(breaks = seq(0, ceiling(max(data.wce$Value)), length.out = 10)) +
  scale_size_continuous(range = c(0.5, 15), breaks = c(5, 10, 30, 50, 75, 100, 150, 200)) +
  facet_wrap(vars(Name)) +
  coord_sf(xlim = c(-70.7, -69.98), ylim = c(41.7, 42.1)) +
  labs(size = 'Expected number\nof countercalls') +
  theme_bw() +
  guides(size = guide_legend(ncol=2)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'right')

ggsave(plot = plot.wce, width = 7.3, height = 2.5, filename = 'ccb/fig/ccbWithinCrossEx.eps')

