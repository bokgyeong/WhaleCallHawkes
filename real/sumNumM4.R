rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable); library(cowplot)
library(sf); library(tigris)


fold = 'real/'
path.data = paste0(fold, 'data/')
path.num = paste0(fold, 'num/')
path.fig = paste0(fold, 'fig/')
path.sum = paste0(fold, 'sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


datai = 'ccb'
fiti = 'LGCPSE'


path.r = paste0(fold, 'src/RFtns.R')
path.cpp = paste0(fold, 'src/RcppFtns.cpp')

# =============================================================================-
# Within- and cross-recorder counter-calls ----
# =============================================================================-
load(paste0(path.data, datai, '.RData'))
load(paste0(path.num, datai, '_', fiti, '_num.RData'))

K = nrow(distmat)
ts = data$ts # unit is minutes
marks = data$marks


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
  dplyr::select(Received, meanBack, meanE) %>% 
  mutate(propE = meanE/sum(meanBack+meanE))



## plot ----
options(tigris_use_cache = TRUE)
counties_sf <- counties(cb = TRUE)
ccb = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')

coord = read_csv(paste0(path.data, "CCB_2010.csv"))
coord = coord %>% dplyr::select(Latitude, Longitude) %>% add_column(MARU = paste0(1:10))

data.wce = data.num %>% 
  dplyr::select(Received, meanWE, meanCE) %>% 
  pivot_longer(!Received, names_to = 'Name', values_to = 'Value')

data.wce$Name = factor(data.wce$Name, levels = c('meanWE', 'meanCE'), labels = c('(a) Within-MARU', '(b) Cross-MARU'))

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
plot.wce

# ggsave(plot = plot.wce, width = 7.3, height = 2.5, filename = paste0(path.fig, datai, 'WithinCrossEx.pdf'))



# =============================================================================-
# Joint posterior distribution of contact calls and coutercalls ----
# =============================================================================-
load(paste0(path.data, datai, '.RData'))
load(paste0(path.num, datai, '_', fiti, '_num.RData'))

K = nrow(distmat)
ts = data$ts # unit is minutes
marks = data$marks

data.num = c()
for(k in 1:K){
  numBackatk = postNum[,k]
  numSEatk = rowSums(sapply(1:K, function(m) postNum[,K*m + k]))
  numTotalatk = numBackatk + numSEatk
  
  data.num = rbind(data.num,
                   data.frame(
                     Iteration = 1:nrow(postNum[,]),
                     MARU = as.factor(paste0('MARU ', k)),
                     numBack = numBackatk, numSE = numSEatk, numTotal = numTotalatk
                   ))
}

data.num.to = data.num %>% 
  group_by(Iteration, MARU) %>% 
  summarise(tonumBack = sum(numBack), tonumSE = sum(numSE), tonumTotal = sum(numTotal))


range.back = range(data.num.to$tonumBack)
range.se = range(data.num.to$tonumSE)

den.joint = list()
for(k in 1:K){
  den.joint[[k]] = data.num.to %>% 
    filter(MARU == paste0('MARU ', k)) %>%
    ggplot(aes(x = tonumBack, y = tonumSE)) +
    stat_density_2d(aes(alpha = after_stat(level)), geom = "polygon", show.legend = FALSE)+                 
    stat_density_2d(geom = "contour", color = "black", alpha = .25, show.legend = FALSE)+
    scale_alpha(range = c(0.05,0.2))+
    guides(colour = FALSE, alpha = FALSE)+
    theme_bw()+
    facet_wrap(~MARU) +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 400)) +
    theme(axis.title = element_blank())
}

y.grob <- textGrob("Expected number of countercalls", rot=90)
x.grob <- textGrob("Expected number of contact calls")

den.joint.all = plot_grid(den.joint[[1]],den.joint[[2]],den.joint[[3]],den.joint[[4]],den.joint[[5]],
                          den.joint[[6]],den.joint[[7]],den.joint[[8]],den.joint[[9]],den.joint[[10]],
                          align = "hv", nrow = 2, axis = "l")

den.joint.all.final = grid.arrange(arrangeGrob(den.joint.all, left = y.grob, bottom = x.grob))


# ggsave(plot = den.joint.all.final, width = 9, height = 3.5, 
#        filename = paste0(path.fig, datai, 'NumJoint.pdf'))

