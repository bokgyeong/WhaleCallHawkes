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
# Create a dataframe ----
# =============================================================================-

burn = 10000
# addburn = 50000
addburn = 0
K = 10

load(paste0('ccb/data/ccb.RData'))
load(paste0('ccb/num/ccbfitLGCPSEnum.RData'))

ts = data$ts # unit is minutes
marks = data$marks

data.num = c()
for(k in 1:K){
  numBackatk = postNum[-(1:addburn),k]
  numSEatk = rowSums(sapply(1:K, function(m) postNum[-(1:addburn),K*m + k]))
  numTotalatk = numBackatk + numSEatk
  
  data.num = rbind(data.num,
                   data.frame(
                     Iteration = burn + addburn + 1:nrow(postNum[-(1:addburn),]),
                     MARU = as.factor(paste0('MARU ', k)),
                     numBack = numBackatk, numSE = numSEatk, numTotal = numTotalatk
                   ))
}



data.num.to = data.num %>% 
  group_by(Iteration, MARU) %>% 
  summarise(tonumBack = sum(numBack), tonumSE = sum(numSE), tonumTotal = sum(numTotal))



# -----------------------------------------------------------------------------=
## Joint distribution per MARU ----
# -----------------------------------------------------------------------------=

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


library(cowplot)
den.joint.all = plot_grid(den.joint[[1]],den.joint[[2]],den.joint[[3]],den.joint[[4]],den.joint[[5]],
                          den.joint[[6]],den.joint[[7]],den.joint[[8]],den.joint[[9]],den.joint[[10]],
                          align = "hv", nrow = 2, axis = "l")

den.joint.all.final = grid.arrange(arrangeGrob(den.joint.all, left = y.grob, bottom = x.grob))


ggsave(plot = den.joint.all.final, width = 9, height = 3.5, device = cairo_ps,
       filename = 'ccb/fig/ccbNumJoint.eps')
       # filename = 'ccb/fig/ccbNumJointFree.eps')





data.num.to %>% 
  filter(MARU %in% paste0('MARU ', 1:9)) %>% 
  ggplot(aes(x = tonumBack, y = tonumSE)) +
  stat_density_2d(aes(alpha = after_stat(level), color = MARU, fill = MARU), geom = "polygon")+                 
  stat_density_2d(aes(color = MARU), geom = "contour", alpha = .25)+
  scale_alpha(range = c(0.05, 0.2))+
  guides(alpha = FALSE) +
  theme_bw()+
  labs(x = "Expected number of contact calls", y = "Expected number of countercalls", color = '', fill = '')


