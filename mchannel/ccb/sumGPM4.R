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

burn = 10000


fit = fits[4]


# -----------------------------------------------------------------------------=
# W(t) vs time ----
# -----------------------------------------------------------------------------=
load(paste0('ccb/fit/ccbfit', fit, '.RData'))


UTC = data.frame(ts = knts) %>%
  left_join(UTC)


W = postWm[-(1:burn),]
Wci = t(sapply(1:ncol(W), function(ii) HPDinterval(as.mcmc(W[,ii]))[1:2]))

data.gp = data.frame(fit = fit, Name = 'w(t)', ts = knts, UTC = UTC$UTC, value = colMeans(W), lb = Wci[,1], ub = Wci[,2])


shade = data.frame(dusk = c(as.POSIXlt("2010-04-02 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), 
                            seq.POSIXt(as.POSIXlt("2010-04-02 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-10 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')), 
                   dawn = c(as.POSIXlt("2010-04-02 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                            seq.POSIXt(as.POSIXlt("2010-04-03 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-11 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')),
                   top = Inf,
                   bottom = -Inf)


plot.gp = data.gp %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray70', alpha = 0.5) +
  geom_ribbon(aes(x = UTC, ymin = lb, ymax = ub), fill = "lightsteelblue", alpha = 0.6) +
  geom_line(aes(x = UTC, y = value)) +
  labs(x = 'Time', y = 'w(t)') +
  theme_bw() +
  # theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
  # theme(strip.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5)) +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")
plot.gp


ggsave(plot = plot.gp, width = 4, height = 2.8, 
       # device = cairo_ps,
       # filename = 'ccb/fig/ccbGP.eps')
       filename = 'ccb/fig/ccbGP.pdf')



# -----------------------------------------------------------------------------=
# delta_k * W(t) vs time ----
# -----------------------------------------------------------------------------=
load(paste0('ccb/fit/ccbfit', fit, '.RData'))

UTC = data.frame(ts = knts) %>%
  left_join(UTC)

data.dgp = c()
for(k in 1:K){
  W = postWm[-(1:burn),]
  delta = postDelta[-(1:burn),k]  
  data.dgp = rbind(data.dgp,
                   data.frame(fit = fit, MARU = factor(paste0('MARU ', k)), Name = 'GP', 
                              ts = knts, UTC = UTC$UTC, value = colMeans(delta * W))  
                   )
}

shade = data.frame(dusk = c(as.POSIXlt("2010-04-02 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), 
                            seq.POSIXt(as.POSIXlt("2010-04-02 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-10 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')), 
                   dawn = c(as.POSIXlt("2010-04-02 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                            seq.POSIXt(as.POSIXlt("2010-04-03 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'day', length.out = 8),
                            as.POSIXlt("2010-04-11 00:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')),
                   top = Inf,
                   bottom = -Inf)


plot.dgp = data.dgp %>% 
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray70', alpha = 0.5) +
  geom_line(aes(x = UTC, y = value)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Time', y = '') +
  facet_wrap(~MARU, nrow = 2) +
  theme_bw() +
  # theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
  # theme(strip.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5)) +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")
plot.dgp

ggsave(plot = plot.dgp, width = 10, height = 3.5, filename = 'ccb/fig/ccbDGP.pdf')


# plot.dgp = data.dgp %>% 
#   ggplot() +
#   geom_rect(data = shade,
#             aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
#             fill = 'gray70', alpha = 0.5) +
#   geom_line(aes(x = UTC, y = value, color = MARU)) +
#   labs(x = 'Time', y = '') +
#   theme_bw() +
#   # theme(strip.text.x = element_text(margin = margin(.02, 0, .02, 0, "cm"))) +
#   # theme(strip.text = element_text(size = 7)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5)) +
#   scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")
# plot.dgp


colMeans(postDelta[-(1:burn),])

