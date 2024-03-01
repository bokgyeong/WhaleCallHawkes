rm(list = ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


bmmean = function(x) { format(round(bm(x)$est, 1), nsmall = 1) }
lb95 = function(x){ HPDinterval(as.mcmc(x), prob = 0.95)[1] }
ub95 = function(x){ HPDinterval(as.mcmc(x), prob = 0.95)[2] }
lb90 = function(x){ HPDinterval(as.mcmc(x), prob = 0.90)[1] }
ub90 = function(x){ HPDinterval(as.mcmc(x), prob = 0.90)[2] }


load('ccb/data/ccb.RData')

burn = 10000

# -----------------------------------------------------------------------------=
# Credible intervals ----
# -----------------------------------------------------------------------------=
load('ccb/fit/ccbfitLGCPSE.RData')


K = nrow(distmat)

postbeta1s = postBeta[[1]][-(1:burn),]



Predictors = paste0('MARU ', 1:K)

dat = data.frame(Mean = colMeans(postbeta1s), 
                 lb95 = apply(postbeta1s, 2, lb95), ub95 = apply(postbeta1s, 2, ub95),
                 lb90 = apply(postbeta1s, 2, lb90), ub90 = apply(postbeta1s, 2, ub90),
                 Predictor = Predictors)



dat$Predictor = factor(dat$Predictor, levels = rev(unique(dat$Predictor)))
# dat = dat %>% mutate(hasZero = ifelse(CIlow <= 0 & CIup >= 0, 'Include zero', 'Do not include zero'))


plot.ci = dat %>% 
  ggplot(aes(y = Predictor)) +
  # geom_errorbar(aes(xmin = lb95, xmax = ub95), width = 0.2) +
  geom_point(aes(x = Mean), size = 2) +
  geom_linerange(aes(xmin = lb95, xmax = ub95)) +
  geom_linerange(aes(xmin = lb90, xmax = ub90), linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.2) +
  # scale_x_break(c(-0.25, -0.13)) +
  # scale_x_break(c(0.21, 1.77), ticklabels = c(1.8)) +
  labs(x = 'HPD interval', y = '', title = '(a) Ambient noise') +
  guides(linetype = guide_legend(title=""),
         shape = guide_legend(title="")) +
  theme(legend.position = 'none') +
  theme(plot.title = element_text(vjust = 0.0001))
plot.ci




# ggsave(plot = plot.ci, width = 4.5, height = 3, file = 'noppUnif/fig/noppCI.eps')


# -----------------------------------------------------------------------------=
# Harmonic effects ----
# -----------------------------------------------------------------------------=

sort(table(data$marks+1))


load('ccb/fig/xbM4.RData')

data.xb$HP = factor(data.xb$HP, levels = paste(1:K), labels = paste0('MARU ', 1:K))

shade = data.frame(dusk = as.POSIXlt("2010-04-02 18:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), 
                   dawn = as.POSIXlt("2010-04-03 06:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                   top = Inf,
                   bottom = -Inf)

plot.xb = data.xb %>% 
  filter(UTC %in% seq.POSIXt(as.POSIXlt("2010-04-02 12:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'min', length.out = 24*60+1)) %>% 
  filter(Name == 'Xk*Bk with harmonics only') %>%
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_ribbon(aes(x = UTC, ymin = lb, ymax = ub), fill = "lightsteelblue", alpha = 0.6) +
  # geom_line(aes(x = UTC, y = ub), linetype = 'dotted') +
  # geom_line(aes(x = UTC, y = lb), linetype = 'dotted') +
  geom_line(aes(x = UTC, y = value)) +
  facet_wrap(~ HP, nrow = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Hour', y = 'Harmonic diurnal effects', title = '(b) Diurnal harmonic effect over time') +
  theme_bw() +
  scale_x_datetime(date_breaks = "4 hours", date_labels = "%H") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))




# -----------------------------------------------------------------------------=
# Combine figures ----
# -----------------------------------------------------------------------------=

plot.both = grid.arrange(plot.ci, plot.xb, nrow = 1, widths = c(1, 2.2))

ggsave(plot = plot.both, width = 9.6, height = 3.1, device = cairo_ps,
       filename = 'ccb/fig/ccbCInXB.eps')













