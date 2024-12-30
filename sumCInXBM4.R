rm(list=ls())
library(coda); library(tidyverse); library(egg); library(grid)
library(batchmeans); library(foreach)
library(xtable)

bmmean = function(x) { format(round(bm(x)$est, 1), nsmall = 1) }
lb95 = function(x){ HPDinterval(as.mcmc(x), prob = 0.95)[1] }
ub95 = function(x){ HPDinterval(as.mcmc(x), prob = 0.95)[2] }
lb90 = function(x){ HPDinterval(as.mcmc(x), prob = 0.90)[1] }
ub90 = function(x){ HPDinterval(as.mcmc(x), prob = 0.90)[2] }


path.data = paste0('data/')
path.fit = paste0('fit/')
path.fig = paste0('fig/')
path.sum = paste0('sum/')

ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)


datasets = c('NHPP', 'LGCP', 'NHPPSE', 'LGCPSE')

runID = 1
datai = paste0('sim', datasets[runID])
fiti = 'LGCPSE'


path.r = paste0('src/RFtns.R')
path.cpp = paste0('src/RcppFtns.cpp')


# burn = 10000 # after seeing the trace plots of alphas
burn = 100


# =============================================================================-
# load ----
# =============================================================================-
load(paste0(path.data, 'ccb_env.RData'))
load(paste0(path.data, datai, '.RData'))
load(paste0(path.fit, datai, '_', fiti, '.RData'))


ts = data$ts[,1] # unit is minutes
marks = data$ts[,2]-1
soundspeed = 1.5 * 60 # km per min
sback = 20
knts = unique(c(0, seq(0, maxT, by = sback), maxT))
m = length(knts) - 1
# rho_beta = max(distmat) / 3 # effective range is the maximum distance between HPs
# rho_w = 60 # effective range is rho_w * 3

UTC = data.frame(ts = knts) %>%
  left_join(UTC)

K = nrow(distmat)
p = ncol(Xm[[1]])



# =============================================================================-
# harmonic effects ----
# =============================================================================-

## compute ----
postBetanew = foreach(k = 1:K) %do% {
  sapply(1:p, function(j) postBeta[[j]][-(1:burn),k])
}

harmInd = c(2:p)
data.xb = c()
for(k in 1:K){
  XBharm = postBetanew[[k]][,harmInd] %*% t(Xm[[k]][,harmInd])
  XBharmci = t(sapply(1:ncol(XBharm), function(ii) HPDinterval(as.mcmc(XBharm[,ii]))[1:2]))

  data.xb = rbind(data.xb, data.frame(fit = fiti, Name = 'Xk*Bk with harmonics only', HP = k, ts = knts, UTC = UTC$UTC, value = colMeans(XBharm), lb = XBharmci[,1], ub = XBharmci[,2]))
}
save(data.xb, file = paste0(path.sum, 'xb.RData'))


## figure ----
data.xb$HP = factor(data.xb$HP, levels = paste(1:K), labels = paste0('MARU ', 1:K))

library(suncalc)
dummy = getSunlightTimes(date = seq.Date(as.Date('2010-04-02'), as.Date('2010-04-03'), by = 1),
                         lat = 42, lon = -70.4, keep = c("sunrise", "sunset"), tz = "EST")
shade = data.frame(dusk = as.POSIXlt(as.character(dummy$sunset[1]), tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                   dawn = as.POSIXlt(as.character(dummy$sunrise[2]), tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'),
                   top = Inf,
                   bottom = -Inf)

plot.xb = data.xb %>% 
  filter(UTC %in% seq.POSIXt(as.POSIXlt("2010-04-02 12:00:00", tz = 'UTC', format = '%Y-%m-%d %H:%M:%S'), by = 'min', length.out = 24*60+1)) %>% 
  filter(fit == fiti, Name == 'Xk*Bk with harmonics only') %>%
  ggplot() +
  geom_rect(data = shade,
            aes(xmin = dusk, xmax = dawn, ymin = bottom, ymax = top),
            fill = 'gray50', alpha = 0.5) +
  geom_ribbon(aes(x = UTC, ymin = lb, ymax = ub), fill = "lightsteelblue", alpha = 0.6) +
  geom_line(aes(x = UTC, y = value)) +
  facet_wrap(~ HP, nrow = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Hour', y = 'Posterior mean estimate', title = '(b) Diel effect over the 24-hour window') +
  theme_bw() +
  scale_x_datetime(date_breaks = "4 hours", date_labels = "%H") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))
plot.xb




# =============================================================================-
# credible intervals ----
# =============================================================================-

postbeta1s = postBeta[[1]][-(1:burn),]

Predictors = paste0('MARU ', 1:K)

dat = data.frame(
  Mean = colMeans(postbeta1s), 
  lb95 = apply(postbeta1s, 2, lb95), ub95 = apply(postbeta1s, 2, ub95),
  lb90 = apply(postbeta1s, 2, lb90), ub90 = apply(postbeta1s, 2, ub90),
  Predictor = Predictors
)


dat$Predictor = factor(dat$Predictor, levels = rev(unique(dat$Predictor)))
# dat = dat %>% mutate(hasZero = ifelse(CIlow <= 0 & CIup >= 0, 'Include zero', 'Do not include zero'))


plot.ci = dat %>% 
  ggplot(aes(y = Predictor)) +
  # geom_errorbar(aes(xmin = lb95, xmax = ub95), width = 0.2) +
  geom_point(aes(x = Mean), size = 2) +
  geom_linerange(aes(xmin = lb95, xmax = ub95)) +
  geom_linerange(aes(xmin = lb90, xmax = ub90), linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.2) +
  # scale_x_break(c(-0.25, -0.13)) +
  # scale_x_break(c(0.21, 1.77), ticklabels = c(1.8)) +
  labs(x = 'HPD interval', y = '', title = '(a) Ambient noise') +
  guides(linetype = guide_legend(title=""), shape = guide_legend(title="")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(vjust = 0.0001)
  )
plot.ci




# =============================================================================-
# combine figures ----
# =============================================================================-

plot.both = grid.arrange(plot.ci, plot.xb, nrow = 1, widths = c(1, 2.2))

# ggsave(plot = plot.both, width = 9.6, height = 3.1,
#        filename = paste0(path.fig, datai, fiti, '_CInXB.pdf'))
# 














