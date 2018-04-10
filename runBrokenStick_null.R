## Written by Hannah Kolus, 09/04/2018 
## Runs broken stick analysis using synthetic data across all relevant (flagged) proxy records

library(geoChronR)
library(segmented)
library(lipdR)

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('gaussianize.R')
source('brokenStick.R')
load('BS_results.RData')

numIt = 1000

for (i in 1:length(TS_BS)) {
  
  print(paste0('RECORD ', i))
  
  TS_BS[[i]]$null_brk_pts = list()
  TS_BS[[i]]$null_brk_ptsErr = list()
  
  synthDat = createSyntheticTimeseries(TS_BS[[i]]$age, TS_BS[[i]]$paleoData_values, nens = numIt)

  # run the broken stick code for all iterations
  for (it in 1:numIt) {
    print(paste0('ITERATION ', it))
    
    results = iterativeBrokenStick(TS_BS[[i]]$age, synthDat[,it], plot.opt = F)
    TS_BS[[i]]$null_brk_pts[[it]] = results$cp
    TS_BS[[i]]$null_brk_ptsErr[[it]] = results$cp.se
  }
  
} # end loop thru records

save(TS_BS, file = 'BS_results_plusNull_1_9.RData')

## ---------------------- HISTOGRAM PLOT ---------------------- ##

histBreaks = seq(100, 11500, by = 200)
yrs = seq(300, 11500, by = 600)
yrStr = yrs / 1000

p=ggplot() + geom_histogram(aes(breakPts), breaks = yrs)
pg = ggplot_build(p)
densityBS = pg$data[[1]]$density
xvals = pg$data[[1]]$x

ggplot() + geom_col(aes(x = xvals, y = densityBS)) + 
  xlab('Event years [BP]') + ylab('Density of events') + 
  ggtitle('Broken Stick Results (simple plot)')


