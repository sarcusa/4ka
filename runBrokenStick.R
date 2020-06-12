library(geoChronR)
library(segmented)
library(lipdR)

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('gaussianize.R')
source('brokenStick.R')
source('createPaths.R')

mainDir = createPaths()
load(file.path(mainDir, 'RData', 'TS_climateInterp_2019.RData'))
#load('RData/TS_climateInterp_complete.RData')

TS_BS = filterTs(TS, 'useBS == 1')

#figDir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/BS_figs/'
figDir = file.path(mainDir, 'broken_stick', 'individual')
dir.create(figDir)

for (i in 1:length(TS_BS)) {
  
  print(paste0('RECORD ', i))

  # run the broken stick code
  results = iterativeBrokenStick(TS_BS[[i]]$age, TS_BS[[i]]$paleoData_values, 
                                 plot.opt = T, plotName = TS_BS[[i]]$dataSetName, figDir = figDir)
  
  if (is.na(results$cp)) {
    
    TS_BS[[i]]$brk_pts = NA
    TS_BS[[i]]$brk_ptsErr = NA
    TS_BS[[i]]$brk_dirs = NA
    
  } else {
    
    ordI = order(results$cp)
    TS_BS[[i]]$brk_pts = results$cp[ordI]
    TS_BS[[i]]$brk_ptsErr = results$cp.se[ordI]
    slopes = results$o$coefficients[3:(length(results$cp)+2)]
    TS_BS[[i]]$brk_dirs = ifelse(slopes > 0, 1, -1)
  }
  
  
} # end loop thru records

fileName = file.path(mainDir, 'RData', 'BS_results_complete.RData')
save(TS_BS, file = fileName)
#save(TS_BS, file = 'BS_results_complete.RData')

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


