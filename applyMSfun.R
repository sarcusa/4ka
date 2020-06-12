## Written by Hannah Kolus, 09/04/2018 
## Runs the mean shift function across all relevant (flagged) proxy records

library(changepoint)
library(lipdR)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('MS_function.R')
source('createPaths.R')

# PLOT OPTIONS: specify whether to plot and directory for figures
figDir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/MS_figs/'
figDir = file.path(createPaths(), 'mean_shift', 'individual')
dir.create(figDir)

load(file.path(createPaths(), 'RData', 'TS_climateInterp_2019.RData')) # Load in the TS structure of all paleoclimate data
TS_MS = filterTs(TS, 'useMS == 1')

for (i in 1:length(TS_MS)) {
  
  print(paste('Record', i))
  
  age = TS_MS[[i]]$age
  vals = TS_MS[[i]]$paleoData_values
  
  output = MS_fun(age, vals, plotOpt = T, figDir = figDir,
                               datNam = TS_MS[[i]]$dataSetName, varNam = TS_MS[[i]]$paleoData_variableName)
  
  TS_MS[[i]]$sig_brks = output$sig_brks
  TS_MS[[i]]$brk_dirs = output$brk_dirs
  
}

save(TS_MS, file = file.path(createPaths(), 'RData', 'MS_results_complete.RData'))
