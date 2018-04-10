## Written by Hannah Kolus, 09/04/2018 
## Runs the mean shift function across all relevant (flagged) proxy records

library(changepoint)
library(lipdR)
source('MS_function.R')

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')

# PLOT OPTIONS: specify whether to plot and directory for figures
figDir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka/MS_figs/'

load('TS_climateInterp.RData') # Load in the TS structure of all paleoclimate data
TS_MS = filterTs(TS, 'useMS == 1')

for (i in 1:length(TS_MS)) {
  
  print(paste('Iteration', i))
  
  age = TS_MS[[i]]$age
  vals = TS_MS[[i]]$paleoData_values
  
  TS_MS[[i]]$sig_brks = MS_fun(age, vals, datNam = TS_MS[[i]]$dataSetName, varNam = TS_MS[[i]]$paleoData_variableName)
  
}

save(TS_MS, file = 'MS_results.RData')