## Written by Hannah Kolus, 09/04/2018 
## Runs mean shift analysis using synthetic data across all relevant (flagged) proxy records

library(changepoint)
library(lipdR)
library(geoChronR)
library(segmented)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('MS_function.R')
source('createPaths.R')
source('createSyntheticTimeseries_v2.R') # method = 'ML' for cases where arima fails

# SET PARAMETERS HERE - (Note: these are default, won't include them as inputs)
maxDiff = 1000    # min segment length between change points
alpha = 0.05      # set confidence level for mean difference testing
gaussOpt = F      # option to gaussianize values before analysis

dataDir = file.path(createPaths(), 'RData')
load(file.path(dataDir, 'MS_results_complete.RData'))

numIt = 1000

for (i in 1:length(TS_MS)) {
  
  print(paste0('RECORD ', i))
  
  TS_MS[[i]]$null_sig_brks = list()
  TS_MS[[i]]$null_brk_dirs = list()
  
  synthDat = try(createSyntheticTimeseries(TS_MS[[i]]$age, TS_MS[[i]]$paleoData_values, nens = numIt), silent = T)
  
  if (class(synthDat) == "try-error") {
    print(paste('Try error, running instead with method = ML'))
    synthDat = try(synth_fun(TS_MS[[i]]$age, TS_MS[[i]]$paleoData_values, nens = numIt))
    
    if (class(synthDat) == "try-error") {
      print(paste('Try error again, skipping this record'))
      TS_MS[[i]]$useMS = -1
      next
    }
  }
  
  # run the mean shift code for all iterations
  for (it in 1:numIt) {
    print(paste0('ITERATION ', it))
    
    output = MS_fun(TS_MS[[i]]$age, synthDat[,it])
    TS_MS[[i]]$null_sig_brks[[it]] = output$sig_brks
    TS_MS[[i]]$null_brk_dirs[[it]] = output$brk_dirs
  }
  
} # end loop thru records

save(TS_MS, file = file.path(dataDir, 'MS_results_plusNull_complete.RData'))
