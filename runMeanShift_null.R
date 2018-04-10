## Written by Hannah Kolus, 09/04/2018 
## Runs mean shift analysis using synthetic data across all relevant (flagged) proxy records

library(changepoint)
library(lipdR)
source('MS_function.R')

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')

# SET PARAMETERS HERE - (Note: these are default, won't include them as inputs)
maxDiff = 1000    # min segment length between change points
alpha = 0.05      # set confidence level for mean difference testing
gaussOpt = F      # option to gaussianize values before analysis

load('MS_results.RData')

numIt = 1000

for (i in 311:length(TS_MS)) {
  
  print(paste0('RECORD ', i))
  
  TS_MS[[i]]$null_sig_brks = list()
  
  synthDat = createSyntheticTimeseries(TS_MS[[i]]$age, TS_MS[[i]]$paleoData_values, nens = numIt)
  
  # run the mean shift code for all iterations
  for (it in 1:numIt) {
    print(paste0('ITERATION ', it))
    
    TS_MS[[i]]$null_sig_brks[[it]] = MS_fun(TS_MS[[i]]$age, synthDat[,it])
  }
  
} # end loop thru records

save(TS_MS, file = 'MS_results_plusNull.RData')
