## Written by Hannah Kolus, 09/04/2018 
## Runs excursion analysis using synthetic data across all relevant (flagged) proxy records

library(lipdR)
library(tidyverse)
library(ggplot2)
library(geoChronR)
library(segmented)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('EX_function.R')

## --------- SET PARAMETERS DEFINING THE ANALYSIS WINDOW AND EVENTS --------- ##
event_yr = 4200     # Set the event year
numIt = 1000

# Load in the TS structure of excursion results, plus event_window and ref_window
load(paste0('EX_results', event_yr/1000, '.RData')) 
TS_EX = filterTs(TS_EX, 'statusEX == 1')

for (i in 1:length(TS_EX)) {
  
  print(paste0('RECORD ', i))
  
  TS_EX[[i]]$null_events = list()
  
  ageInds = which(TS_EX[[i]]$age <= 7500)
  if (length(ageInds) == 0) {
    print(paste('Index', i, ': No ages before 7.5 ka'))
    ageInds = NA
  }
  
  synthDat = createSyntheticTimeseries(TS_EX[[i]]$age, TS_EX[[i]]$paleoData_values, 
                                       nens = numIt, index.to.model = ageInds)
  
  for (it in 1:numIt) {
    print(paste('Iteration', it))
    
    results = EX_fun(TS_EX[[i]]$age, synthDat[,it], event_yr = event_yr, 
                     event_window = event_window, ref_window = ref_window)
    
    TS_EX[[i]]$null_events[[it]] = results[[2]]
  } 
  
} # end loop thru records

save(TS_EX, event_window, file = paste0('EX_results_plusNull2_', event_yr/1000, '.RData'))
