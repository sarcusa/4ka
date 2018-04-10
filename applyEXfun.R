# Hannah Kolus, hrk37@nau.edu
# 
# Determines whether an excursion event has occurred within the specified event window.
# Excursion events are defined as two consecutive values within the event window that
# are more extreme than the avg +/- X std of the reference windows.

library(lipdR)
library(tidyverse)
library(ggplot2)
source('EX_function.R')

load('TS_climateInterp.RData') # Load in the TS structure of all paleoclimate data
TS_EX = filterTs(TS, 'useEX == 1')

## --------- SET PARAMETERS DEFINING THE ANALYSIS WINDOW AND EVENTS --------- ##
event_yr = 4200     # Set the event year
event_window = 600  # (event_yr +/- event_window/2) as the event window
ref_window = 500    # Set the duration of the reference windows 

# NOTE: these are the default values
resCriteria = 50    # Resolution criteria (yr) required to include record
sigNum = 2          # Set how many sigmas below and above the mean define an excursion

# PLOT OPTIONS: specify whether to plot and directory for figures
plotOpt = T
figDir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka/4ka_figs/'

for (i in 1:length(TS_EX)) {
  
  print(paste('Iteration', i))
  
  age = TS_EX[[i]]$age
  vals = TS_EX[[i]]$paleoData_values
  
  results = EX_fun(age, vals, event_yr = event_yr, event_window = event_window, 
         ref_window = ref_window, plotOpt = plotOpt,
         datNam = TS_EX[[i]]$dataSetName, varNam = TS_EX[[i]]$paleoData_variableName,
         units = TS_EX[[i]]$paleoData_units, proxy = TS_EX[[i]]$paleoData_proxy)
  
  TS_EX[[i]]$statusEX = results[[1]]
  TS_EX[[i]]$eventEX = results[[2]]
  
}

fileName = paste0('EX_results', event_yr/1000, '.RData')
save(TS_EX, event_yr, event_window, ref_window, file = fileName)
