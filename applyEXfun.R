# Hannah Kolus, hrk37@nau.edu
# 
# Determines whether an excursion event has occurred within the specified event window.
# Excursion events are defined as two consecutive values within the event window that
# are more extreme than the avg +/- X std of the reference windows.

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
library(lipdR)
library(tidyverse)
library(ggplot2)
library(ggmap)
library(geosphere)
library(gridExtra)
source('EX_function.R')
source('createPaths.R')

## --------- SET PARAMETERS DEFINING THE ANALYSIS WINDOW AND EVENTS --------- ##
eventYrs = seq(1000,11000,by = 400)
event_yr = 3800 #4050     # Set the event year
event_window = 600 #300  # (event_yr +/- event_window/2) as the event window
ref_window = 500 #260    # Set the duration of the reference windows 

# NOTE: these are the default values
resCriteria = 50 #10    # Resolution criteria (yr) required to include record
sigNum = 2          # Set how many sigmas below and above the mean define an excursion

# PLOT OPTIONS: specify whether to plot and directory for figures
plotOpt = T
mainDir = createPaths()
figDir = paste0('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/EX_figs_', floor(event_yr/1000), 'ka_stackPlot/')

for (y in 1:length(eventYrs)) {
  
  load(file.path(mainDir, 'RData', 'TS_climateInterp_2019.RData'))
  TS_EX = filterTs(TS, 'useEX == 1')
  event_yr = eventYrs[y]
  
  figDir = file.path(mainDir, 'excursion', paste0(event_yr/1000, 'ka_individual'))
  dir.create(figDir)
  print(paste('Event year:', event_yr))

for (i in 1:length(TS_EX)) {
  
  print(paste('Record', i))
  
  age = TS_EX[[i]]$age
  vals = TS_EX[[i]]$paleoData_values
  
  results = EX_fun(age, vals, event_yr = event_yr, event_window = event_window, 
         ref_window = ref_window, plotOpt = plotOpt, figDir = figDir,
         datNam = TS_EX[[i]]$dataSetName, varNam = TS_EX[[i]]$paleoData_variableName,
         lat = TS_EX[[i]]$geo_latitude, lon = TS_EX[[i]]$geo_longitude,
         units = TS_EX[[i]]$paleoData_units, proxy = TS_EX[[i]]$paleoData_proxy,
         resCriteria = resCriteria)
  
  TS_EX[[i]]$statusEX = results[[1]]
  TS_EX[[i]]$eventEX = results[[2]]
  TS_EX[[i]]$dirEx = results[[3]]
  TS_EX[[i]]$onsetEX = results[[4]]
  
}
fileName = paste0('EX_results_stackPlot', event_yr/1000, '.RData')
#fileName = paste0('EX_results4.2_short.RData')
fileName = file.path(mainDir, 'Rdata', paste0('EX_results', event_yr/1000, '.RData'))
save(TS_EX, event_yr, event_window, ref_window, file = fileName)

}
