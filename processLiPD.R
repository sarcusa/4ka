## Written by Hannah Kolus, 09/04/2018 
# Process all LiPD files, creating a TS structure
# Remove all missing values from ages and paleodata, reorder by age if necessary
# Filter by those that have a climate interpretation

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')

library(lipdR)
library(tidyverse)
source('createPaths.R')

#dir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/SISAL'
dir = '/Users/hannah/Desktop/LiPDLibrary_HK'
dir = '/Users/hannah/Dropbox/LiPDLibrary/masterDatabase'

D = readLipd(dir)
TS = extractTs(D)

# queryTs(TS, 'interpretation1_variable == M')
climInterp1 = sapply(TS,"[[","interpretation1_variable")
inds1 = which(sapply(climInterp1, function (x) length(x) > 0))

climInterp2 = sapply(TS,"[[","interpretation2_variable")
inds2 = which(sapply(climInterp2, function (x) length(x) > 0))

climInds = union(inds1, inds2)

for (i in 1:length(TS)) {
  
  # set the climate interpretation flag
  if (any(climInds == i)) {
    if (!is.na(TS[[i]]$interpretation1_variable) & TS[[i]]$interpretation1_variable != 'NA') {
      TS[[i]]$climateInterpretation = TRUE
    } else {
      TS[[i]]$climateInterpretation = FALSE
      next
    }
  } else {
    TS[[i]]$climateInterpretation = FALSE
    next
  }
  
  age = TS[[i]]$age
  vals = TS[[i]]$paleoData_values
  
  # Setting all missing values to NA
  age[age == "NaN"] = NA
  vals[vals == "NaN"] = NA
  age[!is.finite(age)] = NA
  vals[!is.finite(vals)] = NA
  vals[vals == -9999] = NA
  vals[vals == -999] = NA
  vals[vals == 999999] = NA
  vals[vals == 99999] = NA
  vals[vals == 9999] = NA
  
  # get rid of NAs
  indsToRemove = which(is.na(vals) | is.na(age))
  if (length(indsToRemove) > 0) {
    age = age[-indsToRemove]
    vals = vals[-indsToRemove]
  }
  
  # reorder ages if out of order
  if (any(diff(age) < 0)) {
    print(paste('OUT OF ORDER AGES: INDEX ', i, 'SITE: ', TS[[i]]$dataSetName))
    ind_order = order(age)
    age = age[ind_order]
    vals = vals[ind_order]
  }
  
  # average overlapping ages if necessary
  if (any(diff(age) == 0)) {
    print(paste('OVERLAPPING AGES: INDEX ', i, 'SITE: ', TS[[i]]$dataSetName))
    inds = which(diff(age) == 0)
    vals[inds] = (vals[inds] + vals[inds+1])/2
    age = age[-(inds+1)]
    vals = vals[-(inds+1)]
  }
  
  # restrict data to Holocene
  inds = which(age <= 11500)
  age = age[inds]
  vals = vals[inds]
  
  # Update TS with processed age and paleodata vectors
  TS[[i]]$age = age
  TS[[i]]$paleoData_values = vals
  
  # Initialize flags for analysis - assume all records are being used, then filter
  # 1 = used in analysis, 0 = dataset doesn't meet criteria for analysis
  TS[[i]]$useEX = 1   # use in excursion
  TS[[i]]$useMS = 1   # use in mean shift
  TS[[i]]$useBS = 1   # use in broken stick
  TS[[i]]$useBC = 1   # use in Bayesian change point
  
  ## ------------------ FILTER RECORDS: ASSIGN FLAGS --------------------- ##
  if (length(vals) == 0 || length(vals) == 1) {
    print(paste('(essentially) EMPTY VECTOR: INDEX ', i, 'SITE: ', TS[[i]]$dataSetName))
    TS[[i]]$useEX = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBS = 0
    TS[[i]]$useBC = 0
    next()
  }
  
  if (length(age) != length(vals)) {
    print(paste('DIFFERENT VECTOR LENGTHS: INDEX ', i, 'SITE: ', TS[[i]]$dataSetName))
    TS[[i]]$useEX = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBS = 0
    TS[[i]]$useBC = 0
    next()
  }
  
  # assume the data is quantized if over half the points are the same as their neighbor
  if (sum(diff(vals) == 0, na.rm = T) > length(vals) / 2) {
    print(paste('STEP FUNCTION/QUANTIZED VALUES: INDEX ', i, 'SITE: ', TS[[i]]$dataSetName))
    TS[[i]]$useEX = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBS = 0
    TS[[i]]$useBC = 0
    next()
  }
  
  if (all.equal(diff(age), rep(50, length(age) - 1)) == TRUE) {
    paste('50-YR INTERPOLATED DATA SET: INDEX ', i, 'SITE: ', TS[[i]]$dataSetName)
    TS[[i]]$useEX = 0
  }
  
  if (max(age) - min(age) < 4000) {
    print(paste('RECORD TOO SHORT FOR BS, MS, BC: INDEX ', i))
    TS[[i]]$useBS = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBC = 0
  }
  
  if (median(diff(age)) > 500) {
    print(paste('RECORD TOO COARSE FOR BS, MS, BC: INDEX ', i))
    TS[[i]]$useBS = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBC = 0
  }
  
  # correct discrepancies in interpretation
  if (length(TS[[i]]$interpretation1_interpDirection) < 1) {
    
    # these existing LiPD files don't have assigned interpretations (can't currently edit, should be +)
    names = c('MV99-GC41','Core17940','MD94-103.Sicre.2005','Nujulla.Larocque.2004','TanaLake')
    
    if (length(TS[[i]]$interpretation1_direction) > 0) {
      TS[[i]]$interpretation1_interpDirection = TS[[i]]$interpretation1_direction
    } else if (any(TS[[i]]$dataSetName == names)) {
      TS[[i]]$interpretation1_interpDirection = 'positive'
    } else {
      #print(paste0('No interp: ', TS[[i]]$dataSetName, ', index = ', i))
      TS[[i]]$interpretation1_interpDirection = 'NA'
    }
  }
  
  ## JUST FOR SISAL
  # if (length(TS[[i]]$age[TS[[i]]$age >= 3400 & TS[[i]]$age <= 5000]) < 11) {
  #   print(paste('NO PTS FOR 4.2 EXCURSION: INDEX ', i))
  #   TS[[i]]$useEX = 0
  # } else if (median(diff(TS[[i]]$age[TS[[i]]$age >= 3400 & TS[[i]]$age <= 5000])) > 50) {
  #   print(paste('RECORD TOO COARSE FOR 4.2 EXCURSION: INDEX ', i))
  #   TS[[i]]$useEX = 0
  # }
  # 
  # if (min(TS[[i]]$age) > 4200 || max(TS[[i]]$age) < 4200) {
  #   TS[[i]]$useMS = 0
  #   TS[[i]]$useBS = 0
  #   print('Out of range')
  # }

} # end loop thru records

newTS = filterTs(TS, 'climateInterpretation == TRUE')
TS = newTS

# need to get back to the code base
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
datPath = file.path(createPaths(), 'RData', 'TS_climateInterp_2019.RData')
save(TS, file = datPath)

## CHECKING SISAL
library(ggplot2)
library(ggmap)
library(gridExtra)

ex = as.numeric(unlist(sapply(TS,"[[","useEX")))
TS_EX = filterTs(TS, 'useEX == 1')
TS_MS = filterTs(TS, 'useMS == 1')
TS_BS = filterTs(TS, 'useBS == 1')

lat_EX = as.numeric(unlist(sapply(TS_EX,"[[","geo_latitude")))
lon_EX = as.numeric(unlist(sapply(TS_EX,"[[","geo_longitude")))

lat_MS = as.numeric(unlist(sapply(TS_MS,"[[","geo_latitude")))
lon_MS = as.numeric(unlist(sapply(TS_MS,"[[","geo_longitude")))

lat_BS = as.numeric(unlist(sapply(TS_BS,"[[","geo_latitude")))
lon_BS = as.numeric(unlist(sapply(TS_BS,"[[","geo_longitude")))

ggplot() + borders("world", colour="black") + 
  geom_point(aes(x = lon_EX, y = lat_EX, color = 'EX'), size = 3, shape = 0, stroke = 1.2) +
  geom_point(aes(x = lon_MS, y = lat_MS, color = 'MS/BS'), size = 3, shape = 4, stroke = 1.2) +
  scale_color_manual(name = '', values = c('EX' = 'blue', 'MS/BS' = 'red'),
                     guide = guide_legend(override.aes = list(shape = c(0,4)))) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle('SISAL sites')
