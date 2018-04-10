## Written by Hannah Kolus, 09/04/2018 
# Process all LiPD files, creating a TS structure
# Remove all missing values from ages and paleodata, reorder by age if necessary
# Filter by those that have a climate interpretation

library(lipdR)
library(tidyverse)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
# dir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/R/LiPD_ex'
dir = '/Users/hannah/Dropbox/LiPDLibrary'

D = readLipd(dir)
TS = extractTs(D)

# queryTs(TS, 'interpretation1_variable == M')
climInterp1 = sapply(TS,"[[","interpretation1_variable")
inds1 = which(sapply(climInterp1, function (x) length(x) > 0))

climInterp2 = sapply(TS,"[[","interpretation2_variable")
inds2 = which(sapply(climInterp2, function (x) length(x) > 0))

climInds = union(inds1, inds2)

for (i in 1:length(TS)) {
  
  if (any(climInds == i)) {
    if (!is.na(TS[[i]]$interpretation1_variable) & TS[[i]]$interpretation1_variable != 'NA') {
      TS[[i]]$climateInterpretation = TRUE
    } else {
      TS[[i]]$climateInterpretation = FALSE
    }
  } else {
    TS[[i]]$climateInterpretation = FALSE
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
    next()
  }
  
  if (max(age) - min(age) < 4000) {
    print(paste('RECORD TOO SHORT FOR BS, MS, BC: INDEX ', i))
    TS[[i]]$useBS = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBC = 0
    next()
  }
  
  if (median(diff(age)) > 500) {
    print(paste('RECORD TOO COARSE FOR BS, MS, BC: INDEX ', i))
    TS[[i]]$useBS = 0
    TS[[i]]$useMS = 0
    TS[[i]]$useBC = 0
    next()
  }

} # end loop thru records

newTS = filterTs(TS, 'climateInterpretation == TRUE')
TS = newTS

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
save(TS, file = 'TS_climateInterp.RData')
