processLiPD <- function(data_in, climInds){
  
  for (i in 1:length(data_in)) {
    
    # set the climate interpretation flag
    if (any(climInds == i)) {
      if (!is.na(data_in[[i]]$interpretation1_variable) & data_in[[i]]$interpretation1_variable != 'NA') {
        data_in[[i]]$climateInterpretation = TRUE
      } else {
        data_in[[i]]$climateInterpretation = FALSE
        next
      }
    } else {
      data_in[[i]]$climateInterpretation = FALSE
      next
    }
    
    age = data_in[[i]]$age
    vals = data_in[[i]]$paleoData_values
    
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
      print(paste('OUT OF ORDER AGES: INDEX ', i, 'SITE: ', data_in[[i]]$dataSetName))
      ind_order = order(age)
      age = age[ind_order]
      vals = vals[ind_order]
    }
    
    # average overlapping ages if necessary
    if (any(diff(age) == 0)) {
      print(paste('OVERLAPPING AGES: INDEX ', i, 'SITE: ', data_in[[i]]$dataSetName))
      inds = which(diff(age) == 0)
      vals[inds] = (vals[inds] + vals[inds+1])/2
      age = age[-(inds+1)]
      vals = vals[-(inds+1)]
    }
    
    # restrict data_in to Holocene
    inds = which(age <= 11500)
    age = age[inds]
    vals = vals[inds]
    
    # Update TS with processed age and paleodata_in vectors
    data_in[[i]]$age = age
    data_in[[i]]$paleoData_values = vals
    
    # Initialize flags for analysis - assume all records are being used, then filter
    # 1 = used in analysis, 0 = data_inset doesn't meet criteria for analysis
    data_in[[i]]$useEX = 1   # use in excursion
    data_in[[i]]$useMS = 1   # use in mean shift
    data_in[[i]]$useBS = 1   # use in broken stick
    data_in[[i]]$useBC = 1   # use in Bayesian change point
    
    ## ------------------ FILTER RECORDS: ASSIGN FLAGS --------------------- ##
    if (length(vals) == 0 || length(vals) == 1) {
      print(paste('(essentially) EMPTY VECTOR: INDEX ', i, 'SITE: ', data_in[[i]]$dataSetName))
      data_in[[i]]$useEX = 0
      data_in[[i]]$useMS = 0
      data_in[[i]]$useBS = 0
      data_in[[i]]$useBC = 0
      next()
    }
    
    if (length(age) != length(vals)) {
      print(paste('DIFFERENT VECTOR LENGTHS: INDEX ', i, 'SITE: ', data_in[[i]]$dataSetName))
      data_in[[i]]$useEX = 0
      data_in[[i]]$useMS = 0
      data_in[[i]]$useBS = 0
      data_in[[i]]$useBC = 0
      next()
    }
    
    # assume the data_in is quantized if over half the points are the same as their neighbor
    if (sum(diff(vals) == 0, na.rm = T) > length(vals) / 2) {
      print(paste('STEP FUNCTION/QUANTIZED VALUES: INDEX ', i, 'SITE: ', data_in[[i]]$dataSetName))
      data_in[[i]]$useEX = 0
      data_in[[i]]$useMS = 0
      data_in[[i]]$useBS = 0
      data_in[[i]]$useBC = 0
      next()
    }
    
    if (all.equal(diff(age), rep(50, length(age) - 1)) == TRUE) {
      paste('50-YR INTERPOLATED data_in SET: INDEX ', i, 'SITE: ', data_in[[i]]$dataSetName)
      data_in[[i]]$useEX = 0
    }
    
    if (max(age) - min(age) < 4000) {
      print(paste('RECORD TOO SHORT FOR BS, MS, BC: INDEX ', i))
      data_in[[i]]$useBS = 0
      data_in[[i]]$useMS = 0
      data_in[[i]]$useBC = 0
    }
    
    if (median(diff(age)) > 500) {
      print(paste('RECORD TOO COARSE FOR BS, MS, BC: INDEX ', i))
      data_in[[i]]$useBS = 0
      data_in[[i]]$useMS = 0
      data_in[[i]]$useBC = 0
    }
    
    # correct discrepancies in interpretation
    if (length(data_in[[i]]$interpretation1_interpDirection) < 1) {
      
      # these existing LiPD files don't have assigned interpretations (can't currently edit, should be +)
      names = c('MV99-GC41','Core17940','MD94-103.Sicre.2005','Nujulla.Larocque.2004','TanaLake')
      
      if (length(data_in[[i]]$interpretation1_direction) > 0) {
        data_in[[i]]$interpretation1_interpDirection = data_in[[i]]$interpretation1_direction
      } else if (any(data_in[[i]]$dataSetName == names)) {
        data_in[[i]]$interpretation1_interpDirection = 'positive'
      } else {
        #print(paste0('No interp: ', TS[[i]]$data_inSetName, ', index = ', i))
        data_in[[i]]$interpretation1_interpDirection = 'NA'
      }
    }
    
  }
  
  data <- filterTs(data_in, 'climateInterpretation == TRUE')
  return(data)
  
}  
  