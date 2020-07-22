source('createPaths.R')
source('inThisCompilation.R')
source('process_LiPD.R')
source('plot_sites.R')
source('EX_function.R')
source('applyEXfun.R')
source('set_parameters.R')
source('runExcursion_null.R')
source('spatialExcursion_direction_null.R')
source('histogram_EX.R')
source('histogram_EX_plot.R')
source('applyMSfun.R')
source('MS_function.R')
source('runMeanShift_null.R')
source('spatialMeanShift_direction_null.R')
source('histogram_MS.R')
source('histogram_MS_plot.R')
source('createSyntheticTimeseries_v2.R')
source('gaussianize.R')
source('brokenStick.R')
source('runBrokenStick.R')
source('runBrokenStick_null.R')
source('spatialBrokenStick_direction_null.R')
source('histogram_BS.R')
source('histogram_BS_plot.R')
source('histogram_netPlot_doubleAxis.R')
source('plotProxyMap.R')

climate_indices <- function(data){
  climInterp1 = sapply(data,"[[","interpretation1_variable")
  inds1 = which(sapply(climInterp1, function (x) length(x) > 0))
  
  climInterp2 = sapply(data,"[[","interpretation2_variable")
  inds2 = which(sapply(climInterp2, function (x) length(x) > 0))
  
  climInds = union(inds1, inds2)
  
  return(climInds)
  
}

loading <- function(dataset){
  file <- load(dataset)
  return(mget(file))
}
