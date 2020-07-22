plan_prep <- function(data_in = data_all){
  
  # Data preparation
  data_index  = inThisCompilation(TS = data_in,compName = "HoloceneAbruptChange", compVers = "0_9_0")
  data_comp = data_all[unlist(data_index)]
  climate_index = climate_indices(data_comp)
  data = processLiPD(data_in = data_comp, climInds = climate_index)
  out = {
    datPath = file.path(createPaths(), 'RData', 'TS_climateInterp_2020.RData')
    save(data, file = datPath)
  }
  
  
  return(data)
  
}