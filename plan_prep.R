plan_prep <- function(data_in = data_all, param){
  
  # Data preparation
  data_index  = inThisCompilation(TS = data_in,compName = param$CName, compVers = param$CVers)
  data_comp = data_all[unlist(data_index)]
  climate_index = climate_indices(data_comp)
  data = processLiPD(data_in = data_comp, climInds = climate_index, 
                     detrend = param$det)
  out = {
    datPath = file.path(createPaths(), 'RData', param$OutDat)
    save(data, file = datPath)
  }
  
  
  return(data)
  
}
