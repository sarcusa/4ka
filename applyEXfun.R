applyEX <- function(data_in, param){
  for (y in 1:length(param$eventYrs)) {
    
    data_EX = filterTs(data_in, 'useEX == 1')
    event_yr = param$eventYrs[y]
    
    figDir = file.path(mainDir, 'excursion', 
                       paste0(event_yr/1000, 'ka_individual'))
    dir.create(figDir)
    print(paste('Event year:', event_yr))
    
    for (i in 1:length(data_EX)) {
      
      print(paste('Record', i))
      
      age = data_EX[[i]]$age
      vals = data_EX[[i]]$paleoData_values
      
      results = EX_fun(age, vals, 
                       event_yr = event_yr, event_window = param$event_window, 
                       ref_window = param$ref_window, 
                       plotOpt = param$plotOpt, figDir = figDir,
                       datNam = data_EX[[i]]$dataSetName, 
                       varNam = data_EX[[i]]$paleoData_variableName,
                       lat = data_EX[[i]]$geo_latitude, 
                       lon = data_EX[[i]]$geo_longitude,
                       units = data_EX[[i]]$paleoData_units, 
                       proxy = data_EX[[i]]$paleoData_proxy,
                       resCriteria = param$resCriteria)
      
      data_EX[[i]]$statusEX = results[[1]]
      data_EX[[i]]$eventEX = results[[2]]
      data_EX[[i]]$dirEx = results[[3]]
      data_EX[[i]]$onsetEX = results[[4]]
      
    }
    
    fileName = file.path(mainDir, 'RData', paste0('EX_results', event_yr/1000, '.RData'))
    event_window = param$event_window
    ref_window = param$ref_window
    save(data_EX, event_yr, event_window, ref_window, file = fileName)
    
  }
  
  return(data_EX)
}