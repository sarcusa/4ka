applyMS <- function(data_in, param){
  
  figDir = file.path(createPaths(), 'mean_shift', 'individual')
  dir.create(figDir)
  
  data_MS = filterTs(data_in, 'useMS == 1')
  
  for (i in 1:length(data_MS)) {
    
    print(paste('Record', i))
    
    age = data_MS[[i]]$age
    vals = data_MS[[i]]$paleoData_values
    
    output = MS_fun(age, vals, plotOpt = T, figDir = figDir, 
                    maxDiff = param$maxDiff,
                    datNam = data_MS[[i]]$dataSetName, 
                    varNam = data_MS[[i]]$paleoData_variableName)
    
    data_MS[[i]]$sig_brks = output$sig_brks
    data_MS[[i]]$brk_dirs = output$brk_dirs
    
  }
  
  return(data_MS)
  
}