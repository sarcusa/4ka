MeanShift_null <- function(data_in, param){
  
  dataDir = file.path(createPaths(), 'RData')
  
  data_MS <- data_in
  
  for (i in 1:length(data_MS)) {
    
    print(paste0('RECORD ', i))
    
    data_MS[[i]]$null_sig_brks = list()
    data_MS[[i]]$null_brk_dirs = list()
    
    synthDat = try(createSyntheticTimeseries(data_MS[[i]]$age, 
                                             data_MS[[i]]$paleoData_values, 
                                             nens = numIt), silent = T)
    
    if (class(synthDat) == "try-error") {
      print(paste('Try error, running instead with method = ML'))
      synthDat = try(synth_fun(data_MS[[i]]$age, 
                               data_MS[[i]]$paleoData_values, 
                               nens = numIt))
      
      if (class(synthDat) == "try-error") {
        print(paste('Try error again, skipping this record'))
        data_MS[[i]]$useMS = -1
        next
      }
    }
    
    # run the mean shift code for all iterations
    for (it in 1:numIt) {
      print(paste0('ITERATION ', it))
      
      output = MS_fun(data_MS[[i]]$age, synthDat[,it])
      data_MS[[i]]$null_sig_brks[[it]] = output$sig_brks
      data_MS[[i]]$null_brk_dirs[[it]] = output$brk_dirs
    }
    
  }
  
  save(data_MS, file = file.path(dataDir, 'MS_results_plusNull_complete.RData'))
  return(data_MS)
  
}