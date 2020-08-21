MeanShift_null <- function(data_in, param){
  
  dataDir = file.path(createPaths(), 'RData')
  
  data_MS <- data_in
  
  for (i in 1:length(data_MS)) {
   
    
    #print(paste0('RECORD ', i))
    
    data_MS[[i]]$null_sig_brks = list()
    data_MS[[i]]$null_brk_dirs = list()
    
    synthDat = try(createSyntheticTimeseries(time = data_MS[[i]]$age, 
                                    values = data_MS[[i]]$paleoData_values, 
                                             n.ens = param$numIt), silent = T)
    
    if (class(synthDat) == "try-error") {
      print(paste('Try error, running instead with method = ML'))
      synthDat = try(synth_fun(data_MS[[i]]$age, 
                               data_MS[[i]]$paleoData_values, 
                               nens = param$numIt))
      
      if (class(synthDat) == "try-error") {
        print(paste('Try error again, skipping this record'))
        data_MS[[i]]$useMS = -1
        next
      }
    }
    
    registerDoParallel(cores = param$ncores)
    # run the mean shift code for all iterations
    #for (it in 1:numIt) {
    out <-foreach(it=1:param$numIt,
                  .verbose=F,.errorhandling = "pass") %dopar% { 
                    
                    #print(paste0('MS null ITERATION ', it))
                    
                    output = MS_fun(data_MS[[i]]$age, synthDat[,it])
                    
                    return(list(sig = output$sig_brks, brk = output$brk_dirs))
                  }
    stopImplicitCluster()
    
      #data_MS[[i]]$null_sig_brks[[it]] = output$sig_brks
      #data_MS[[i]]$null_brk_dirs[[it]] = output$brk_dirs
    data_MS[[i]]$null_sig_brks = lapply(out, `[[`, 1)
    data_MS[[i]]$null_brk_dirs = lapply(out, `[[`, 2)
    
  }
  
  save(data_MS, file = file.path(dataDir, 'MS_results_plusNull_complete.RData'))
  return(data_MS)
  
}