Excursion_null <- function(data_in, param){
  
  dataDir = file.path(createPaths(), 'RData')
  
  for (y in 1:length(param$eventYrs)) {
    
    remove(data_EX)
    event_yr = param$eventYrs[y]
    print(paste('Event year:', event_yr))
    
    # Load in the TS structure of excursion results, plus event_window and ref_window
    list2env(loading(file.path(dataDir, paste0('EX_results', event_yr/1000, '.RData'))),envir=.GlobalEnv)
    #load(file.path(dataDir, paste0('EX_results', event_yr/1000, '.RData')))
    data_EX = filterTs(data_EX, 'statusEX == 1')
    
    if(!length(data_EX) == 0){
    
    for (i in 1:length(data_EX)) {
            
      #print(paste0('RECORD ', i))
      
      data_EX[[i]]$null_events = list()
      
      ageInds = which(data_EX[[i]]$age <= 7500)
      if (length(ageInds) == 0) {
        print(paste('Index', i, ': No ages before 7.5 ka'))
        ageInds = NA
      }
      
      synthDat = try(createSyntheticTimeseries(time = data_EX[[i]]$age, 
                                      values = data_EX[[i]]$paleoData_values,
                                      n.ens = param$numIt, 
                                    index.to.model = ageInds),silent = T)
      
      if (class(synthDat) == "try-error") {
        print(paste('Try error, running instead with method = ML'))
        synthDat = try(synth_fun(data_EX[[i]]$age, 
                                 data_EX[[i]]$paleoData_values, 
                                 nens = param$numIt,index.to.model = ageInds))
        
        if (class(synthDat) == "try-error") {
          print(paste('Try error again, skipping this record'))
          data_EX[[i]]$useEX = -999
          next
        }
      }
      #registerDoParallel(cores = param$ncores)
      
      registerDoParallel(cores = Sys.getenv("SLURM_CPUS_PER_TASK"))
      
      #for (it in 1:param$numIt) {
      out <-foreach(it=1:param$numIt,
                    .verbose=F,.errorhandling = "pass") %dopar% { 
                      
                      #print(paste('Iteration', it))
                      
                                         
                      results = EX_fun(data_EX[[i]]$age, synthDat[,it], 
                                       event_yr = event_yr, 
                                       event_window = param$event_window, 
                                       ref_window = param$ref_window)
                      
                      return(list(null_2 = results[[2]], null_3 = results[[3]]))
                    }
      stopImplicitCluster()
        
      data_EX[[i]]$null_events = lapply(out, `[[`, 1)
      data_EX[[i]]$null_events_dir = unlist(lapply(out, `[[`, 2))
      
      #data_EX[[i]]$null_events[[it]] = results[[2]]
      #data_EX[[i]]$null_events_dir[[it]] = results[[3]]
       
      
    } # end loop thru records
    event_window = param$event_window
    save(data_EX, event_window, file = file.path(dataDir, paste0('EX_results_plusNull_complete_', event_yr/1000, '.RData')))
    
    }else{
      print("no records")
      event_window = param$event_window
      save(event_window, file = file.path(dataDir, paste0('EX_results_plusNull_complete_', event_yr/1000, '.RData')))
      
    }
  }
  
  return(data_EX)

}