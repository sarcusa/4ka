TrendChanges <- function(input_data, input_param, input_var){
  
  analysis_3a = BrokenStick(input_data, input_param)
  
  out3a = {
    datPath = file.path(createPaths(), 'RData', 'TrendChanges.RData')
    save(analysis_3a, file = datPath)
  }
  
  #cl <- parallel::makeCluster(2L)
  #print("calling workers")
  #future::plan(cluster, workers = cl)
  #analysis_3b = future({BrokenStick_null(analysis_3a[[1]], input_param)},
                       lazy = TRUE)
  #future::plan(sequential)
  #print("going back to sequential")
  


  analysis_3c = spatialBrokenStick(analysis_3b[[1]], input_param, 
                                   climateVar = input_var)
  
    
  analysis_3d = hist_BS(data_in = analysis_3b[[1]], param = input_param,
                        climateVar = input_var)
  
  out_3b = {
    datPath = file.path(createPaths(), 'RData', 'BS_results_plusNull_complete.RData')
    save(analysis_3b, file = datPath)
  }
  
  out3d = {
    datPath = file.path(createPaths(), 'RData', 'TrendChange_histogram.RData')
    save(analysis_3d, file = datPath)
  }
  
  analysis_3e = hist_BS_plot(data_in = analysis_3d, param = input_param,
                             climateVar = input_var)
    
  #output description
  # analysis_3a[[2]] = histogram
  # analysis_3b[[2]] = histogram
  # analysis_3e[[1]][[1]] = line plot
  # analysis_3e[[1]][[2]] = net histogram
  output  <- list(analysis_3a[[2]], analysis_3b[[2]],
                  analysis_3e[[1]][[1]], analysis_3e[[1]][[2]], analysis_3b[[1]], analysis_3d)
  return(output)
}