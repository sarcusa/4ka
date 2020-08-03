TrendChanges <- function(input_data, input_param, input_var){
  
  analysis_3c = spatialBrokenStick(input_data, input_param, 
                                   climateVar = input_var)
  
    
  analysis_3d = hist_BS(data_in = input_data, param = input_param,
                        climateVar = input_var)
  
  out3d = {
    datPath = file.path(createPaths(), 'RData', 'TrendChange_histogram.RData')
    save(analysis_3d, file = datPath)
  }
  
  analysis_3e = hist_BS_plot(data_in = analysis_3d, param = input_param,
                             climateVar = input_var)
    
  #output description
  # analysis_3e[[1]][[1]] = line plot
  # analysis_3e[[1]][[2]] = net histogram
  output  <- list(analysis_3e[[1]][[1]], analysis_3e[[1]][[2]], analysis_3d)
  return(output)
}