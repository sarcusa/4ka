MeanShift <- function(input_data,input_param, input_var){
  
  analysis_2c = spatialMeanShift_null(data_in = input_data, 
                                      param = input_param, climateVar = input_var)
  
  analysis_2d = histogram_MS(data_in = input_data, 
                             param = input_param, climateVar = input_var)
  
  out2d = {
    datPath = file.path(createPaths(), 'RData', 'MeanShift_histogram.RData')
    save(analysis_2d, file = datPath)
  }
  
  
  analysis_2e = hist_MS_plot(data_in = analysis_2d, 
                             param = input_param, climateVar = input_var)
  
  
  output = list(analysis_2e, analysis_2d)
  
  return(output)
  
}