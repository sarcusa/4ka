Excursion <- function(input_data,input_param, input_var){
  
  analysis_1a = applyEX(data_in = input_data, param = input_param)
  
  out1a = {
    datPath = file.path(createPaths(), 'RData', 'Excursion.RData')
    save(analysis_1a, file = datPath)
  }
  
  analysis_1b = Excursion_null(data_in = analysis_1a, param = input_param)
  
  out1b = {
    datPath = file.path(createPaths(), 'RData', 'Excursion_plusnull.RData')
    save(analysis_1b, file = datPath)
  }
  
  analysis_1c = spatialExcursion_null(data_in = analysis_1b, 
                                      param = input_param, climateVar = input_var)
  
  analysis_1d = histEX(data_in = analysis_1b, param = input_param, 
                       climateVar = input_var)
  
  
  analysis_1e = histEX_plot(data_in = analysis_1d, param = input_param, 
                            climateVar = input_var)
  
  out1e = {
    datPath = file.path(createPaths(), 'RData', 'Excursion_histogram.RData')
    save(analysis_1e, file = datPath)
  }
  
  output = list(analysis_1e, analysis_1b, analysis_1d)  
  
  return(output)
  
}

