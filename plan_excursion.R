Excursion <- function(input_data,input_param, input_var){
  
  #analysis_1c = spatialExcursion_null(data_in = input_data, 
  #                                    param = input_param, climateVar = input_var)
  
  analysis_1d = histEX(data_in = input_data, param = input_param, 
                       climateVar = input_var)
  
  
  analysis_1e = histEX_plot(data_in = analysis_1d, param = input_param, 
                            climateVar = input_var)
  
  out1e = {
    datPath = file.path(createPaths(), 'RData', 'Excursion_histogram.RData')
    save(analysis_1e, file = datPath)
  }
  
  output = list(analysis_1e, analysis_1d)  
  
  return(output)
  
}

