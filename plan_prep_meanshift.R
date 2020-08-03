MeanShift_prep <- function(input_data,input_param){
  
  analysis_2a = applyMS(data_in = input_data)
  
  out2a = {
    datPath = file.path(createPaths(), 'RData', 'MeanShift.RData')
    save(analysis_2a, file = datPath)
  }
  
  analysis_2b = MeanShift_null(data_in = analysis_2a, param = input_param)
  
  out_2b = {
    datPath = file.path(createPaths(), 'RData', 'MS_results_plusNull_complete.RData')
    save(analysis_2b, file = datPath)
  }
  
  return(analysis_2b)
  
}