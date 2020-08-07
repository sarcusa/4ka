MeanShift_prep <- function(input_data,input_param){
  
  print("Start apply meanshift function")
  analysis_2a = applyMS(data_in = input_data)
  
  out2a = {
    datPath = file.path(createPaths(), 'RData', 'MeanShift.RData')
    save(analysis_2a, file = datPath)
  }
  print("Completed apply meanshift function")
  
  print("Start meanshift null function")
  analysis_2b = MeanShift_null(data_in = analysis_2a, param = input_param)
  
  out_2b = {
    datPath = file.path(createPaths(), 'RData', 'MS_results_plusNull_complete.RData')
    save(analysis_2b, file = datPath)
  }
  print("Completed meanshift null function")
  
  return(analysis_2b)
  
}