TrendChanges_prep <- function(input_data, input_param){
  
  analysis_3a = BrokenStick(input_data, input_param)
  
  out3a = {
    datPath = file.path(createPaths(), 'RData', 'TrendChanges.RData')
    save(analysis_3a, file = datPath)
  }
  
  analysis_3b = BrokenStick_null(analysis_3a, input_param)
  
  out_3b = {
    datPath = file.path(createPaths(), 'RData', 'BS_results_plusNull_complete.RData')
    save(analysis_3b, file = datPath)
  }

  return(analysis_3b)
}