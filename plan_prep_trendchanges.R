TrendChanges_prep <- function(input_data, input_param){
  
  print("Start brocken stick function")
  analysis_3a = BrokenStick(input_data, input_param)
  
  out3a = {
    datPath = file.path(createPaths(), 'RData', 'TrendChanges.RData')
    save(analysis_3a, file = datPath)
  }
  print("completed brocken stick function")
  
  print("start brockenstick null function")
  analysis_3b = BrokenStick_null(analysis_3a, input_param)
  
  out_3b = {
    datPath = file.path(createPaths(), 'RData', 'BS_results_plusNull_complete.RData')
    save(analysis_3b, file = datPath)
  }
  print("completed brockenstick null function")

  
  return(analysis_3b)
}