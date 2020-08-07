Excursion_prep <- function(input_data,input_param){
  
  print("start apply excursion function")
    analysis_1a = applyEX(data_in = input_data, param = input_param)
    
    out1a = {
      datPath = file.path(createPaths(), 'RData', 'Excursion.RData')
      save(analysis_1a, file = datPath)
    }
  print("completed apply excursion function")  
  
  print("start excursion null function")
    analysis_1b = Excursion_null(data_in = analysis_1a, param = input_param)
    
    out1b = {
      datPath = file.path(createPaths(), 'RData', 'Excursion_plusnull.RData')
      save(analysis_1b, file = datPath)
    }
    print("completed excursion null function")
    
  return(analysis_1b)
  
}
