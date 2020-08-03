Plotting <- function(prep1,prep2,analysis1, analysis2, analysis3,input_param, input_var){
  
  plots_1 = histogram_net(EX = analysis1, MS = analysis2, 
                          BS = analysis3, param = input_param)
  
  plots_2 = ProxyMap(data_MS = prep2, data_EX = prep1, 
                     param = input_param, climate_Var = input_var)
  
  output  <- list(plots_1[[1]], plots_1[[2]], plots_2[[1]], plots_2[[2]])
  
  return(output)
  
}

