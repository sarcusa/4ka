my_plan <- drake_plan(
  
  trace = TRUE,
  
  max_expand = 2,
    
  data = plan_prep(),
  
  map = plot_sites(data),
  
  parameters = list(eventYrs = eventYrs, event_window = event_window, 
                    ref_window = ref_window, resCriteria = resCriteria, 
                    sigNum = sigNum, plotOpt = plotOpt, 
                    mainDir = mainDir, numIt = numIt, 
                    res = res, radius = radius, 
                    eventWindow = eventWindow, climateVar = climateVar,
                    eventDetector = eventDetector),
  
  analysis_1 = Excursion(input_data = data,input_param = parameters,
                         input_var = parameters$climateVar),
                      
  
  analysis_2 = MeanShift(input_data = data, input_param = parameters,
                                input_var = parameters$climateVar),
                        
  analysis_3 = TrendChanges(input_data = data,input_param = parameters,
                                   input_var = parameters$climateVar),
    
  #results_1 = histogram_net(data_EX = analysis_1, 
  #                        data_MS = analysis_2,
  #                        data_BS = analysis_3, param = parameters),
  
                     
  results_2 = proxyMap(data_MS = analysis_2,
                              data_EX = analysis_1, param = parameters),
     
  report  = rmarkdown::render(
      knitr_in("report.Rmd"),
      output_file = file_out("report.html"),
      quiet = TRUE)
  
)

#analysis_2 = target(MeanShift(input_data = data[[1]], input_param = parameters,
#                              input_var = climate_var),
#                    transform = cross(
#                      climate_var = c("M", "T")),
#  

