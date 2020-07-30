my_plan <- drake_plan(
  
  trace = TRUE,
  
  data = target(plan_prep(), hpc = FALSE),
                 
  map = target(plot_sites(data), hpc = FALSE),
                 
  parameters = target(list(eventYrs = eventYrs, event_window = event_window, 
                    ref_window = ref_window, resCriteria = resCriteria, 
                    sigNum = sigNum, plotOpt = plotOpt, 
                    mainDir = mainDir, numIt = numIt, 
                    res = res, radius = radius, 
                    eventWindow = eventWindow, 
                    eventDetector = eventDetector), hpc = FALSE),
  
  analysis_1 = target(Excursion(input_data = data,input_param = parameters,
                         input_var = climate_var),
                                         transform = cross(
                                           climate_var = c("M", "T", "all")),
                      resources = list(cores = 4)),
  
  
  analysis_2 = target(MeanShift(input_data = data, input_param = parameters,
                         input_var = climate_var),
                      transform = cross(climate_var = c("M", "T", "all")),
                      resources = list(cores = 4)),
  
  analysis_3 = target(TrendChanges(input_data = data,input_param = parameters,
                            input_var = climate_var), 
                      transform = cross(climate_var = c("M", "T")),
                      resources = list(cores = 16)),
  
  results_1 = target(Plotting(analysis1 = analysis_1, analysis2 = analysis_2,
                              analysis3 = analysis_3, 
                              input_param = parameters, input_var = climate_var),
                     transform = map(analysis_1, analysis_2, analysis_3),
                     resources = list(cores = 1)),
  
  report  = target(rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
   quiet = TRUE), hpc = FALSE)
  
)

