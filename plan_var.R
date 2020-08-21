my_plan <- drake_plan(
  
  trace = TRUE,
                 
  parameters = target(list(eventYrs = eventYrs, event_window = event_window, 
                    ref_window = ref_window, resCriteria = resCriteria, 
                    sigNum = sigNum, plotOpt = plotOpt, 
                    mainDir = mainDir, numIt = numIt, 
                    res = res, radius = radius, 
                    eventWindow = eventWindow, 
                    CName = CName, CVers = CVers, 
		    eventDetector = eventDetector, OutDat = OutDat,
	            ncores = ncores), hpc = FALSE),

  data = target(plan_prep(param = parameters), hpc = FALSE),

  map = target(plot_sites(data), hpc = FALSE),

  #prep_1 = target(Excursion_prep(input_data = data,input_param = parameters),
  #                resources = list(cores = 16)),  
  
  #analysis_1 = target(Excursion(input_data = prep_1,input_param = parameters,
  #                       input_var = climate_var),
  #                    transform = cross(
  #                      climate_var = c("T","M","All")),
  #                    resources = list(cores = 16)),
  
  #prep_2 = target(MeanShift_prep(input_data = data,input_param = parameters),
  #                resources = list(cores = 32)), 
  
  #analysis_2 = target(MeanShift(input_data = prep_2, input_param = parameters,
  #                       input_var = climate_var),
  #                    transform = cross(climate_var = c("T","M","All")),
  #                    resources = list(cores = 16)),
  
  prep_3a = target(TrendChanges_prep1(input_data = data,input_param = parameters),
                  resources = list(cores = 16)),
  prep_3b = target(TrendChanges_prep2(input_data = prep_3a,
                                      input_param = parameters),
                   resources = list(cores = 32)),
  
  #test = target(BrokenStick_null_test(data_in = prep_3a, param = parameters),
  #              resources = list(cores = 16))
  
  #analysis_3 = target(TrendChanges(input_data = prep_3b,input_param = parameters,
  #                          input_var = climate_var), 
  #                    transform = cross(climate_var = c("T", "M")),
  #                    resources = list(cores = 16)),
  
  #results_1 = target(ProxyMap(prep1 = prep_1, prep2 = prep_2,
  #                            param = parameters, input_var = climate_var),
  #                   hpc = FALSE,
  #                   transform = cross(climate_var = c("T", "M"))),
  
  #results_2 = target(histogram_net(EX = analysis_1, MS = analysis_2, 
  #                                 BS = analysis_3, param = parameters), 
  #                   transform = map(analysis_1, analysis_2, analysis_3),
  #                   hpc = FALSE),
  
  #report  = target(rmarkdown::render(
  #  knitr_in("report.Rmd"),
  #  output_file = file_out("report.html"),
  # quiet = TRUE), hpc = FALSE)
  
)

