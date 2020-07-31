BrokenStick_null <- function(data_in, param){
  
  print("Running brockenStick_null now")
  
  TS_BS <- data_in[[1]]
  
  for (i in 1:length(TS_BS)) {
  
    print(paste0('RECORD ', i))
    
    TS_BS[[i]]$null_brk_pts = list()
    TS_BS[[i]]$null_brk_ptsErr = list()
    
    synthDat = try(createSyntheticTimeseries(TS_BS[[i]]$age,
                                         TS_BS[[i]]$paleoData_values,
                                         nens = param$numIt))
    if(class(synthDat) == "try-error"){
      
      print(paste0("Had to skip record = ",TS_BS[[i]]$dataSetName))
      TS_BS[[i]]$null_brk_pts = rep(list(NA), param$numIt)
      TS_BS[[i]]$null_brk_ptsErr = rep(list(NA), param$numIt)
      
    }else{
    
      registerDoParallel(cores = 16)
      
    cp.out <-foreach(it=1:param$numIt,
                     .verbose=TRUE,.errorhandling = "pass") %dopar% {
                       
                       results = iterativeBrokenStick(TS_BS[[i]]$age, 
                                                      synthDat[,it], 
                                                      plot.opt = F)
                     
                       return(list(cp = results$cp, cp.se = results$cp.se))
                     }    
    stopImplicitCluster()
    }
    
    TS_BS[[i]]$null_brk_pts = lapply(cp.out, `[[`, 1)
    TS_BS[[i]]$null_brk_ptsErr = lapply(cp.out, `[[`, 2)
    
  } #end of for loop
  save(TS_BS, file = 'BS_results_plusNull.RData')
  
  histBreaks = seq(100, 11500, by = 200)
  yrs = seq(300, 11500, by = 600)
  yrStr = yrs / 1000
  
  p = ggplot() + geom_histogram(aes(histBreaks), breaks = yrs) # histBreaks used to be breakPts
  pg = ggplot_build(p)
  densityBS = pg$data[[1]]$density
  xvals = pg$data[[1]]$x
  
  s  <- ggplot() + geom_col(aes(x = xvals, y = densityBS)) + 
    xlab('Event years [BP]') + ylab('Density of events') + 
    ggtitle('Broken Stick Results (simple plot)')
  
  output  <- list(TS_BS,s)
  return(output)
  
}
