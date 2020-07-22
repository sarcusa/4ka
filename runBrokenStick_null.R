BrokenStick_null <- function(data_in, param){
  
  TS_BS <- data_in[[1]]
  
  for (i in 1:length(TS_BS)) {
    
    print(paste0('RECORD ', i))
    
    TS_BS[[i]]$null_brk_pts = list()
    TS_BS[[i]]$null_brk_ptsErr = list()
    
    tryCatch({synthDat = createSyntheticTimeseries(TS_BS[[i]]$age, 
                                         TS_BS[[i]]$paleoData_values, 
                                         nens = param$numIt)
    
    # run the broken stick code for all iterations
    for (it in 1:param$numIt) {
      print(paste0('ITERATION ', it))
      
      results = iterativeBrokenStick(TS_BS[[i]]$age, synthDat[,it], 
                                     plot.opt = F)
      
      TS_BS[[i]]$null_brk_pts[[it]] = results$cp
      TS_BS[[i]]$null_brk_ptsErr[[it]] = results$cp.se
    }}, 
    error = function(e){cat("Error:", conditionMessage(e), 
                            paste0(" and had to skip record", i))})
    
  }
  
  save(TS_BS, file = 'BS_results_plusNull.RData')
  
  histBreaks = seq(100, 11500, by = 200)
  yrs = seq(300, 11500, by = 600)
  yrStr = yrs / 1000
  
  p = ggplot() + geom_histogram(aes(histBreaks), breaks = yrs) # breakPts
  pg = ggplot_build(p)
  densityBS = pg$data[[1]]$density
  xvals = pg$data[[1]]$x
  
  s  <- ggplot() + geom_col(aes(x = xvals, y = densityBS)) + 
    xlab('Event years [BP]') + ylab('Density of events') + 
    ggtitle('Broken Stick Results (simple plot)')
  
  output  <- list(TS_BS,s)
  return(output)
  
}