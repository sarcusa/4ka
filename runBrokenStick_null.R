BrokenStick_null <- function(data_in, param){
  
  print("Running brockenStick_null now")
  
  TS_BS <- data_in
  skipped_records  <- vector()
  
  #plan(cluster)
  
  for (i in 1:length(TS_BS)) {
    #for (i in 1:50) {
    #print(paste0('RECORD ', i))
    
    TS_BS[[i]]$null_brk_pts = list()
    TS_BS[[i]]$null_brk_ptsErr = list()
    
    # added on 18-08-20 SA
    TS_BS[[i]]$null_brk_dirs = list()
    
    synthDat = try(createSyntheticTimeseries(time = TS_BS[[i]]$age,
                                             values = TS_BS[[i]]$paleoData_values,
                                         n.ens = param$numIt))
    if(class(synthDat) == "try-error"){
      
      synthDat = try(createSyntheticTimeseries(time = TS_BS[[i]]$age,
                                       values = TS_BS[[i]]$paleoData_values,
                                               n.ens = param$numIt))
      if(class(synthDat) == "try-error"){
        
        print(paste0("Had to skip record = ",TS_BS[[i]]$dataSetName))
        skipped_records[i]  <- TS_BS[[i]]$dataSetName
        TS_BS[[i]]$useBS = -1
        
        next 
      }
    }
    
    registerDoParallel(cores = param$ncores)
      
    #cp  <- list()
    #cp.se <- list()
    #brk_dirs  <- list()
      
    cp.out <-foreach(it=1:param$numIt,
                     .verbose=F,.errorhandling = "pass") %do% {
                       
                       #print(paste0('BS null ITERATION ', it))
                       #for (it in 1:param$numIt) {
                       results = iterativeBrokenStick(TS_BS[[i]]$age, 
                                                      synthDat[,it], 
                                                      plot.opt = F)
                     #  TS_BS[[i]]$null_brk_pts[[it]] = results$cp
                     #  TS_BS[[i]]$null_brk_ptsErr[[it]] = results$cp.se
                     #}
                       
                       # added on 18-08-20 SA
                       if(is.character(results$cp) | is.na(results$cp)){
                         cp = NA 
                         cp.se = NA
                         brk_dirs = NA
                       
                       }else{
                         
                         cp = results$cp
                         cp.se = results$cp.se
                         
                         slopes = results$o$coefficients[3:(length(results$cp)+2)]
                         brk_dirs = ifelse(slopes > 0, 1, -1)
                         
                         
                       }
                                              
                       return(list(cp = cp, 
                                   cp.se = cp.se, 
                                   brk_dirs = brk_dirs))
                     }    
    stopImplicitCluster()
    
    
    TS_BS[[i]]$null_brk_pts = lapply(cp.out, `[[`, "cp")
    TS_BS[[i]]$null_brk_ptsErr = lapply(cp.out, `[[`, "cp.se")
    # added on 18-08-20 SA
    TS_BS[[i]]$null_brk_dirs = lapply(cp.out, `[[`, "brk_dirs")
    
    
  } #end of for loop
  
  fileName = file.path(mainDir, 'RData', 'BS_results_plusNull.RData')
  save(TS_BS, file = fileName)
  write.table(skipped_records, file = file.path(mainDir,"skipped_records_BS.txt"))
  
  return(TS_BS)
  
}
