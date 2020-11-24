dir = "/projects/pd_lab/sha59/4ka"
setwd(dir)

#envir <- new.env(parent = globalenv())

source('packages.R')
source('functions.R')
source('set_parameters.R')
source('plan_prep.R')
source('plan_prep_excursion.R')
source('plan_prep_meanshift.R')
source('plan_prep_trendchanges.R')
source('plan_excursion.R')
source('plan_meanshift.R')
source('plan_trendchanges.R')
source('plan_var.R')

library(doParallel)

parameters = list(eventYrs = eventYrs, event_window = event_window, 
                         ref_window = ref_window, resCriteria = resCriteria, 
                         sigNum = sigNum, plotOpt = plotOpt, 
                         mainDir = mainDir, numIt = numIt, 
                         res = res, radius = radius, 
                         eventWindow = eventWindow, 
                         CName = CName, CVers = CVers, 
                         eventDetector = eventDetector, OutDat = OutDat,
                         ncores = ncores)

param = parameters

skipped_records  <- vector()

#Manual
list2env(loading("/projects/pd_lab/sha59/4ka/RData/BS_results_complete.RData"),envir=.GlobalEnv)

for (i in 801:850) {
  print(paste0('RECORD ', i))
  
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
  
  registerDoParallel(cores = Sys.getenv("SLURM_CPUS_PER_TASK"))
  #registerDoParallel(cores = param$ncores)
  
  #cp  <- list()
  #cp.se <- list()
  #brk_dirs  <- list()
  
  cp.out <-foreach(it=1:param$numIt,
                   .verbose=T,.errorhandling = "pass") %dopar% {
                     
                     #print(paste0('BS null ITERATION ', it))
                     #for (it in 1:param$numIt) {
                     results = iterativeBrokenStick(TS_BS[[i]]$age, 
                                                    synthDat[,it], 
                                                    plot.opt = F )
                     #  TS_BS[[i]]$null_brk_pts[[it]] = results$cp
                     #  TS_BS[[i]]$null_brk_ptsErr[[it]] = results$cp.se
                     #}
                     
                     # added on 18-08-20 SA
                     if(is.character(results$cp) | is.na(results$cp)){
                       CP = NA #was cp
                       CP.SE = NA #was cp.se
                       brk_dirs = NA
                       
                     }else{
                       
                       CP = results$cp #was cp
                       CP.SE = results$cp.se #was cp.se
                       
                       slopes = results$o$coefficients[3:(length(results$cp)+2)]
                       brk_dirs = ifelse(slopes > 0, 1, -1)
                       
                       
                     }
                     
                     return(list(c = CP,  #was cp for both
                                 SE = CP.SE,  #was cp.se for both
                                 brk = brk_dirs)) #was brk_dirs
                   }    
  stopImplicitCluster()
  
  
  TS_BS[[i]]$null_brk_pts = lapply(cp.out, `[[`, "c") #was "cp" instead of 1
  TS_BS[[i]]$null_brk_ptsErr = lapply(cp.out, `[[`, "SE") #was "cp.se"
  # added on 18-08-20 SA
  TS_BS[[i]]$null_brk_dirs = lapply(cp.out, `[[`, "brk") #was "brk_dirs"
  
  
} #end of for loop

fileName = file.path(mainDir, 'RData', 'BS_results_plusNull_17.RData')
save(TS_BS, file = fileName)
write.table(skipped_records, file = file.path(mainDir,"skipped_records_BS.txt"))
