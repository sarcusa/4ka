#############
# Combine all parameter choices for each detector (rather than the central parameters) for spatial maps
# December 2020 Arcusa

dir = "/projects/pd_lab/sha59/4ka" #set the directory of the input scripts
setwd(dir)

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
#source('plan_plotting.R')
source('plan_var.R')

parameters = list(eventYrs = eventYrs, event_window = event_window, 
                  ref_window = ref_window, resCriteria = resCriteria, 
                  sigNum = sigNum, plotOpt = plotOpt, 
                  mainDir = mainDir, numIt = numIt, 
                  res = res, radius = radius, 
                  eventWindow = eventWindow, 
                  CName = CName, CVers = CVers, 
                  eventDetector = eventDetector, OutDat = OutDat,
                  ncores = ncores, maxDiff = maxDiff)

param = parameters

# lists folders in directory, and separates for each detector
listdirs <- list.dirs(path = dir, full.names = T,recursive = F)
seldirs <- listdirs[grep(pattern = "Signif", x = listdirs)]
seldirs <- seldirs[-grep(pattern = "random", seldirs)]
seldirs <- seldirs[-grep(pattern = "seq", seldirs)]
seldirsEX <- c(seldirs[grep(pattern = "EX", seldirs)])
seldirsMS <- c(seldirs[grep(pattern = "MS", seldirs)])

#' Creates the output needed for MS spatial maps
#' @description This function takes the output of the null analysis and calculates the spatial grids
#' @param D list of directories
#' @param param a list of parameters used in the analysis
#' @param climateVar climate variables (T or M)

GridCreate <- function(D, param, climateVar) {
  
  for(k in 1:length(D)){
    
    list2env(loading(paste0(D[k], '/RData/MS_results_plusNull_complete.RData')),envir=.GlobalEnv)
    #load(paste0(D[k],"RData/MS_results_plusNull_complete.RData")
    
    for (y in 1:length(param$eventYrs[1:25])) {
      
      TS_MS = data_MS
      eventYr = param$eventYrs[y]
      
      for (i in 1:length(TS_MS)) {
        
        # Filter records that don't contain event year in age range
        if (min(TS_MS[[i]]$age) > eventYr || max(TS_MS[[i]]$age) < eventYr) {
          TS_MS[[i]]$useMS = 0
          #print('Out of range')
        }
        
        # Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
        if (length(TS_MS[[i]]$interpretation1_seasonalityGeneral) > 0) {
          if (tolower(TS_MS[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
                tolower(TS_MS[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
            TS_MS[[i]]$useMS = -1
            print(TS_MS[[i]]$interpretation1_seasonalityGeneral)
          }
        }
      }
      TS_MS = filterTs(TS_MS, 'useMS == 1')
      
      # Assign event occurrence
      for (i in 1:length(TS_MS)) {
        
        TS_MS[[i]]$eventMS = 0
        TS_MS[[i]]$dirMS = 0
        if (!is.na(TS_MS[[i]]$sig_brks) && any(TS_MS[[i]]$sig_brks >= eventYr -  param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr +  param$eventWindow)) {
          TS_MS[[i]]$eventMS = 1
          eve_i = which(TS_MS[[i]]$sig_brks >= eventYr -  param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr +  param$eventWindow)
          TS_MS[[i]]$dirMS = TS_MS[[i]]$brk_dirs[eve_i]
        }
        
      }
      
      # isolate only the records corresponding to the chosen climate interpretation
      interps = unlist(sapply(TS_MS,"[[","interpretation1_variable"))
      if (climateVar == 'M') {
        inds = which(interps == 'M' | interps == 'P' | interps == 'P-E' | interps ==  'P/E')
      } else {
        inds = which(interps == 'T' | interps == 'TM')
      }
      
      allTsLat = as.numeric(sapply(TS_MS,"[[","geo_latitude"))[inds]
      allTsLon = as.numeric(sapply(TS_MS,"[[","geo_longitude"))[inds]
      events = as.numeric(sapply(TS_MS,"[[","eventMS"))[inds]
      interps = interps[inds]
      
      # calculate the climate event direction based on the proxy climate dir and event dir
      dirs = unlist(sapply(TS_MS,"[[","interpretation1_interpDirection"))[inds]
      dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
      dirs[dirs == 'negative'] = -1
      dirs[dirs == 'NA' | is.na(dirs)] = 0
      dirs = as.numeric(dirs)
      dirEvents = unlist(sapply(TS_MS,"[[","dirMS"))[inds]  #(0, 1, -1):(no, positive, negative) event
      dirChange = dirs * dirEvents #(0, 1, -1):(no, positive, negative) climate event
      
      #### THIS IS NEW, NOT WORKING
      dirChange <- matrix(unlist(dirChange), ncol = 9)
      ####
      
      
      ## Spatial null
      
      # Make a grid of lat/lon center points
      latitude = seq(-90 + param$res/2, 90 - param$res/2, param$res)
      longitude = seq(-180 + param$res/2, 180 - param$res/2, param$res)
      
      # Grids that will contain fraction of passes
      gridEvents = matrix(NA, nrow = length(longitude), ncol = length(latitude))
      gridNumRecs = matrix(NA, nrow = length(longitude), ncol = length(latitude))
      gridNullEvents = array(NA, dim = c(length(longitude), 
                                         length(latitude), param$numIt))
      totNullEvents = matrix(0, nrow = length(inds), ncol = param$numIt)
      gridPercentEvents_null = matrix(NA, nrow = length(longitude), 
                                      ncol = length(latitude))
      
      # Assign event occurrence
      for (i in 1:length(inds)) {
        
        nullBreaks = TS_MS[[inds[i]]]$null_sig_brks
        nullDirs = TS_MS[[inds[i]]]$null_brk_dirs
        
        # if no breaks are found in any iteration, nullBreaks is just an empty list and will break code
        if (length(nullBreaks) > 0) {
          for (j in 1:param$numIt) {
            
            eventInd = which(nullBreaks[[j]] >= eventYr - param$eventWindow & nullBreaks[[j]] <= eventYr + param$eventWindow)
            
            if (length(eventInd) > 0) {
              totNullEvents[i,j] = nullDirs[[j]][eventInd[1]] * dirs[i]
            }
            
          }  # end it loop (synth data)
        } 
      }  # end event occurrence assignment
      
      for (i in 1:length(longitude)) {
        
        print(paste('iteration i = ', i))
        
        for (j in 1:length(latitude)) {
          
          indsLoc = which(distGeo(cbind(allTsLon, allTsLat), c(longitude[i],latitude[j])) <= param$radius*1000)
          
          if (length(indsLoc) > 0) { 
            # number of events = # pos events - # neg events, must have > 1 record contributing to grid cell
            gridNumRecs[i, j] = length(indsLoc)
            gridEvents[i, j] = ifelse(length(indsLoc) > 1, sum(dirChange[indsLoc,] == 1) - sum(dirChange[indsLoc,] == -1), NA)
          }
          
          if (length(indsLoc) == 1) {
            gridNullEvents[i, j, ] = totNullEvents[indsLoc, ]
          } else {
            gridNullEvents[i, j, ] = apply(totNullEvents[indsLoc, ], 2, function (x) sum(x == 1) - sum(x == -1))
          }
          
          if (length(indsLoc) > 0) {
            if (is.na(gridEvents[i, j])) {
              next
            } else if (gridEvents[i, j] < 0) {
              gridPercentEvents_null[i,j] = -sum(gridNullEvents[i,j,] > gridEvents[i, j]) / param$numIt
            } else if (gridEvents[i, j] > 0) {
              gridPercentEvents_null[i,j] = sum(gridNullEvents[i,j,] < gridEvents[i, j]) / param$numIt
            } else {
              gridPercentEvents_null[i,j] = 0
            }
          }
        }
        
      } # end location grid loops
      
      OP <- list(gridNumRecs,gridEvents,gridNullEvents,totNullEvents,allTsLon, allTsLat, dirChange)
      save(OP, file = file.path(D[k], paste0("RData/results_spatial_MS_", 
                                             climateVar, "_", 
                                             eventYr/1000, ".RData")))         
      
    }# end of eventYr loop
    
    
  }# end of D loop
  
  
}

future::plan(batchtools_slurm, template = file.path(dir,"slurm_batchtools_spatial.tmpl"), resources = list(cores=Sys.getenv("SLURM_CPUS_PER_TASK")))

future( {
  GridCreate(D = seldirsMS[1], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[1], param = param, climateVar = "T")
)}

future( {
  GridCreate(D = seldirsMS[2], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[2], param = param, climateVar = "T")
)}

future( {
  GridCreate(D = seldirsMS[3], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[3], param = param, climateVar = "T")
)}

future( {
  GridCreate(D = seldirsMS[4], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[4], param = param, climateVar = "T")
)}
future( {
  GridCreate(D = seldirsMS[6], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[6], param = param, climateVar = "T")
)}
future( {
  GridCreate(D = seldirsMS[7], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[7], param = param, climateVar = "T")
)}
future( {
  GridCreate(D = seldirsMS[8], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[8], param = param, climateVar = "T")
)}
future( {
  GridCreate(D = seldirsMS[9], param = param, climateVar = "M")
  GridCreate(D = seldirsMS[9], param = param, climateVar = "T")
)}

