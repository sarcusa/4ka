########
# Stability test
# 23-11-20
# 

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

load("/projects/pd_lab/sha59/4ka/RData/TS_climateInterp_2020.RData")

analysis_2a = applyMS(data_in = data, param = param)

param$eventYrs = param$eventYrs[1:25]
allEvents = rep(NA, length(param$eventYrs))
posEvents = rep(NA, length(param$eventYrs))
negEvents = rep(NA, length(param$eventYrs))
diffEvents = rep(NA, length(param$eventYrs))

for (y in 1:length(param$eventYrs)) {
  eventYr = param$eventYrs[y]
  TS_MS = analysis_2a
  
  for (i in 1:length(TS_MS)) {
  # Filter records that don't contain event year in age range
    if (min(TS_MS[[i]]$age) > eventYr || max(TS_MS[[i]]$age) < eventYr) {
      TS_MS[[i]]$useMS = 0
      #print('Out of range')
    }
    
    #Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
    if (length(TS_MS[[i]]$interpretation1_seasonalityGeneral) > 0) {
      if (tolower(TS_MS[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
            tolower(TS_MS[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
        TS_MS[[i]]$useMS = -1
        print(TS_MS[[i]]$interpretation1_seasonalityGeneral)
      }
    }
  }
  TS_MS = filterTs(TS_MS, 'useMS == 1')
  
  for (i in 1:length(TS_MS)) {
    
    
    TS_MS[[i]]$eventMS = 0
    TS_MS[[i]]$dirMS = 0
    if (!is.na(TS_MS[[i]]$sig_brks) && any(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow )) {
      TS_MS[[i]]$eventMS = 1
      eve_i = which(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow  & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow)
      TS_MS[[i]]$dirMS = TS_MS[[i]]$brk_dirs[eve_i]
    }
    
  }
  TS = TS_MS
  
  interps = unlist(sapply(TS,"[[","interpretation1_variable"))
  if (climateVar == 'M') {
    inds = which(interps == 'M' | interps == 'P' | interps == 'P-E' | interps ==  'P/E')
  } else if (climateVar == 'T') {
    inds = which(interps == 'T' | interps == 'TM')
  } else {
    inds = 1:length(interps)
  }
  
  events = as.numeric(sapply(TS,"[[","eventMS"))[inds]
  interps = interps[inds]
  
  dirs = unlist(sapply(TS,"[[","interpretation1_interpDirection"))[inds]
  dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
  dirs[dirs == 'negative'] = -1
  dirs[dirs == 'NA' | is.na(dirs)] = 0
  dirs = as.numeric(dirs)
  dirEvents = unlist(sapply(TS,"[[","dirMS"))[inds]
  dirChange = dirs * dirEvents 
  
  allEvents[y] = sum(events == 1) / length(events)
  posEvents[y] = sum(dirChange == 1) / length(events)
  negEvents[y] = sum(dirChange == -1) / length(events)
  diffEvents[y] = posEvents[y] - negEvents[y]
}

posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
negCol = ifelse(climateVar == 'M', '#543005', '#053061')
posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')

posDiff = diffEvents
negDiff = diffEvents
posDiff[diffEvents < 0] = 0
negDiff[diffEvents > 0] = 0

ggplot() +
  geom_col(aes(x = param$eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = param$eventYrs, y = negDiff), fill = negCol)

a = cbind(posDiff,negDiff)
b = cbind(posDiff,negDiff)
c = cbind(posDiff,negDiff)

# Data are stable when run this way, but not when run in full mode -- why?

name <- c("A","B","C","D","E","F","G","H","I","J")

stab <- list()

for(i in 1:length(name)){
  
  n <- name[i]
  
  load(paste0("/projects/pd_lab/sha59/4ka/Stability_test/",n,"/RData/MS_results_plusNull_complete.RData"))
  
  stab[[i]]  <- data_MS
  
}



# How many iterations until null results are stable?


eventDetector == 'MS'
stab_res <- list()

for(k in 1:length(name)){
  
  allNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  posNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  negNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  diffNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  allEvents = rep(NA, length(param$eventYrs))
  posEvents = rep(NA, length(param$eventYrs))
  negEvents = rep(NA, length(param$eventYrs))
  diffEvents = rep(NA, length(param$eventYrs))
  
  recordCounts = matrix(NA, nrow = length(param$eventYrs), ncol = 4) # years, all records, + events, - events
  recordCounts[,1] = param$eventYrs
  
  for (y in 1:length(param$eventYrs)) {
    
    eventYr = param$eventYrs[y]
    print(paste('Year:', eventYr))
    
    ## MEAN SHIFT PROCESSING
    
    TS_MS = stab[[k]]
    
    for (i in 1:length(TS_MS)) {
      # Filter records that don't contain event year in age range
      if (min(TS_MS[[i]]$age) > eventYr || max(TS_MS[[i]]$age) < eventYr) {
        TS_MS[[i]]$useMS = 0
        #print('Out of range')
      }
      
      #Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
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
      if (!is.na(TS_MS[[i]]$sig_brks) && any(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow )) {
        TS_MS[[i]]$eventMS = 1
        eve_i = which(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow  & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow)
        TS_MS[[i]]$dirMS = TS_MS[[i]]$brk_dirs[eve_i]
      }
      
    }
    TS = TS_MS
    
    # isolate only the records corresponding to the chosen climate interpretation
    interps = unlist(sapply(TS,"[[","interpretation1_variable"))
    if (climateVar == 'M') {
      inds = which(interps == 'M' | interps == 'P' | interps == 'P-E' | interps ==  'P/E')
    } else if (climateVar == 'T') {
      inds = which(interps == 'T' | interps == 'TM')
    } else {
      inds = 1:length(interps)
    }
        
    events = as.numeric(sapply(TS,"[[","eventMS"))[inds]
    
    interps = interps[inds]
    
    # calculate the climate event direction based on the proxy climate dir and event dir: checks if the interpretation direction is the same as in the new
    dirs = unlist(sapply(TS,"[[","interpretation1_interpDirection"))[inds]
    dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
    dirs[dirs == 'negative'] = -1
    dirs[dirs == 'NA' | is.na(dirs)] = 0
    dirs = as.numeric(dirs)
    dirEvents = unlist(sapply(TS,"[[","dirMS"))[inds] # (0, 1, -1): (no, positive, negative) event
    dirChange = dirs * dirEvents# (0, 1, -1): (no, positive, negative) climate event
    
    # add these for printing table of results by site
    names = unlist(sapply(TS, '[[', 'dataSetName'))[inds]
    SiteName = unlist(sapply(TS, '[[', 'geo_siteName'))[inds]
    lat = unlist(sapply(TS, '[[', 'geo_latitude'))[inds]
    lon = unlist(sapply(TS, '[[', 'geo_longitude'))[inds]
    ArchiveType = unlist(sapply(TS, '[[', 'archiveType'))[inds]
    ProxyType = unlist(sapply(TS, '[[', 'paleoData_proxyGeneral'))[inds]
    for (i in 1:length(TS)) {
      if (length(TS[[i]]$pub1_doi) == 0) {
        TS[[i]]$pub1_doi = NA
      }
    }
    for (i in 1:length(TS)) {
      if (length(TS[[i]]$originalDataUrl) == 0) {
        TS[[i]]$originalDataUrl = NA
      }
    } 
    for(i in 1:length(TS)){ # problem with this specific reference
      if(TS[[i]]$paleoData_TSid == "RyLPf9ScmVz"){
        TS[[i]]$pub1_doi = NA
      }
    }
    doi = unlist(sapply(TS, '[[', 'pub1_doi'))[inds]
    TSid = unlist(sapply(TS, '[[', 'paleoData_TSid'))[inds]
    Original = unlist(sapply(TS, '[[', 'originalDataUrl'))[inds]
    
    # save table of results by site
    event_by_site_df = data.frame(TVerseID = TSid, Dataset = names, Site = SiteName, Lat = lat, Lon = lon, Event = dirChange,Year = rep(eventYr, length(names)), Var = rep(climateVar, length(names)), Archive = ArchiveType, Proxy = ProxyType, OriginalURL = Original, DOI = doi)
        
    # store real events summary for the year
    allEvents[y] = sum(events == 1) / length(events)
    posEvents[y] = sum(dirChange == 1) / length(events)
    negEvents[y] = sum(dirChange == -1) / length(events)
    diffEvents[y] = posEvents[y] - negEvents[y]
    recordCounts[y,2:4] = c(length(events), posEvents[y], negEvents[y])
    
    # Assign event occurrence
    totNullEvents = matrix(0, nrow = length(inds), ncol = param$numIt) 
    for (i in 1:length(inds)) {
      
    nullBreaks = TS[[inds[i]]]$null_sig_brks
    nullDirs = TS[[inds[i]]]$null_brk_dirs
      
      for (j in 1:param$numIt) {
        
        #sensitivity testing
        
        eventInd = which(nullBreaks[[j]] >= eventYr - param$eventWindow  & nullBreaks[[j]] <= eventYr + param$eventWindow )
        
        if (length(eventInd) > 0) {
          
          if (eventDetector == 'MS') {
            totNullEvents[i,j] = nullDirs[[j]][eventInd[1]] * dirs[i]
          } else {
            # Metric combining events occurring in event window (event takes the sign of the sum)
            dir = sum(nullDirs[[j]][eventInd]) / abs(sum(nullDirs[[j]][eventInd]))
            totNullEvents[i,j] = ifelse(is.na(dir), 0, dir) * dirs[i]
          }
          
        }
        
      } # end it loop
      
    } # end record loop
    
    allNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x != 0)) / length(inds)
    posNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == 1)) / length(inds)
    negNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == -1)) / length(inds)
    
  } # end event year loop 
  
  allNullQuants = apply(allNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  posNullQuants = apply(posNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  negNullQuants = apply(negNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  diffNullQuants = apply(posNullEvents - negNullEvents, 1, function(x) quantile(x, probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99)))
  
  res <- list(
    allEvents = allEvents,
    posEvents = posEvents,
    negEvents = negEvents,
    diffEvents = diffEvents,
    recordCounts = recordCounts,
    allNullEvents = allNullEvents,
    posNullEvents = posNullEvents,
    negNullEvents = negNullEvents,
    allNullQuants = allNullQuants,
    posNullQuants = posNullQuants,
    negNullQuants = negNullQuants,
    diffNullQuants = diffNullQuants
    )
  stab_res[[k]] <- res
  
} #end of folder loop

eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))
eventDetectorStr = ifelse(eventDetector == 'MS', 'Mean shift', 'Broken stick')

posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
negCol = ifelse(climateVar == 'M', '#543005', '#053061')
posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')

eventYrs = param$eventYrs[1:25]

m  <- 1
posDiff = stab_res[[m]]$diffEvents
negDiff = stab_res[[m]]$diffEvents
posDiff[stab_res[[m]]$diffEvents < 0] = 0
negDiff[stab_res[[m]]$diffEvents > 0] = 0
diffNullQuants <- stab_res[[m]]$diffNullQuants

ggplot() +
  geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = eventYrs, y = negDiff), fill = negCol)

ggplot() +
  geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3])   

# Checks for stability of the data results first, then check for the stability of number of iterations

posNullQuants = apply(posNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
negNullQuants = apply(negNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
diffNullQuants = apply(posNullEvents - negNullEvents, 1, function(x) quantile(x, probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99)))

DE_all <- apply(stab_res,"[[","diffEvents")
PNE_all <- apply(stab_res, "[[", "posNullEvents")
NNE_all <- apply(stab_res, "[[", "negNullEvents")

# Measure of significant events

sig <- matrix(NA,nrow = length(eventYrs), ncol = length(name))

for(o in 1:length(name)){
  
  m <- seq(1000,length(name)*1000,1000)
  DE <- DE_all[,o]
  posDiff = DE
  negDiff = DE
  posDiff[DE < 0] = 0
  negDiff[DE > 0] = 0
    
  DNQ <- apply(PNE_all[,1:m[o]] - NNE_all[,1:m[o]], 1, function(x) quantile(x, probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99)))
  
  for(i in 1:length(eventYrs)){
    
    if(DE[i] >= DNQ[6,i] | DE[i] <= DNQ[3,i]){
      
      sig[i,o] = 1
      
    }else{
    
      sig[i,o] = 0
      
    }
    
  }
  
}

tot <- colSums(sig)

plot(m, tot, ylab = "Total # events (1-10.6 ka)", xlab = "# iterations")

# Individual iterations
# 1000 iterations
DNQ1 <- apply(PNE_all[,1:1000] - NNE_all[,1:1000], 1, function(x) quantile(x, probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99))) #not sure on the size
DE1 <- DE_all[,1]

posDiff = DE1
negDiff = DE1
posDiff[DE1 < 0] = 0
negDiff[DE1 > 0] = 0
diffNullQuants <- DNQ1

ggplot() +
  geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3]) 

