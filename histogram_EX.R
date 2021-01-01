histEX <- function(data_in, param, climateVar){
  
  dataDir = file.path(createPaths(), 'RData')
  
  param$eventYrs = param$eventYrs[1:25]
  
  allNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  posNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  negNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  #diffNullEvents = matrix(NA, nrow = length(param$eventYrs), ncol = param$numIt)
  allEvents = rep(NA, length(param$eventYrs))
  posEvents = rep(NA, length(param$eventYrs))
  negEvents = rep(NA, length(param$eventYrs))
  diffEvents = rep(NA, length(param$eventYrs))
  
  recordCounts = matrix(NA, nrow = length(param$eventYrs), ncol = 4) # years, all records, + events, - events
  recordCounts[,1] = param$eventYrs
  
  #these are for the composite
  everyEvent <- list()
  everyDirChange <- list()
  everyTotNullEvent <- list()
  
  for (y in 1:length(param$eventYrs)) {
    
    remove(data_EX)
    eventYr = param$eventYrs[y]
    print(paste('Year:', eventYr))
    
    list2env(loading(file.path(dataDir, paste0('EX_results_plusNull_complete_', eventYr/1000, '.RData'))),envir=.GlobalEnv)
    
    if(!length(data_EX) == 0 & !length(data_EX) == 1){
    
    for (i in 1:length(data_EX)) {
      
    # Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
      if (length(data_EX[[i]]$interpretation1_seasonalityGeneral) > 0) {
        if (tolower(data_EX[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
              tolower(data_EX[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
          data_EX[[i]]$useEX = -1
          print(data_EX[[i]]$interpretation1_seasonalityGeneral)
        }
      }
    }
    TS = filterTs(data_EX, 'useEX == 1')
    
    # isolate only the records corresponding to the chosen climate interpretation
    interps = unlist(sapply(TS,"[[","interpretation1_variable"))
    if (climateVar == 'M') {
      inds = which(interps == 'M' | interps == 'P' | interps == 'P-E' | interps ==  'P/E')
    } else if (climateVar == 'T') {
      inds = which(interps == 'T' | interps == 'TM')
    } else {
      inds = 1:length(interps)
    }
    
    if(!length(inds) == 0){
    
    events = as.numeric(sapply(TS,"[[","eventEX"))[inds]
    interps = interps[inds]
    
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
    doi = unlist(sapply(TS, '[[', 'pub1_doi'))[inds]
    TSid = unlist(sapply(TS, '[[', 'paleoData_TSid'))[inds]
    Original = unlist(sapply(TS, '[[', 'originalDataUrl'))[inds]
    
    # calculate the climate event direction based on the proxy climate dir and event dir
    dirs = unlist(sapply(TS,"[[","interpretation1_interpDirection"))[inds]
    dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
    dirs[dirs == 'negative'] = -1
    dirs[dirs == 'NA' | is.na(dirs)] = 0
    dirs = as.numeric(dirs)
    dirEvents = unlist(sapply(TS,"[[","dirEx"))[inds]  #(0, 1, -1):(no, positive, negative) event
    dirChange = dirs * dirEvents  #(0, 1, -1):(no, positive, negative) climate event
    
    everyEvent[[y]] <- events
    everyDirChange[[y]] <- dirChange
    
    # save table of results by site
    event_by_site_df = data.frame(TVerseID = TSid, Dataset = names, Site = SiteName, Lat = lat, Lon = lon, Event = dirChange,Year = rep(eventYr, length(names)), Var = rep(climateVar, length(names)), Archive = ArchiveType, Proxy = ProxyType, OriginalURL = Original, DOI = doi)
    write.table(event_by_site_df, file = file.path(dataDir, paste0('event_each_site_EX_',climateVar,'.csv')), append = T,row.names = F, col.names = F, sep = ',')
    
    # store real events summary for the year
    allEvents[y] = sum(events == 1) / length(events)
    posEvents[y] = sum(dirChange == 1) / length(events)
    negEvents[y] = sum(dirChange == -1) / length(events)
    diffEvents[y] = posEvents[y] - negEvents[y]
    recordCounts[y,2:4] = c(length(events), posEvents[y], negEvents[y])
    
    # Assign event occurrence
    totNullEvents = matrix(0, nrow = length(inds), ncol = param$numIt) 
    for (i in 1:length(inds)) {
      
      if(length(TS[[inds[i]]]$null_events_dir) != param$numIt){
        
        add_NA <- rep(NA,param$numIt - length(TS[[inds[i]]]$null_events_dir))
        TS[[inds[i]]]$null_events_dir  <- c(TS[[inds[i]]]$null_events_dir,add_NA)
        totNullEvents[i,] = TS[[inds[i]]]$null_events_dir * dirs[i]
      }
      
      totNullEvents[i,] = TS[[inds[i]]]$null_events_dir * dirs[i]
      
    }
    
    allNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x != 0)) / length(inds)
    posNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == 1)) / length(inds)
    negNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == -1)) / length(inds)
    
    } else {
      #added so code will not break when no event is found
      print(paste0("no ",climateVar, " variable found"))
      allEvents[y] = 0
      posEvents[y] = 0
      negEvents[y] = 0
      diffEvents[y] = posEvents[y] - negEvents[y]
      allNullEvents[y,] = 0
      posNullEvents[y,] = 0
      negNullEvents[y,] = 0 
      recordCounts[y,2:4] = c(0, 0, 0)
    }# end of second ifelse
    
      
    }else{
      
      allEvents[y] = 0
      posEvents[y] = 0
      negEvents[y] = 0
      diffEvents[y] = posEvents[y] - negEvents[y]
      allNullEvents[y,] = 0
      posNullEvents[y,] = 0
      negNullEvents[y,] = 0
      recordCounts[y,2:4] = c(0, 0, 0)
      
    } # end of first ifelse 
      
   everyTotNullEvent[[y]] <- totNullEvents

    
  } # end event year loop
   
  allNullQuants = apply(allNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  posNullQuants = apply(posNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  negNullQuants = apply(negNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  diffNullQuants = apply(posNullEvents - negNullEvents, 1, function(x) quantile(x, probs = c(0.1, 0.05, 0.01, 0.9, 0.95, 0.99)))
  
  output <- list(
    everyEvent = everyEvent,
    everyDirChange = everyDirChange,
    allEvents = allEvents,
    posEvents = posEvents,
    negEvents = negEvents,
    diffEvents = diffEvents,
    everyTotNullEvent = everyTotNullEvent,
    allNullEvents = allNullEvents,
    posNullEvents = posNullEvents,
    negNullEvents = negNullEvents,
    recordCounts = recordCounts,
    allNullQuants = allNullQuants,
    posNullQuants = posNullQuants,
    negNullQuants = negNullQuants, 
    diffNullQuants = diffNullQuants)
  
  save(output, file = file.path(dataDir, paste0('results_',climateVar,'_EX.RData')))
  
  return(output)
}