spatialMeanShift_null <- function(data_in, param, climateVar){
  
  msDir = file.path(createPaths(), 'mean_shift')
  dir.create(file.path(msDir, 'spatial_MS_M'))
  dir.create(file.path(msDir, 'spatial_MS_T'))
  
  TS_MS_orig = data_in
  
  # Manual
  load("/projects/pd_lab/sha59/4ka/RData/MS_results_plusNull_complete.RData")
  TS_MS_orig = analysis_2b
  
  for (y in 1:length(param$eventYrs)) {
    
    # reset the TS_MS dataset for the next event year
    TS_MS = TS_MS_orig
    
    eventYr = param$eventYrs[y]
    print(eventYr)
    
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
      if (!is.na(TS_MS[[i]]$sig_brks) && any(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow)) {
        TS_MS[[i]]$eventMS = 1
        eve_i = which(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow)
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
          gridEvents[i, j] = ifelse(length(indsLoc) > 1, sum(dirChange[indsLoc] == 1) - sum(dirChange[indsLoc] == -1), NA)
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
    
    # pull out non NA values
    locs = which(!is.na(gridEvents), arr.ind = T)
    locs_na = which(is.na(gridEvents), arr.ind = T)
    percentEvents_NULL = gridPercentEvents_null[locs]
    
    bins = c('<= 0.05', '0.05 - 0.1', '0.1 - 0.25', '0.25 - 0.5', '> 0.5', '0.5 - 0.25', '0.25 - 0.1', '0.1 - 0.05', ' <= 0.05')
    locs_binned = matrix(nrow = nrow(locs) + length(bins), ncol = 2)
    locs_binned[1:nrow(locs),] = locs
    percentEvents_NULL_binned = percentEvents_NULL
    percentEvents_NULL_binned[percentEvents_NULL >= 0.95] = bins[1]
    percentEvents_NULL_binned[percentEvents_NULL < 0.95 & percentEvents_NULL >= 0.9] = bins[2]
    percentEvents_NULL_binned[percentEvents_NULL < 0.9 & percentEvents_NULL >= 0.75] = bins[3]
    percentEvents_NULL_binned[percentEvents_NULL < 0.75 & percentEvents_NULL >= 0.5] = bins[4]
    percentEvents_NULL_binned[percentEvents_NULL < 0.5 & percentEvents_NULL > -0.5] = bins[5]
    percentEvents_NULL_binned[percentEvents_NULL > -0.75 & percentEvents_NULL <= -0.50] = bins[6]
    percentEvents_NULL_binned[percentEvents_NULL > -0.9 & percentEvents_NULL <= -0.75] = bins[7]
    percentEvents_NULL_binned[percentEvents_NULL > -0.95 & percentEvents_NULL <= -0.9] = bins[8]
    percentEvents_NULL_binned[percentEvents_NULL <= -0.95] = bins[9]
    
    dummy = length(percentEvents_NULL_binned) + 1
    dummy_i = dummy:(dummy+length(bins)-1)
    for (i in 1:length(bins)) {
      percentEvents_NULL_binned[dummy] = bins[i]
      locs_binned[dummy,] = locs_na[10*i,]
      dummy = dummy + 1
    }
    
    percentEvents_NULL_binned = factor(percentEvents_NULL_binned)
    percentEvents_NULL_binned = factor(percentEvents_NULL_binned,levels(percentEvents_NULL_binned)[c(2,4,6,8,3,9,7,5,1)])
        
    if (climateVar == 'M') {
      
      myCol = c('#543005','#bf812d','#dfc27d','#f6e8c3','snow2','#c7eae5','#80cdc1','#35978f','#003c30')
      
      tryCatch({s <-baseMAP(lon=0,lat = 0,projection = "mollweide",global = TRUE,map.type = "line",restrict.map.range = F, country.boundaries = F) + 
        geom_tile(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]],width = 5.5,height = 5.5, fill = as.factor(percentEvents_NULL_binned))) + 
        borders("world", colour="grey70") + 
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 2, shape = 21, fill = "grey50") +
        geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'), size = 3, shape = 24, fill = "blue") +
        geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'), size = 3, shape = 25, fill = "tomato4") +
        scale_color_manual(name = '', values = c('no event' = 'black', '+ event' = 'white', '- event' = 'white'),breaks = c('+ event', '- event', 'no event'),guide = guide_legend(override.aes = list(shape = c(24, 25, 21), fill = c('blue','tomato4','grey50'),color = c('black','black','black')))) +
        scale_fill_manual(name = '', values = rev(myCol))+
        ggtitle(paste0('MS: Fraction of null events < real event #\n', eventYr/1000,'+/-',param$eventWindow/2/1000, 'ka events'))+
        #geom_rect(aes(xmax=180.1,xmin=-180.1,ymax=90.1,ymin=-90.1),fill=NA, colour="black")     
      
      pdf(file.path(msDir, 'spatial_MS_M', 
                    paste0(y,'-', eventYr/1000, '.pdf')))
      print(s)
      dev.off() }, error = function(e){cat("Error:", conditionMessage(e), " may not have plotted")})   
      
    }
    
    if (climateVar == 'T') {
      
      myCol = c('#67001f','#d6604d','#f4a582','#fddbc7','snow2','#d1e5f0','#92c5de','#4393c3','#053061')
      
      tryCatch({ s <- baseMAP(lon=0,lat = 0,projection = "mollweide",global = TRUE,map.type = "line",restrict.map.range = F, country.boundaries = F) + 
        geom_tile(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]],width = 5.5,height = 5.5, fill = as.factor(percentEvents_NULL_binned))) + 
        borders("world", colour="grey70") + 
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 2, fill = "grey50", shape = 21) +
        geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'), size = 3, fill = "red", shape = 24) +
        geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'), size = 3, fill = "royalblue", shape = 25) +
        scale_color_manual(name = '', values = c('no event' = 'black', '+ event' = 'white', '- event' = 'white'),
                           breaks = c('+ event', '- event', 'no event'),
                           guide = guide_legend(override.aes = list(shape = c(24, 25, 21), fill = c('red','royalblue','grey50'),color = c('black','black','black')))) +
        scale_fill_manual(name = '', values = myCol) +
        ggtitle(paste0('MS: Fraction of null events < real event #\n', eventYr/1000,'+/-',param$eventWindow/2/1000, 'ka events'))+
        #geom_rect(aes(xmax=180.1,xmin=-180.1,ymax=90.1,ymin=-90.1),fill=NA, colour="black")
      
      pdf(file.path(msDir, 'spatial_MS_T', 
                       paste0(y,'-', eventYr/1000, '.pdf')))
      print(s)
      dev.off()}, error = function(e){cat("Error:", conditionMessage(e), " may not have plotted")})
      
      
      
    }
    
  } 
  
  
}