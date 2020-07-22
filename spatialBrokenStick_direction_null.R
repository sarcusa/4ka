spatialBrokenStick <- function(data_in,param, climateVar){
  
  bsDir = file.path(createPaths(), 'broken_stick')
  dir.create(file.path(bsDir, 'spatial_BS_M'))
  dir.create(file.path(bsDir, 'spatial_BS_T'))
  
  eventYrs = param$eventYrs
  
  TS_BS_orig = data_in[[1]]
  
  for (y in 1:length(eventYrs)) {
    
    TS_BS = TS_BS_orig
    
    eventYr = eventYrs[y]
    print(eventYr)
    
    for (i in 1:length(TS_BS)) {
      
      # Filter records that don't contain event year in age range
      if (min(TS_BS[[i]]$age) > eventYr || max(TS_BS[[i]]$age) < eventYr) {
        TS_BS[[i]]$useBS = 0
        print('Out of range')
      }
      
      
      #Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
      if (length(TS_BS[[i]]$interpretation1_seasonalityGeneral) > 0) {
        if (tolower(TS_BS[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
              tolower(TS_BS[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
          TS_BS[[i]]$useBS = -1
          print(TS_BS[[i]]$interpretation1_seasonalityGeneral)
        }
      }
      
    }
    TS_BS = filterTs(TS_BS, 'useBS == 1')
    
    # Assign event occurrence
    for (i in 1:length(TS_BS)) {
      
      TS_BS[[i]]$eventBS = 0
      TS_BS[[i]]$dirBS = 0
      if (!is.na(TS_BS[[i]]$brk_pts) && any(TS_BS[[i]]$brk_pts >= eventYr - param$eventWindow & TS_BS[[i]]$brk_pts <= eventYr + param$eventWindow)) {
        eve_i = which(TS_BS[[i]]$brk_pts >= eventYr - param$eventWindow & TS_BS[[i]]$brk_pts <= eventYr + param$eventWindow)
        
        # Metric combining events occurring in event window (event takes the sign of the sum)
        dir = sum(TS_BS[[i]]$brk_dirs[eve_i]) / abs(sum(TS_BS[[i]]$brk_dirs[eve_i]))
        TS_BS[[i]]$eventBS = ifelse(is.na(dir), 0, 1)
        TS_BS[[i]]$dirBS = ifelse(is.na(dir), 0, dir)
      }
      
    }
    
    # isolate only the records corresponding to the chosen climate interpretation
    interps = unlist(sapply(TS_BS,"[[","interpretation1_variable"))
    if (climateVar == 'M') {
      inds = which(interps == 'M' | interps == 'P')
    } else {
      inds = which(interps == 'T' | interps == 'TM')
    }
    
    allTsLat = as.numeric(sapply(TS_BS,"[[","geo_latitude"))[inds]
    allTsLon = as.numeric(sapply(TS_BS,"[[","geo_longitude"))[inds]
    events = as.numeric(sapply(TS_BS,"[[","eventBS"))[inds]
    interps = interps[inds]
    
    # calculate the climate event direction based on the proxy climate dir and event dir
    dirs = unlist(sapply(TS_BS,"[[","interpretation1_interpDirection"))[inds]
    dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
    dirs[dirs == 'negative'] = -1
    dirs[dirs == 'NA' | is.na(dirs)] = 0
    dirs = as.numeric(dirs)
    dirEvents = unlist(sapply(TS_BS,"[[","dirBS"))[inds]  # (0, 1, -1): (no, positive, negative) event
    dirChange = dirs * dirEvents  # (0, 1, -1): (no, positive, negative) climate event
    
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
      
      nullBreaks = TS_BS[[inds[i]]]$null_brk_pts
      nullDirs = TS_BS[[inds[i]]]$null_brk_dirs
      
      for (j in 1:param$numIt) {
        if(length(nullBreaks) < j){
          stop("uh oh, j is too big")
        }
       
        eventInd = which(nullBreaks[[j]] >= eventYr - param$eventWindow & nullBreaks[[j]] <= eventYr + param$eventWindow)
        
        if (length(eventInd) > 0) {
          
          # Metric combining events occurring in event window (event takes the sign of the sum)
          dir = sum(nullDirs[[j]][eventInd]) / abs(sum(nullDirs[[j]][eventInd]))
          totNullEvents[i,j] = ifelse(is.na(dir), 0, dir) * dirs[i]
          
        }
        
      }
      
    }
    
    for (i in 1:length(longitude)) {
      
      print(paste('iteration i = ', i))
      
      for (j in 1:length(latitude)) {
        
        indsLoc = which(distGeo(cbind(allTsLon, allTsLat), 
                                c(longitude[i],latitude[j])) <= param$radius*1000)
        
        if (length(indsLoc) > 0) { 
          # number of events = # pos events - # neg events, must have > 1 record contributing to grid cell
          gridNumRecs[i, j] = length(indsLoc)
          gridEvents[i, j] = ifelse(length(indsLoc) > 1, 
                                    sum(dirChange[indsLoc] == 1) - 
                                      sum(dirChange[indsLoc] == -1), NA)
        }
        
        if (length(indsLoc) == 1) {
          gridNullEvents[i, j, ] = totNullEvents[indsLoc, ]
        } else {
          gridNullEvents[i, j, ] = apply(totNullEvents[indsLoc, ], 2, 
                                         function (x) sum(x == 1) - sum(x == -1))
        }
        
        if (length(indsLoc) > 0) {
          if (is.na(gridEvents[i, j])) {
            next
          } else if (gridEvents[i, j] < 0) {
            gridPercentEvents_null[i,j] = -sum(gridNullEvents[i,j,] > 
                                                 gridEvents[i, j]) / param$numIt
          } else if (gridEvents[i, j] > 0) {
            gridPercentEvents_null[i,j] = sum(gridNullEvents[i,j,] < 
                                                gridEvents[i, j]) / param$numIt
          } else {
            gridPercentEvents_null[i,j] = 0
          }
        }
      }
      
    }
    
    
    locs = which(!is.na(gridEvents), arr.ind = T)
    locs_na = which(is.na(gridEvents), arr.ind = T)
    percentEvents_NULL = gridPercentEvents_null[locs]
    
    bins = c('<= 0.05', '0.05 - 0.1', '0.1 - 0.25', '0.25 - 0.5', '> 0.5', '0.5 - 0.25', '0.25 - 0.1', '0.1 - 0.05', ' <= 0.05')
    locs_binned = matrix(nrow = nrow(locs) + length(bins), ncol = 2)
    locs_binned[1:nrow(locs),] = locs
    percentEvents_NULL_binned = percentEvents_NULL
    percentEvents_NULL_binned[percentEvents_NULL >= 0.95] = bins[1]
    percentEvents_NULL_binned[percentEvents_NULL < 0.95 & 
                                percentEvents_NULL >= 0.9] = bins[2]
    percentEvents_NULL_binned[percentEvents_NULL < 0.9 & 
                                percentEvents_NULL >= 0.75] = bins[3]
    percentEvents_NULL_binned[percentEvents_NULL < 0.75 & 
                                percentEvents_NULL >= 0.5] = bins[4]
    percentEvents_NULL_binned[percentEvents_NULL < 0.5 & 
                                percentEvents_NULL > -0.5] = bins[5]
    percentEvents_NULL_binned[percentEvents_NULL > -0.75 & 
                                percentEvents_NULL <= -0.50] = bins[6]
    percentEvents_NULL_binned[percentEvents_NULL > -0.9 & 
                                percentEvents_NULL <= -0.75] = bins[7]
    percentEvents_NULL_binned[percentEvents_NULL > -0.95 & 
                                percentEvents_NULL <= -0.9] = bins[8]
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
          
      ## BINNED, more distinct color for > 95%
      myCol = c('#543005','#bf812d','#dfc27d','#f6e8c3','snow2','#c7eae5','#80cdc1','#35978f','#003c30')
     tryCatch({ ggplot() + geom_raster(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]], fill = as.factor(percentEvents_NULL_binned))) + 
        geom_tile(aes(x = longitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),1]],
                      y = latitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),2]]), 
                  height=5, width=6, fill = 'white') +
        borders("world", colour="black") + 
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'),fill='grey50',shape = 21, size = 2) +
        geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='wet event'),fill='blue',shape = 24, size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='dry event'),fill='tomato4',shape = 25, size = 3) +
        scale_color_manual(name = '', values = c('no event' = 'black', 'wet event' = 'white', 'dry event' = 'white'),
                           breaks = c('wet event', 'dry event', 'no event'),
                           guide = guide_legend(override.aes = list(shape = c(24, 25, 21), fill = c('blue','tomato4','grey50'),
                                                                    color = c('black','black','black')))) +
        scale_fill_manual(name = '', values = rev(myCol)) +
        theme_bw() + xlab('Longitude') + ylab('Latitude') +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(-180, 180) + ylim(-90, 90) +
        ggtitle(paste0(eventYr/1000,' ka broken stick: Moisture\nProbability that event density is simulated in the null'))
      ggsave(file.path(bsDir, 'spatial_BS_M', paste0(y,'-', eventYr/1000, '.pdf')))
      
    }, error = function(e){cat("Error:", conditionMessage(e), " may not have plotted")})
      
    }
    
    if (climateVar == 'T') {
      
     ## BINNED, more distinct color for > 95%
      myCol = c('#67001f','#d6604d','#f4a582','#fddbc7','snow2','#d1e5f0','#92c5de','#4393c3','#053061')
     tryCatch({ ggplot() + geom_raster(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]], fill = as.factor(percentEvents_NULL_binned))) + 
        geom_tile(aes(x = longitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),1]],
                      y = latitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),2]]), 
                  height=5, width=6, fill = 'white') +
        borders("world", colour="black") + 
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), fill='grey50',shape = 21,size = 2) +
        geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='warm event'),fill='red',shape = 24,size = 3) +
        geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='cold event'),fill='royalblue',shape = 25, size = 3) +
        scale_color_manual(name = '', values = c('no event' = 'black', 'warm event' = 'white', 'cold event' = 'white'),
                           breaks = c('warm event', 'cold event', 'no event'),
                           guide = guide_legend(override.aes = list(shape = c(24, 25, 21), fill = c('red','royalblue','grey50'),
                                                                    color = c('black','black','black')))) +
        scale_fill_manual(name = '', values = myCol) +
        theme_bw() + xlab('Longitude') + ylab('Latitude') +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(-180, 180) + ylim(-90, 90) +
        ggtitle(paste0(eventYr/1000,' ka broken stick: Temperature\nProbability that event density is simulated in the null'))
      ggsave(file.path(bsDir, 'spatial_BS_T', paste0(y,'-', eventYr/1000, '.pdf')))
       
     }, error = function(e){cat("Error:", conditionMessage(e), " may not have plotted")})
      
     }
    
  }
  
}