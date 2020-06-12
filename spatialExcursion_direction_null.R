## Written by Hannah Kolus, 07/06/2018 
## Map the results of excursion analysis on a global grid. Test the spatial significance of events 
## with a null model (analysis performed on an ensemble of synthetic data records with the same 
## autocorrelation, trend, and temporal resolution as the orginal records)
  
  library(lipdR)
  library(ggplot2)
  library(ggmap)
  #library(ggalt)
  library(geosphere)
  library(gridExtra)
  setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
  source('createPaths.R')
  
  # parameters
  res = 5            # grid resolution, deg
  radius = 2000      # search radius, km
  eventYr = 8200     # event year, yr BP
  numIt = 1000       # number of iterations in null model
  climateVar = 'M'   # M for moisture, T for temperature
  
  eventYrs = seq(1000,10600,by = 400)
  exDir = file.path(createPaths(), 'excursion')
  dir.create(file.path(exDir, 'spatial_EX_M'))
  dir.create(file.path(exDir, 'spatial_EX_T'))
  
  for (y in 1:length(eventYrs)) {
    
    eventYr = eventYrs[y]
  
  # Load in relevant results: TS_EX
  #load('EX_results_plusNull_complete_4.2_short.RData')
  load(file.path(createPaths(), 'RData', paste0('EX_results_plusNull_complete_', eventYr/1000, '.RData')))
  
  for (i in 1:length(TS_EX)) {
    # Only include temp12k records for temperature analysis
    if (climateVar == 'T') {
      # if the flag is defined, but not temp12k, exclude
      if (length(TS_EX[[i]]$paleoData_inCompilation) > 0) {
        if (tolower(TS_EX[[i]]$paleoData_inCompilation) != 'temp12k') {
          TS_EX[[i]]$useEX = -1
        }
      } else {
        # if the flag isn't defined, exclude
        TS_EX[[i]]$useEX = -1
      }
    }
    
    # Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
    if (length(TS_EX[[i]]$interpretation1_seasonalityGeneral) > 0) {
      if (tolower(TS_EX[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
          tolower(TS_EX[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
        TS_EX[[i]]$useEX = -1
        print(TS_EX[[i]]$interpretation1_seasonalityGeneral)
      }
    }
  }
  TS_EX = filterTs(TS_EX, 'useEX == 1')
  
  # isolate only the records corresponding to the chosen climate interpretation
  interps = unlist(sapply(TS_EX,"[[","interpretation1_variable"))
  if (climateVar == 'M') {
    inds = which(interps == 'M' | interps == 'P')
  } else {
    inds = which(interps == 'T' | interps == 'TM')
  }
  
  allTsLat = as.numeric(sapply(TS_EX,"[[","geo_latitude"))[inds]
  allTsLon = as.numeric(sapply(TS_EX,"[[","geo_longitude"))[inds]
  events = as.numeric(sapply(TS_EX,"[[","eventEX"))[inds]
  interps = interps[inds]
  
  # calculate the climate event direction based on the proxy climate dir and event dir
  dirs = unlist(sapply(TS_EX,"[[","interpretation1_interpDirection"))[inds]
  dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
  dirs[dirs == 'negative'] = -1
  dirs[dirs == 'NA' | is.na(dirs)] = 0
  dirs = as.numeric(dirs)
  dirEvents = unlist(sapply(TS_EX,"[[","dirEx"))[inds]  # (0, 1, -1): (no, positive, negative) event
  dirChange = dirs * dirEvents                          # (0, 1, -1): (no, positive, negative) climate event
  
  # Array that will contain fraction of passes 
  gridNullEvents = matrix(NA, nrow = length(inds), ncol = 1)
  
  # Calculate fraction of null model not a climate event in the real event's direction
  for (i in 1:length(inds)) {
    
    nullEvents = TS_EX[[inds[i]]]$null_events_dir * dirs[i]
    gridNullEvents[i] = dirChange[i] * sum(nullEvents != dirChange[i]) / numIt
    
  }
  
  #i_pos = which(dirChange == 1)
  #i_neg = which(dirChange == -1)
  #i_no = which(dirChange == 0)
  
  # bins = c('>= 95%','-(>= 95%)','90 - 95%','-(90 - 95%)','75 - 90%','-(75 - 90%)','50 - 75%',
  #          '-(50 - 75%)','< 50%')
  # 
  # gridNullEvents_binned = gridNullEvents
  # gridNullEvents_binned[gridNullEvents >= 0.95] = bins[1]
  # gridNullEvents_binned[gridNullEvents <= -0.95] = bins[2]
  # gridNullEvents_binned[gridNullEvents < 0.95 & gridNullEvents >= 0.9] = bins[3]
  # gridNullEvents_binned[gridNullEvents > -0.95 & gridNullEvents <= -0.9] = bins[4]
  # gridNullEvents_binned[gridNullEvents < 0.9 & gridNullEvents >= 0.75] = bins[5]
  # gridNullEvents_binned[gridNullEvents > -0.9 & gridNullEvents <= -0.75] = bins[6]
  # gridNullEvents_binned[gridNullEvents < 0.75 & gridNullEvents >= 0.5] = bins[7]
  # gridNullEvents_binned[gridNullEvents > -0.75 & gridNullEvents <= -0.5] = bins[8]
  # gridNullEvents_binned[gridNullEvents < 0.5 & gridNullEvents > -0.5] = bins[9]
  # 
  # dummy = length(gridNullEvents_binned) + 1
  # dummy_i = dummy:(dummy+length(bins)-1)
  # for (i in 1:length(bins)) {
  #   gridNullEvents_binned[dummy] = bins[i]
  #   allTsLon[dummy] = allTsLon[i_pos[1]]
  #   allTsLat[dummy] = allTsLat[i_pos[1]]
  #   dummy = dummy + 1
  # }
  # 
  # gridNullEvents_binned = factor(gridNullEvents_binned)
  # gridNullEvents_binned = factor(gridNullEvents_binned,levels(gridNullEvents_binned)[c(6,9,8,7,5,2,3,4,1)])
  # 
  # if (climateVar == 'M') {
  #   
  #   ggplot() + borders("world", colour="black") +
  #     geom_point(aes(x = allTsLon[i_no], y = allTsLat[i_no], color = 'no event'), shape = 19, size = 3) +
  #     geom_point(aes(x = allTsLon[i_pos], y = allTsLat[i_pos], fill = gridNullEvents[i_pos]), shape = 24, size = 3) +
  #     geom_point(aes(x = allTsLon[i_neg], y = allTsLat[i_neg], fill = gridNullEvents[i_neg]), shape = 25, size = 3) +
  #     scale_fill_gradient2(name = 'Event significance', low = 'saddlebrown', mid = 'white', high = 'forestgreen', midpoint = 0) +
  #     scale_color_manual(name = '', values = c('no event' = 'grey50')) +
  #     theme_bw() + xlab('Longitude') + ylab('Latitude') +
  #     xlim(-180, 180) + ylim(-90, 90) 
  #   
  #   ggplot() + borders("world", colour="black") +
  #     geom_point(aes(x = allTsLon[i_no], y = allTsLat[i_no]), color = 'black', shape = 19, size = 4) +
  #     geom_point(aes(x = allTsLon[i_no], y = allTsLat[i_no], color = 'no event'), shape = 19, size = 3) +
  #     geom_point(aes(x = allTsLon[dummy_i], y = allTsLat[dummy_i], fill = gridNullEvents_binned[dummy_i]), shape = 24, size = 3) +
  #     geom_point(aes(x = allTsLon[i_pos], y = allTsLat[i_pos], fill = gridNullEvents_binned[i_pos]), shape = 24, size = 3) +
  #     geom_point(aes(x = allTsLon[i_neg], y = allTsLat[i_neg], fill = gridNullEvents_binned[i_neg]), shape = 25, size = 3) +
  #     scale_fill_brewer(name = 'Event significance', palette = 'RdBu', direction = -1, guide = guide_legend(override.aes = list(
  #       shape = c(rep(24, 4), 5, rep(25, 4))))) +
  #     scale_color_manual(name = '', values = c('no event' = 'grey50')) +
  #     theme_bw() + xlab('Longitude') + ylab('Latitude') +
  #     xlim(-180, 180) + ylim(-90, 90) +
  #     ggtitle(paste0(eventYr/1000, ' excursion: Moisture'))
  #   
  # }
  # 
  # if (climateVar == 'T') {
  #   
  #  ggplot() + borders("world", colour="black") +
  #     geom_point(aes(x = allTsLon[i_no], y = allTsLat[i_no], fill = gridNullEvents[i_no]), shape = 19, size = 3) +
  #     geom_point(aes(x = allTsLon[i_pos], y = allTsLat[i_pos], fill = gridNullEvents[i_pos]), shape = 24, size = 3) +
  #     geom_point(aes(x = allTsLon[i_neg], y = allTsLat[i_neg], fill = gridNullEvents[i_neg]), shape = 25, size = 3) +
  #     scale_fill_gradient2(name = '', low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  #     theme_bw() + xlab('Longitude') + ylab('Latitude') +
  #     xlim(-180, 180) + ylim(-90, 90) 
  #   
  #   ggplot() + borders("world", colour="black") +
  #     geom_point(aes(x = allTsLon[i_no], y = allTsLat[i_no]), color = 'black', shape = 19, size = 4) +
  #     geom_point(aes(x = allTsLon[i_no], y = allTsLat[i_no], color = 'no event'), shape = 19, size = 3) +
  #     geom_point(aes(x = allTsLon[dummy_i], y = allTsLat[dummy_i], fill = gridNullEvents_binned[dummy_i]), shape = 24, size = 3) +
  #     geom_point(aes(x = allTsLon[i_pos], y = allTsLat[i_pos], fill = gridNullEvents_binned[i_pos]), shape = 24, size = 3) +
  #     geom_point(aes(x = allTsLon[i_neg], y = allTsLat[i_neg], fill = gridNullEvents_binned[i_neg]), shape = 25, size = 3) +
  #     scale_fill_brewer(name = 'Event significance', palette = 'RdBu', guide = guide_legend(override.aes = list(
  #       shape = c(rep(24, 4), 5, rep(25, 4))))) +
  #     scale_color_manual(name = '', values = c('no event' = 'grey50')) +
  #     theme_bw() + xlab('Longitude') + ylab('Latitude') +
  #     xlim(-180, 180) + ylim(-90, 90) +
  #     ggtitle(paste0(eventYr/1000, ' excursion: Temperature'))
  #   
  # }
  
  ## Spatial null
  
  # Make a grid of lat/lon center points
  latitude = seq(-90 + res/2, 90 - res/2, res)
  longitude = seq(-180 + res/2, 180 - res/2, res)
  
  # Grids that will contain fraction of passes
  gridEvents = matrix(NA, nrow = length(longitude), ncol = length(latitude))
  gridNumRecs = matrix(NA, nrow = length(longitude), ncol = length(latitude))
  gridNullEvents = array(NA, dim = c(length(longitude), length(latitude), numIt))
  totNullEvents = matrix(NA, nrow = length(inds), ncol = numIt)
  gridPercentEvents_null = matrix(NA, nrow = length(longitude), ncol = length(latitude))
  
  # Assign event occurrence
  for (i in 1:length(inds)) {
    
    totNullEvents[i,] = TS_EX[[inds[i]]]$null_events_dir * dirs[i]
    
  }
  
  for (i in 1:length(longitude)) {
    
    print(paste('iteration i = ', i))
    
    for (j in 1:length(latitude)) {
      
      indsLoc = which(distGeo(cbind(allTsLon, allTsLat), c(longitude[i],latitude[j])) <= radius*1000)
      
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
          gridPercentEvents_null[i,j] = -sum(gridNullEvents[i,j,] > gridEvents[i, j]) / numIt
        } else if (gridEvents[i, j] > 0) {
          gridPercentEvents_null[i,j] = sum(gridNullEvents[i,j,] < gridEvents[i, j]) / numIt
        } else {
          gridPercentEvents_null[i,j] = 0
        }
      }
    }
    
  } # end location grid loops
  
  #save(gridEvents, gridNullEvents, allTsLat, allTsLon, dirChange, file = 'EX_3.8event.RData')
  
  # pull out non NA values
  locs = which(!is.na(gridEvents), arr.ind = T)
  locs_na = which(is.na(gridEvents), arr.ind = T)
  percentEvents_NULL = gridPercentEvents_null[locs]
  
  #bins = c('>= 95%','90 - 95%','75 - 90%','50 - 75%','< 50%','-(50 - 75%)','-(75 - 90%)','-(90 - 95%)','-(>= 95%)')
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
    
    # ggplot() + geom_raster(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]], fill = as.factor(percentEvents_NULL_binned))) + 
    #   geom_tile(aes(x = longitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),1]],
    #                 y = latitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),2]]), 
    #             height=5, width=6, fill = 'white') +
    #   borders("world", colour="black") + 
    #   geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
    #   geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white',shape = 24, size = 2) +
    #   geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white',shape = 25, size = 2) +
    #   geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 2) +
    #   geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'),fill='blue',shape = 24, size = 2, stroke=1.3) +
    #   geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'),fill='tomato4',shape = 25, size = 2, stroke=1.3) +
    #   scale_color_manual(name = '', values = c('no event' = 'grey50', '+ event' = 'blue', '- event' = 'tomato4'),
    #                      guide = guide_legend(override.aes = list(shape = c(25, 24, 19), fill = c('tomato4','blue','grey50'),
    #                                                               color = c('black','black','black')))) +
    #   scale_fill_manual(name = '', values = rev(myCol)) +
    #   theme_bw() + xlab('Longitude') + ylab('Latitude') +
    #   xlim(-180, 180) + ylim(-90, 90) +
    #   ggtitle(paste0(eventYr/1000,' ka excursion:\nPercent of null events < real event'))
    
    ## BINNED, more distinct color for > 95%
    myCol = c('#543005','#bf812d','#dfc27d','#f6e8c3','snow2','#c7eae5','#80cdc1','#35978f','#003c30')
    ggplot() + geom_raster(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]], fill = as.factor(percentEvents_NULL_binned))) + 
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
      ggtitle(paste0(eventYr/1000,' ka excursion: Moisture\nProbability that event density is simulated in the null'))
    ggsave(file.path(exDir, 'spatial_EX_M', paste0(y,'-', eventYr/1000, '.pdf')))
    
  }
  
  if (climateVar == 'T') {
    
    # ggplot() + geom_raster(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]], fill = as.factor(percentEvents_NULL_binned))) + 
    #   geom_tile(aes(x = longitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),1]],
    #                 y = latitude[locs_binned[(length(percentEvents_NULL)+1):length(percentEvents_NULL_binned),2]]), 
    #             height=5, width=6, fill = 'white') +
    #   borders("world", colour="black") + 
    #   geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
    #   geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white',shape = 24, size = 3) +
    #   geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white',shape = 25, size = 3) +
    #   geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 2) +
    #   geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'),shape = 24,size = 2, stroke = 1.3) +
    #   geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'),shape = 25, size = 2, stroke = 1.3) +
    #   scale_color_manual(name = '', values = c('no event' = 'grey50', '+ event' = 'red', '- event' = 'blue'),
    #                      guide = guide_legend(override.aes = list(shape = c(25, 24, 19)))) +
    #   scale_fill_manual(name = '', values = myCol) +
    #   theme_bw() + xlab('Longitude') + ylab('Latitude') +
    #   xlim(-180, 180) + ylim(-90, 90) +
    #   ggtitle(paste0(eventYr/1000,' ka excursion:\nProbability that event density is simulated in the null'))
    
    ## BINNED, more distinct color for > 95%
    myCol = c('#67001f','#d6604d','#f4a582','#fddbc7','snow2','#d1e5f0','#92c5de','#4393c3','#053061')
    ggplot() + geom_raster(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]], fill = as.factor(percentEvents_NULL_binned))) + 
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
      ggtitle(paste0(eventYr/1000,' ka excursion: Temperature\nProbability that event density is simulated in the null'))
    #ggsave(paste0('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/EX_T_plots/4.05.pdf'))
    ggsave(file.path(exDir, 'spatial_EX_T', paste0(y,'-', eventYr/1000, '.pdf')))
  }
  
  }
  
  
