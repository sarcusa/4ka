## Written by Hannah Kolus, 07/06/2018 
## Map the results of excursion analysis on a global grid. Test the spatial significance of events 
## with a null model (analysis performed on an ensemble of synthetic data records with the same 
## autocorrelation, trend, and temporal resolution as the orginal records)

library(lipdR)
library(ggplot2)
library(ggmap)
library(geosphere)
library(gridExtra)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('createPaths.R')

#eventYrs = seq(1000,10000,by = 200)
eventYrs = seq(1000,10600,by = 400)
msDir = file.path(createPaths(), 'mean_shift')
dir.create(file.path(msDir, 'spatial_MS_M'))
dir.create(file.path(msDir, 'spatial_MS_T'))

# parameters
res = 5            # grid resolution, deg
radius = 2000      # search radius, km
eventYr = 8200     # event year, yr BP
eventWindow = 199  # total event window = eventYr +/- eventWindow
numIt = 1000       # number of iterations in null model
climateVar = 'M'   # M for moisture, T for temperature

# Load in relevant results: TS_MS
load(file.path(createPaths(), 'RData', 'MS_results_plusNull_complete.RData'))
TS_MS_orig = TS_MS

for (y in 1:length(eventYrs)) {
  
  # reset the TS_MS dataset for the next event year
  TS_MS = TS_MS_orig
  
  eventYr = eventYrs[y]
  print(eventYr)

for (i in 1:length(TS_MS)) {
  
  # Filter records that don't contain event year in age range
  if (min(TS_MS[[i]]$age) > eventYr || max(TS_MS[[i]]$age) < eventYr) {
    TS_MS[[i]]$useMS = 0
    #print('Out of range')
  }
  
  # Only include temp12k records for temperature analysis
  if (climateVar == 'T') {
    # if the flag is defined, but not temp12k, exclude
    if (length(TS_MS[[i]]$paleoData_inCompilation) > 0) {
      if (tolower(TS_MS[[i]]$paleoData_inCompilation) != 'temp12k') {
        TS_MS[[i]]$useMS = -1
      }
    } else {
      # if the flag isn't defined, exclude
      TS_MS[[i]]$useMS = -1
    }
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
  if (!is.na(TS_MS[[i]]$sig_brks) && any(TS_MS[[i]]$sig_brks >= eventYr - eventWindow & TS_MS[[i]]$sig_brks <= eventYr + eventWindow)) {
    TS_MS[[i]]$eventMS = 1
    eve_i = which(TS_MS[[i]]$sig_brks >= eventYr - eventWindow & TS_MS[[i]]$sig_brks <= eventYr + eventWindow)
    TS_MS[[i]]$dirMS = TS_MS[[i]]$brk_dirs[eve_i]
  }

}

# isolate only the records corresponding to the chosen climate interpretation
interps = unlist(sapply(TS_MS,"[[","interpretation1_variable"))
if (climateVar == 'M') {
  inds = which(interps == 'M' | interps == 'P')
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
dirEvents = unlist(sapply(TS_MS,"[[","dirMS"))[inds]  # (0, 1, -1): (no, positive, negative) event
dirChange = dirs * dirEvents                          # (0, 1, -1): (no, positive, negative) climate event

## Spatial null

# Make a grid of lat/lon center points
latitude = seq(-90 + res/2, 90 - res/2, res)
longitude = seq(-180 + res/2, 180 - res/2, res)

# Grids that will contain fraction of passes
gridEvents = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNumRecs = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNullEvents = array(NA, dim = c(length(longitude), length(latitude), numIt))
totNullEvents = matrix(0, nrow = length(inds), ncol = numIt)
gridPercentEvents_null = matrix(NA, nrow = length(longitude), ncol = length(latitude))

# Assign event occurrence
for (i in 1:length(inds)) {
  
  nullBreaks = TS_MS[[inds[i]]]$null_sig_brks
  nullDirs = TS_MS[[inds[i]]]$null_brk_dirs
  
  # if no breaks are found in any iteration, nullBreaks is just an empty list and will break code
  if (length(nullBreaks) > 0) {
    for (j in 1:numIt) {
      
      eventInd = which(nullBreaks[[j]] >= eventYr - eventWindow & nullBreaks[[j]] <= eventYr + eventWindow)
      
      if (length(eventInd) > 0) {
        totNullEvents[i,j] = nullDirs[[j]][eventInd[1]] * dirs[i]
      }
      
    }  # end it loop (synth data)
  } 
}  # end event occurrence assignment

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

#save(gridEvents, gridNullEvents, allTsLat, allTsLon, dirChange, file = 'MS_3.8event.RData')

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
#percentEvents_NULL_binned = factor(percentEvents_NULL_binned,levels(percentEvents_NULL_binned)[c(6,9,8,7,5,2,3,4,1)])

if (climateVar == 'M') {
  ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = percentEvents_NULL)) + 
    borders("world", colour="black") + 
    scale_fill_distiller(name = 'Fraction\nof events', palette = 'RdBu', direction = 1) + 
    geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 2) +
    geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white', size = 2) +
    geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white', size = 2) +
    geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 1) +
    geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'), size = 1) +
    geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'), size = 1) +
    scale_color_manual(name = '', values = c('no event' = 'grey50', '+ event' = 'blue', '- event' = 'red')) +
    theme_bw() + xlab('Longitude') + ylab('Latitude') +
    xlim(-180, 180) + ylim(-90, 90) +
    ggtitle(paste0('MS: Fraction of null events < real event #\n', eventYr/1000,'+/-',eventWindow/2/1000, 'ka events'))
  
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
    ggtitle(paste0(eventYr/1000,' ka mean shift: Moisture\nProbability that event density is simulated in the null'))
  #ggsave(paste0('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/MS_M_plots/',y,'-', eventYr/1000, '.pdf'))
  ggsave(file.path(msDir, 'spatial_MS_M', paste0(y,'-', eventYr/1000, '.pdf')))
}

if (climateVar == 'T') {
  ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = percentEvents_NULL)) + 
    borders("world", colour="black") + 
    scale_fill_distiller(name = 'Fraction\nof events', palette = 'RdBu', direction = -1) + 
    geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 2) +
    geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white', size = 2) +
    geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white', size = 2) +
    geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 1) +
    geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'), size = 1) +
    geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'), size = 1) +
    scale_color_manual(name = '', values = c('no event' = 'grey50', '+ event' = 'red', '- event' = 'blue')) +
    theme_bw() + xlab('Longitude') + ylab('Latitude') +
    xlim(-180, 180) + ylim(-90, 90) +
    ggtitle(paste0('MS: Fraction of null events < real event #\n', eventYr/1000,'+/-',eventWindow/2/1000, 'ka events'))
  
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
    ggtitle(paste0(eventYr/1000,' ka mean shift: Temperature\nProbability that event density is simulated in the null'))
  #ggsave(paste0('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/MS_T_plots/',y,'-', eventYr/1000, '.pdf'))
  ggsave(file.path(msDir, 'spatial_MS_T', paste0(y,'-', eventYr/1000, '.pdf')))
}

} # end loop thru event years
