## Written by Hannah Kolus, 09/04/2018 
## Map the results of broken stick analysis on a global grid. Test the spatial significance of events 
## with a null model (analysis performed on an ensemble of synthetic data records with the same 
## autocorrelation, trend, and temporal resolution as the orginal records)

library(lipdR)
library(ggplot2)
library(ggmap)
library(geosphere)

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')

# Load in relevant results: TS_BS
load('RData/BS_results_plusNull.RData')

# parameters
res = 5            # grid resolution, deg
radius = 1500      # search radius, km
eventYr = 8200     # event year, yr BP
eventWindow = 500  # total event window = eventYr +/- eventWindow
numIt = 1000       # number of iterations to create the null model

# Make a grid of lat/lon center points
latitude = seq(-90 + res/2, 90 - res/2, res)
longitude = seq(-180 + res/2, 180 - res/2, res)

# Filter records that don't contain event year in age range
for (i in 1:length(TS_BS)) {
  if (min(TS_BS[[i]]$age) > eventYr || max(TS_BS[[i]]$age) < eventYr) {
    TS_BS[[i]]$useBS = 0
    print('Out of range')
  }
}
TS_BS = filterTs(TS_BS, 'useBS == 1')

# Grids that will contain fraction of passes
gridEvents = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNumRecs = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNullEvents = array(NA, dim = c(length(longitude), length(latitude), numIt))
totNullEvents = matrix(NA, nrow = length(TS_BS), ncol = numIt)

# Assign event occurrence
for (i in 1:length(TS_BS)) {

  TS_BS[[i]]$eventOcc = 0
  if (!is.na(TS_BS[[i]]$brk_pts) && any(TS_BS[[i]]$brk_pts >= eventYr - eventWindow & TS_BS[[i]]$brk_pts <= eventYr + eventWindow)) {
    TS_BS[[i]]$eventOcc = 1
  }
  
  nullBrks = TS_BS[[i]]$null_brk_pts
  totNullEvents[i,] = sapply(nullBrks, function(x) sum(x >= eventYr - eventWindow & x <= eventYr + eventWindow, na.rm = T))
}

allTsLat = as.numeric(sapply(TS_BS,"[[","geo_latitude"))
allTsLon = as.numeric(sapply(TS_BS,"[[","geo_longitude"))
events = sapply(TS_BS,"[[","eventOcc")

for (i in 1:length(longitude)) {
  
  print(paste('iteration i = ', i))
  
  for (j in 1:length(latitude)) {
    
    inds = which(distGeo(cbind(allTsLon, allTsLat), c(longitude[i],latitude[j])) <= radius*1000)
    
    if (length(inds) > 0) { 
      gridEvents[i, j] = sum(events[inds])
      gridNumRecs[i, j] = length(inds)
    }
    
    if (length(inds) == 1) {
      gridNullEvents[i, j, ] = totNullEvents[inds, ]
    } else {
      gridNullEvents[i, j, ] = apply(totNullEvents[inds, ], 2, sum)
    }
    
  }
  
}

# Calculate the fraction of sites within the radius that pass the test
gridEvents[gridNumRecs == 1] = NA
gridNullEvents_95 = apply(gridNullEvents, c(1,2), function (x) quantile(x, 0.95))
#gridNullEvents_95[gridNullEvents_95 == 0] = .Machine$double.xmin # ensure no dividing by 0
gridNullEvents_95[gridNumRecs == 1] = NA
gridNumRecs[gridNumRecs == 1] = NA
gridFracEvents = gridEvents / gridNumRecs
gridFracEvents_NULL = gridEvents / gridNullEvents_95
#gridFracEvents_NULL[gridFracEvents_NULL > 1e10] = 0

# Calculate fraction of null events lower than real # of events
gridPercentEvents_null = matrix(NA, nrow = length(longitude), ncol = length(latitude))
for (i in 1:length(longitude)) {
  for (j in 1:length(latitude)) {
    
    gridPercentEvents_null[i,j] = sum(gridNullEvents[i,j,] < gridEvents[i,j]) / numIt
    
  }
}

# pull out non NA values
locs = which(!is.na(gridFracEvents), arr.ind = T)
fracEvents = gridFracEvents[locs]
numEvents = gridEvents[locs]
nullEvents = gridNullEvents_95[locs]
fracEvents_NULL = gridFracEvents_NULL[locs]
percentEvents_NULL = gridPercentEvents_null[locs]

## ------------------------------------ PLOTTING ------------------------------------ ##
# Note: work in progress
ggplot() + borders("world", colour="gray50", fill="gray50") + 
  geom_point(aes(x=longitude[locs[,1]], y=latitude[locs[,2]], color=fracEvents), size=3) +
  scale_color_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon, y = allTsLat), color = 'black', size = 1) +
  ggtitle(paste0('Fraction of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events')) +
  xlim(-180, 180) + ylim(-90, 90) +
  ylab('Latitude') + xlab('Longitude')
  
ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = fracEvents)) + 
  borders("world", colour="grey50") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Fraction of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))

ggplot() + borders("world", colour="gray50", fill="gray50") + 
  geom_raster(aes(x=longitude[locs[,1]], y=latitude[locs[,2]], fill=numEvents)) +
  scale_fill_distiller(name = 'Number\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon, y = allTsLat), color = 'black', size = 1) +
  ggtitle(paste0('Number of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events')) +
  xlim(-180, 180) + ylim(-90, 90) +
  ylab('Latitude') + xlab('Longitude')
  
# Do this for the null as well
ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = nullEvents)) + 
  borders("world", colour="black") + scale_fill_distiller(name = '# of null\n events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Null: 95% Q of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = fracEvents_NULL)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Real events / null events for', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = percentEvents_NULL)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)], color='no event'), size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)], color='event'), size = 1) +
  scale_color_manual(name = '', values = c('no event' = 'blue', 'event' = 'red')) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('BS: Fraction of null events < real event #\n', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = percentEvents_NULL)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'YlOrRd', direction = 1) + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)], color='no event'), size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)], color='event'), size = 1) +
  scale_color_manual(name = '', values = c('no event' = 'blue', 'event' = 'red')) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('BS: Fraction of null events < real event #\n', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))


binary_frac = fracEvents_NULL
binary_frac[binary_frac > 1] = 1
binary_frac[binary_frac <= 1] = 0
ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = factor(binary_frac))) + 
  borders("world", colour="black") + 
  scale_fill_manual(name = '95% Confidence\n          Level', values = c('0' = 'grey70', '1' = 'steelblue3'), labels = c('insignificant','significant')) + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)], color = 'no event'), size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)], color = 'event'), size = 1) +
  scale_color_manual(name = '', values = c('no event' = 'blue', 'event' = 'red')) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Broken Stick: Signficant ', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))


