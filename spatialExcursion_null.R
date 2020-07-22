## Written by Hannah Kolus, 09/04/2018 
## Map the results of excursion analysis on a global grid. Test the spatial significance of events 
## with a null model (analysis performed on an ensemble of synthetic data records with the same 
## autocorrelation, trend, and temporal resolution as the orginal records)

library(lipdR)
library(ggplot2)
library(ggmap)
library(geosphere)
#library(viridis)
#library(viridisLite)

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')

# parameters
res = 5          # grid resolution, deg
radius = 2000      # search radius, km
eventYr = 8200     # event year, yr BP
numIt = 1000

# Load in relevant results: TS_EX
load(paste0('RData/EX_results_plusNull2_', eventYr/1000, '.RData'))

# Make a grid of lat/lon center points
latitude = seq(-90 + res/2, 90 - res/2, res)
longitude = seq(-180 + res/2, 180 - res/2, res)

# Grids that will contain fraction of passes
gridEvents = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNumRecs = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNullEvents = array(NA, dim = c(length(longitude), length(latitude), numIt))
totNullEvents = matrix(NA, nrow = length(TS_EX), ncol = numIt)

# Assign event occurrence
for (i in 1:length(TS_EX)) {
  
  totNullEvents[i,] = sapply(TS_EX[[i]]$null_events, sum)
  
}

allTsLat = as.numeric(sapply(TS_EX,"[[","geo_latitude"))
allTsLon = as.numeric(sapply(TS_EX,"[[","geo_longitude"))
events = as.numeric(sapply(TS_EX,"[[","eventEX"))

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
gridEvents[gridEvents == 0] = .Machine$double.xmin
gridNullEvents_95 = apply(gridNullEvents, c(1,2), function (x) quantile(x, 0.95))
#gridNullEvents_95[gridNullEvents_95 == 0] = .Machine$double.xmin # ensure no dividing by 0
gridNullEvents_95[gridNumRecs == 1] = NA
gridNumRecs[gridNumRecs == 1] = NA
gridFracEvents = gridEvents / gridNumRecs
gridFracEvents_NULL = gridEvents / gridNullEvents_95

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
df = data.frame(lon = longitude[locs[,1]], lat = latitude[locs[,2]],frac = fracEvents)
ggplot(df) + borders("world", colour="gray50", fill="gray50") + 
  geom_point(aes(x=lon, y=lat, color=frac), size=3) +
  scale_color_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  ggtitle(paste0('Fraction of ', eventYr/1000,'+/-',event_window/1000, 'ka events')) +
  xlim(-180, 180) + ylim(-90, 90) +
  ylab('Latitude') + xlab('Longitude')

## using geoChronR, ggplot2, and ggmap
## baseMap(c(-180,180),c(-90,90),map.type = 'line',global = TRUE,projection='mollweide',f=0)

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = fracEvents)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Fraction of ', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = numEvents)) + 
  borders("world", colour="black") + 
  scale_fill_distiller(name = 'Number\nof events', palette = 'Spectral', limits = c(0, max(nullEvents))) + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Number of ', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))
  
# Do this for the null as well
ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = nullEvents)) + 
  borders("world", colour="black") + 
  scale_fill_distiller(name = '# of null\n events', palette = 'Spectral', limits = c(0, max(nullEvents))) + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Null: 95% Q of ', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = fracEvents_NULL)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Real events / null events for', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = percentEvents_NULL)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)], color='no event'), size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)], color='event'), size = 1) +
  scale_color_manual(name = '', values = c('no event' = 'blue', 'event' = 'red')) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('EX: Fraction of null events < real event #\n', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = percentEvents_NULL)) + 
  borders("world", colour="black") + 
  scale_fill_distiller(name = 'Fraction\nof events', palette = 'YlOrRd', direction = 1) + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color='white', size = 2) +
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)], color='no event'), size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)], color='event'), size = 1) +
  scale_color_manual(name = '', values = c('no event' = 'blue', 'event' = 'red')) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('EX: Fraction of null events < real event #\n', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))


binary_frac = fracEvents_NULL
binary_frac[binary_frac <= 1] = 0
binary_frac[binary_frac > 1] = 1
ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = factor(binary_frac))) + 
  borders("world", colour="black") + 
  scale_fill_manual(name = '95% Confidence\n          Level', values = c('0' = 'grey70', '1' = 'steelblue3'), labels = c('insignificant','significant')) + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)], color = 'no event'), size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)], color = 'event'), size = 1) +
  scale_color_manual(name = '', values = c('no event' = 'blue', 'event' = 'red')) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Excursion: Signficant ', eventYr/1000,'+/-',event_window/2/1000, 'ka events'))

