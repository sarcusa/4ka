## Written by Hannah Kolus, 09/04/2018 
## Map the results of excursion analysis on a global grid. Plot the number of events occurring within
## a certain radius of each grid cell. Plot the fraction of events relative to available records
## within a certain radius of each grid cell.

library(lipdR)
library(ggplot2)
library(ggmap)
library(geosphere)

setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')

# parameters
res = 5          # grid resolution, deg
radius = 1500      # search radius, km
eventYr = 8200     # event year, yr BP
eventWindow = 300  # total event window = eventYr +/- eventWindow

# Load in relevant results: TS_EX
load(paste0('RData/EX_results', eventYr/1000, '.RData'))
TS_EX = filterTs(TS_EX, 'statusEX == 1') # only use records that passed all criteria

# Make a grid of lat/lon center points
latitude = seq(-90 + res/2, 90 - res/2, res)
longitude = seq(-180 + res/2, 180 - res/2, res)

# Grid the will contain fraction of passes
gridEvents = matrix(NA, nrow = length(longitude), ncol = length(latitude))
gridNumRecs = matrix(NA, nrow = length(longitude), ncol = length(latitude))

allTsLat = as.numeric(sapply(TS_EX,"[[","geo_latitude"))
allTsLon = as.numeric(sapply(TS_EX,"[[","geo_longitude"))
events = as.numeric(sapply(TS_EX,"[[","eventEX"))

for (i in 1:length(longitude)) {
  
  print(paste('iteration i = ', i))
  
  for (j in 1:length(latitude)) {
    
    inds = which(distGeo(cbind(allTsLon, allTsLat), c(longitude[i],latitude[j])) <= radius*1000)
    
    gridEvents[i, j] = sum(gridEvents[i, j], events[inds], na.rm = T)
    gridNumRecs[i, j] = sum(gridNumRecs[i, j], length(inds), na.rm = T)
    
  }
  
}

# Calculate the fraction of sites within the radius that pass the test
gridFracEvents = gridEvents / gridNumRecs
locs = which(!is.na(gridFracEvents), arr.ind = T)
fracEvents = gridFracEvents[locs]
numEvents = gridEvents[locs]

# Plotting
# Note: work in progress
df = data.frame(lon = longitude[locs[,1]], lat = latitude[locs[,2]],frac = fracEvents)
ggplot(df) + borders("world", colour="gray50", fill="gray50") + 
  geom_point(aes(x=lon, y=lat, color=frac), size=3) +
  scale_color_distiller(name = 'Fraction\nof events', palette = 'Spectral') + 
  ggtitle(paste0('Fraction of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events')) +
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
  ggtitle(paste0('Fraction of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))

ggplot() + geom_raster(aes(x = longitude[locs[,1]], y = latitude[locs[,2]], fill = numEvents)) + 
  borders("world", colour="black") + scale_fill_distiller(name = 'Number\nof events', palette = 'Spectral') + 
  geom_point(aes(x = allTsLon[which(events == 0)], y = allTsLat[which(events == 0)]), color = 'blue', size = 1) +
  geom_point(aes(x = allTsLon[which(events == 1)], y = allTsLat[which(events == 1)]), color = 'red', size = 1) +
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  xlim(-180, 180) + ylim(-90, 90) +
  ggtitle(paste0('Number of ', eventYr/1000,'+/-',eventWindow/1000, 'ka events'))
  
# Do this for the null as well

