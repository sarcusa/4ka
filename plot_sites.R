plot_sites <- function(data){
  
  ex = as.numeric(unlist(sapply(data,"[[","useEX")))
  data_EX = filterTs(data, 'useEX == 1')
  data_MS = filterTs(data, 'useMS == 1')
  data_BS = filterTs(data, 'useBS == 1')
  
  lat_EX = as.numeric(unlist(sapply(data_EX,"[[","geo_latitude")))
  lon_EX = as.numeric(unlist(sapply(data_EX,"[[","geo_longitude")))
  
  lat_MS = as.numeric(unlist(sapply(data_MS,"[[","geo_latitude")))
  lon_MS = as.numeric(unlist(sapply(data_MS,"[[","geo_longitude")))
  
  lat_BS = as.numeric(unlist(sapply(data_BS,"[[","geo_latitude")))
  lon_BS = as.numeric(unlist(sapply(data_BS,"[[","geo_longitude")))
  
  p <- ggplot() + borders("world", colour="black") + 
    geom_point(aes(x = lon_EX, y = lat_EX, color = 'EX'), size = 3, 
               shape = 0, stroke = 1.2) +
    geom_point(aes(x = lon_MS, y = lat_MS, color = 'MS/BS'), size = 3, 
               shape = 4, stroke = 1.2) +
    scale_color_manual(name = '', values = c('EX' = 'blue', 'MS/BS' = 'red'),
                       guide = guide_legend(override.aes = list(shape = c(0,4)))) +
    theme_bw() + xlab('Longitude') + ylab('Latitude') +
    xlim(-180, 180) + ylim(-90, 90) +
    ggtitle('Sites')
  
  return(p)
  
}