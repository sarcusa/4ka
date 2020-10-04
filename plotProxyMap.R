ProxyMap  <- function(prep1, prep2, param, input_var){
  
  datDir = file.path(createPaths(), 'RData')
  figDir = file.path(createPaths(), 'histograms')
  
  climateVar = input_var
  
  #Manual load
  load("/projects/pd_lab/sha59/4ka/RData/MS_results_plusNull_complete.RData")
  TS_MS = analysis_2b
  
  # Uncheck for automatic
  #TS_MS = prep1
  
  # avoid double-counting sites if they have annual and seasonal records
  for (i in 1:length(TS_MS)) {
    
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
  
  # isolate only the records corresponding to the chosen climate interpretation
  interps = unlist(sapply(TS_MS,"[[","interpretation1_variable"))
  if (climateVar == 'M') {
    inds = which(interps == 'M' | interps == 'P'| interps == 'P-E' | interps ==  'P/E')
  } else {
    inds = which(interps == 'T' | interps == 'TM')
  
  }
  TS_all  <- TS_MS
  for(i in 1:length(TS_all)){
    TS_all[[i]]$plot <- -1
  }
  for(i in 1:length(inds)){
  TS_all[[inds[i]]]$plot <- 1
  }
  
  lats = as.numeric(unlist(sapply(TS_MS,"[[","geo_latitude"))[inds])
  lons = as.numeric(unlist(sapply(TS_MS,"[[","geo_longitude"))[inds])
  names = unlist(sapply(TS_MS,"[[","dataSetName"))[inds]
  archives = unlist(sapply(TS_MS,"[[","archiveType"))[inds]
  archives[which(archives == 'marnine sediment' | archives == 'marine' | archives == 'TBD')] = 'marine sediment'
   
  missing_proxies = vector()
  for (i in 1:length(TS_MS)) {
    if(length(TS_MS[[i]]$paleoData_proxyGeneral) == 0){
      missing_proxies[i] = TS_MS[[i]]$dataSetName
    }
  }
  
  for (i in 1:length(TS_MS)) {
    print(TS_MS[[i]]$paleoData_proxyGeneral)
    if (length(TS_MS[[i]]$paleoData_proxyGeneral) == 0) {
      TS_MS[[i]]$paleoData_proxyGeneral = NA
    }
  }
  proxies = unlist(sapply(TS_MS,"[[","paleoData_proxyGeneral"))[inds]
  maxYrs = unlist(sapply(TS_MS, "[[", "maxYear"))[inds]
  
  missing_maxYr  <- vector()
  for (i in 1:length(TS_MS)) {
    if(length(TS_MS[[i]]$maxYear) == 0){
      missing_maxYr[i] = TS_MS[[i]]$dataSetName
    }
  }
  
  
  
  allDF = data.frame(name = names, lat = lats, lon = lons,
                     archive = archives, proxy = proxies, maxYr = maxYrs)
  
  eventYrs = param$eventYrs
  
  for (y in 1:length(eventYrs)) {
    
    eventYr = eventYrs[y]
    
    # must uncomment for automatic
    #TS = prep2
     
    #use to plot manually
    load(file.path(datDir, paste0('EX_results_plusNull_complete_', eventYr/1000, '.RData')))
    TS = data_EX
    ##
    
    
    # isolate only the records corresponding to the chosen climate interpretation
    interp = unlist(sapply(TS,"[[","interpretation1_variable"))
    if (climateVar == 'M') {
      inds_ex = which(interp == 'M' | interp == 'P'| interp == 'P-E' | interp ==  'P/E')
    } else if (climateVar == 'T')  {
      inds_ex = which(interp == 'T' | interp == 'TM')
    } else {
      inds_ex = which(interp == 'T' | interp == 'TM' | interp == 'M' | interp == 'P'| interp == 'P-E' | interp ==  'P/E')
    }
      
    test <- which(inds_ex %in% inds)
    if(!length(test) == length(inds_ex)){
      for(i in 1:length(inds_ex)){
        if(!inds_ex[i] %in% inds){
          TS_all[[inds_ex[i]]]$plot <- 1 
          print(paste0("adding ",TS_all[[inds_ex[i]]]$dataSetName))
        }
      }
    }
    
    lat = as.numeric(unlist(sapply(TS,"[[","geo_latitude"))[inds_ex])
    lon = as.numeric(unlist(sapply(TS,"[[","geo_longitude"))[inds_ex])
    name = unlist(sapply(TS,"[[","dataSetName"))[inds_ex]
    archive = unlist(sapply(TS,"[[","archiveType"))[inds_ex]
    
    more_proxies  <- vector()
    for (i in 1:length(TS)) {
      if(length(TS[[i]]$paleoData_proxyGeneral) == 0){
        more_proxies[i] = TS[[i]]$dataSetName
      }
    }
    
    for (i in 1:length(TS)) {
      print(TS[[i]]$paleoData_proxyGeneral)
      if (length(TS[[i]]$paleoData_proxyGeneral) == 0) {
        TS[[i]]$paleoData_proxyGeneral = NA
      }
    }
    proxy = unlist(sapply(TS,"[[","paleoData_proxyGeneral"))[inds_ex]
    
    maxYr = unlist(sapply(TS_MS, "[[", "maxYear"))[inds_ex]
    
    more_maxYrs  <- vector()
    for (i in 1:length(TS_MS)) {
      if(length(TS_MS[[i]]$maxYear) == 0){
        missing_maxYr[i] = TS_MS[[i]]$dataSetName
      }
    }
    
    dfAdd = data.frame(name, lat, lon, archive, proxy, maxYr)
    allDF = rbind(allDF, dfAdd)
        
    additions = which(!name %in% names)
    if (length(additions) > 0) {
      lats = c(lats, lat[additions])
      lons = c(lons, lon[additions])
      names = c(names, name[additions])
      archives = c(archives, archive[additions])
      proxies = c(proxies, proxy[additions]) # added 09-18-20 
      missing_proxies = c(missing_proxies, more_proxies[additions])#added 09-18-20
      maxYrs = c(maxYrs, more_maxYrs[additions]) #added 09-18-20
    }
    
    
  }
  
  unique(missing_proxies)
  write.csv(x = unique(missing_proxies), paste0(datDir,"/missing_proxyGeneral_", climateVar,".csv"))  
  
  unique(missing_maxYr)
  write.csv(x = unique(missing_maxYr), paste0(datDir,"/missing_maxYr_", climateVar,".csv")) 
  
  archives[which(archives == 'MarineSediment')] = 'marine sediment'
  archives[which(archives == 'LakeSediment')] = 'lake sediment'
  archives[which(archives == 'Peat')] = 'peat'
  archives[which(archives == 'Wood')] = 'wood'
  archives[which(archives == 'Midden')] = 'midden'
  archives[which(archives == 'Speleothem')] = 'speleothem'
  archives[which(archives == 'GlacierIce')] = 'glacier ice'
  archives[which(archives == 'Ice-other')] = 'ice-other'
  
  lons = lons[which(!is.na(archives))]
  lats = lats[which(!is.na(archives))]
  archives = archives[which(!is.na(archives))]
  maxYear = maxYrs[which(!is.na(archives))]
  unique(archives)
  
  df = data.frame(lats, lons, archives, maxYear)
  col_m = c('#1c9099', '#08519c', '#810f7c', '#a63603', '#006d2c', '#c6dbef', '#f16913', '#74c476', '#8c96c6')
  sh_m = c(seq(0,length(unique(archives))-1,1))
  
  a  <- ggplot() + borders("world", colour="black") + 
    geom_point(aes(x = lons, y = lats, fill = archives, shape = archives)) +
    scale_shape_manual(values = sh_m)
  
  pdf(file.path(figDir, 'All_archives.pdf'))
  print(a)
  dev.off()
  
    
  if (climateVar == 'M') {
  
    
    
    
  ## Moisture
  sss = c(seq(-90,90), seq(90,-90, by = -1))
  boundcirc = data.frame(y = sss, 
                         x = c(rep(-180, length.out = length(sss)/2), 
                               rep(180, length.out=length(sss)/2)))
  tryCatch({t  <- baseMAP(c(-180,180), c(-90,90),map.type='line',
          global=T,projection='mollweide',restrict.map.range=F, 
          country.boundaries = F) + 
    geom_path(aes(x = boundcirc$x, y = sss)) +
      geom_point(aes(x = lons[which(archives == 'midden')], 
                     y = lats[which(archives == 'midden')], 
                     fill = 'midden'), shape = 23, size = 2, color = "white") +
      geom_point(aes(x = lons[which(archives == 'marine sediment')], 
                     y = lats[which(archives == 'marine sediment')], 
                     fill = 'marine sediment'),  shape = 21, size = 2, 
                 color = "white") +
      geom_point(aes(x = lons[which(archives == 'lake sediment')],
                     y = lats[which(archives == 'lake sediment')], 
                     fill = 'lake sediment'), shape = 22, size = 2, 
                 color = "white") +
      geom_point(aes(x = lons[which(archives == 'speleothem')], 
                     y = lats[which(archives == 'speleothem')], 
                     fill = 'speleothem'), shape = 24, size = 2, 
                 color = "white") +
      geom_point(aes(x = lons[which(archives == 'glacier ice')], 
                     y = lats[which(archives == 'glacier ice')], 
                     fill = 'glacier ice'), shape = 21, size = 2,
                 color = "white") +
      geom_point(aes(x = lons[which(archives == 'peat')], 
                     y = lats[which(archives == 'peat')], 
                     fill = 'peat'), shape =22, size = 2, color = "white") +
      geom_point(aes(x = lons[which(archives == 'wood')], 
                     y = lats[which(archives == 'wood')], 
                     fill = 'wood'), shape = 24, size = 2, color = "white") +
      scale_fill_manual(name = '', values = c('midden' = col_m[3],
                                              'marine sediment' = col_m[1],
                                              'lake sediment' = col_m[2],
                                              'wood' = col_m[8],
                                              'speleothem' = col_m[4],
                                              'glacier ice' = col_m[6],
                                              'peat' = col_m[7]), 
    guide = guide_legend(override.aes = list(shape = c(21,22,21,23,22,24,24), fill = col_m[c(6,2,1,3,7,4,8)],size = rep(3,7), color = rep('white',7)))) +
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(), 
            axis.text.y = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            panel.border = element_blank()) +
      ggtitle('Moisture proxies')
      
    #theme_bw() + xlab('Longitude') + ylab('Latitude') +
    #theme(plot.title = element_text(hjust = 0.5)) +
    #xlim(-180, 180) + ylim(-90, 90) +
    #ggtitle('Moisture proxies')
  
  pdf(file.path(figDir, 'moistureProxyMap_v1.pdf'))
  print(t)
  dev.off()
  
  #file.copy(from = file.path(figDir, 'moistureProxyMap_v1.pdf'),overwrite = T, 
  #           to = file.path(param$mainDir,"report_files/figure-html/"))
  
  }, error = function(e){cat("Error:", conditionMessage(e), " no plot was produced ")})
  
  Time.plot <- plotTimeAvailabilityTs(TS = filterTs(TS_all, expression = 'plot == 1'), group.var = "archiveType")+
    scale_fill_manual(name = '', values = c('Midden' = col_m[3],
                                            'MarineSediment' = col_m[1],
                                            'LakeSediment' = col_m[2],
                                            'Wood' = col_m[8],
                                            'Speleothem' = col_m[4],
                                            'GlacierIce' = col_m[6],
                                            'Peat' = col_m[7]), 
                      guide = guide_legend(override.aes = list(fill = col_m[c(6,2,1,3,7,4,8)])))+
    theme_bw()
  
  pdf(file.path(figDir, 'MoistureArchiveTime.pdf'))
  print(Time.plot)
  dev.off()
  
  
  }
  
  if (climateVar == 'T') {
  
    
    ## Temperature
  sss = c(seq(-90,90), seq(90,-90, by = -1))
  boundcirc = data.frame(y = sss, x = c(rep(-180, 
                                            length.out = length(sss)/2), 
                                        rep(180, length.out=length(sss)/2)))
  
  tryCatch({t  <- baseMAP(c(-180,180), c(-90,90),map.type='line',global=T,
          projection='mollweide',restrict.map.range=F, country.boundaries = F) + 
    geom_path(aes(x = boundcirc$x, y = sss)) +
    geom_point(aes(x = lons[which(archives == 'midden')], 
                   y = lats[which(archives == 'midden')], 
                   fill = 'midden'), shape = 23, size = 2, color = "white") +
    geom_point(aes(x = lons[which(archives == 'marine sediment')], 
                   y = lats[which(archives == 'marine sediment')], 
                   fill = 'marine sediment'),  shape = 21, size = 2, 
               color = "white") +
    geom_point(aes(x = lons[which(archives == 'lake sediment')],
                   y = lats[which(archives == 'lake sediment')], 
                   fill = 'lake sediment'), shape = 22, size = 2, 
               color = "white") +
    geom_point(aes(x = lons[which(archives == 'speleothem')], 
                   y = lats[which(archives == 'speleothem')], 
                   fill = 'speleothem'), shape = 24, size = 2, color = "white") +
    geom_point(aes(x = lons[which(archives == 'glacier ice')], 
                   y = lats[which(archives == 'glacier ice')], 
                   fill = 'glacier ice'), shape = 21, size = 2, color = "white") +
    geom_point(aes(x = lons[which(archives == 'peat')], 
                   y = lats[which(archives == 'peat')], 
                   fill = 'peat'), shape =22, size = 2, color = "white") +
    geom_point(aes(x = lons[which(archives == 'wood')], 
                   y = lats[which(archives == 'wood')], 
                   fill = 'wood'), shape = 24, size = 2, color = "white") +
    geom_point(aes(x = lons[which(archives == 'ice-other')], 
                   y = lats[which(archives == 'ice-other')], 
                   fill = 'ice-other'), shape = 23, size = 2, color = "white") +
    scale_fill_manual(name = '', values = c('midden' = col_m[3],
                                            'marine sediment' = col_m[1],
                                            'lake sediment' = col_m[2],
                                            'wood' = col_m[8],
                                            'speleothem' = col_m[4],
                                            'glacier ice' = col_m[6],
                                            'ice-other'=col_m[9],
                                            'peat' = col_m[7]), 
                      guide = guide_legend(override.aes = list(shape = c(21,23,22,21,23,22,24,24), fill = col_m[c(6,9,2,1,3,7,4,8)], size = rep(3,8), color = rep('white',8)))) +
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(), 
            axis.text.y = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            panel.border = element_blank()) +
      ggtitle('Temperature proxies')
    
    #theme_bw() + 
    #xlab('Longitude') + ylab('Latitude') +
    #theme(plot.title = element_text(hjust = 0.5)) +
    #xlim(-180, 180) + ylim(-90, 90) +
    #ggtitle('Temperature proxies')
  
  pdf(file.path(figDir, 'temperatureProxyMap_v2.pdf'))
  print(t)
  dev.off()
  
  #file.copy(from = file.path(figDir, 'temperatureProxyMap_v2.pdf'),overwrite = T, 
  #           to = file.path(param$mainDir,"report_files/figure-html/"))
  
  }, error = function(e){cat("Error:", conditionMessage(e), " no plot was produced ")})
  
  Time.plot <- plotTimeAvailabilityTs(TS = filterTs(TS_all, expression = 'plot == 1'), group.var = "archiveType")+
    scale_fill_manual(name = '', values = c('Midden' = col_m[3],
                                            'MarineSediment' = col_m[1],
                                            'LakeSediment' = col_m[2],
                                            'Wood' = col_m[8],
                                            'Speleothem' = col_m[4],
                                            'GlacierIce' = col_m[6],
                                            'ice-other'=col_m[9],
                                            'Peat' = col_m[7]), 
   guide = guide_legend(override.aes = list(fill = col_m[c(6,9,2,1,3,7,4,8)])))+
    theme_bw()
  
  pdf(file.path(figDir, 'temperatureArchiveTime.pdf'))
  print(Time.plot)
  dev.off()
  
  }
  
  # Sum these categories in the output file for temperature and moisture
  print(paste('wood:', length(archives[which(archives == 'wood')])))
  print(paste('marine sediment:', 
              length(archives[which(archives == 'marine sediment')])))
  print(paste('lake sediment:', 
              length(archives[which(archives == 'lake sediment')])))
  print(paste('speleothem:', length(archives[which(archives == 'speleothem')])))
  print(paste('ice-other:', length(archives[which(archives == 'ice-other')])))
  print(paste('glacier ice:', length(archives[which(archives == 'glacier ice')])))
  print(paste('midden:', length(archives[which(archives == 'midden')])))
  print(paste('peat:', length(archives[which(archives == 'peat')])))
  print(paste('Total records:', length(archives)))
  print(paste('Total sites:', length(unique(names))))
  
  archs = c('wood','marine sediment','lake sediment','speleothem','ice-other','glacier ice','midden','peat','records','sites')
  num_archs = c(length(archives[which(archives == 'wood')]),
                length(archives[which(archives == 'marine sediment')]),
                length(archives[which(archives == 'lake sediment')]),
                length(archives[which(archives == 'speleothem')]),
                length(archives[which(archives == 'ice-other')]),
                length(archives[which(archives == 'glacier ice')]),
                length(archives[which(archives == 'midden')]),
                length(archives[which(archives == 'peat')]),
                length(archives),length(unique(names)))
  rec_df = data.frame(archs, num_archs)
  
  # ONLY UNIQUE SITES
  names_unique_inds = order(names)[!duplicated(sort(names))]
  names_unique = names[names_unique_inds]
  archives2 = archives[names_unique_inds]
  num_archs2 = c(length(archives2[which(archives2 == 'wood')]),
                 length(archives2[which(archives2 == 'marine sediment')]),
                 length(archives2[which(archives2 == 'lake sediment')]),
                 length(archives2[which(archives2 == 'speleothem')]),
                 length(archives2[which(archives2 == 'ice-other')]),
                 length(archives2[which(archives2 == 'glacier ice')]),
                 length(archives2[which(archives2 == 'midden')]),
                 length(archives2[which(archives2 == 'peat')]),
                 length(archives2),length(unique(names)))
  rec_df2 = data.frame(archs, num_archs2)
  
  if (climateVar == 'T') {
    write.csv(rec_df2, file.path(figDir, 'temperature_records_unique.csv'))
  } else if (climateVar == "M") {
    write.csv(rec_df2, file.path(figDir, 'moisture_records_unique.csv'))
  } else {
    write.csv(rec_df, file.path(figDir, "all_records_unique.csv"))
  }
  
  # Proxies
  
  lons = lons[which(!is.na(proxies))]
  lats = lats[which(!is.na(proxies))]
  proxies = proxies[which(!is.na(proxies))]
  unique(proxies)
  
  df = data.frame(lats, lons, proxies)
  col_m = c('#1c9099', '#08519c', '#810f7c', '#a63603', '#006d2c', '#c6dbef', '#f16913', '#74c476', '#8c96c6')
  sh_m = c(seq(0,length(unique(proxies))-1,1))
  
  z  <- ggplot() + borders("world", colour="black") + 
    geom_point(aes(x = lons, y = lats, fill = proxies, shape = proxies)) +
    scale_shape_manual(values = sh_m)
  
  pdf(file.path(figDir, 'All_proxies.pdf'))
  print(z)
  dev.off()
  
  if (climateVar == 'M') {
    
    ## Moisture
    sss = c(seq(-90,90), seq(90,-90, by = -1))
    boundcirc = data.frame(y = sss, 
                           x = c(rep(-180, length.out = length(sss)/2), 
                                 rep(180, length.out=length(sss)/2)))
    tryCatch({tt  <- baseMAP(c(-180,180), c(-90,90),map.type='line',
                            global=T,projection='mollweide',restrict.map.range=F, 
                            country.boundaries = F) + 
                geom_path(aes(x = boundcirc$x, y = sss)) +
                geom_point(aes(x = lons[which(proxies == 'other microfossil')], 
                               y = lats[which(proxies == 'other microfossil')], 
                               fill = 'other microfossil'), shape = 23, size = 2, color = "white") +
                geom_point(aes(x = lons[which(proxies == 'alkenone')], 
                               y = lats[which(proxies == 'alkenone')], 
                               fill = 'alkenone'),  shape = 21, size = 2, 
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'Mg/Ca')],
                               y = lats[which(proxies == 'Mg/Ca')], 
                               fill = 'Mg/Ca'), shape = 22, size = 2, 
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'pollen')], 
                               y = lats[which(proxies == 'pollen')], 
                               fill = 'pollen'), shape = 24, size = 2, 
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'chironomid')], 
                               y = lats[which(proxies == 'chironomid')], 
                               fill = 'chironomid'), shape = 21, size = 2,
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'isotope')], 
                               y = lats[which(proxies == 'isotope')], 
                               fill = 'isotope'), shape =22, size = 2, color = "white") +
                geom_point(aes(x = lons[which(proxies == 'other ice')], 
                               y = lats[which(proxies == 'other ice')], 
                               fill = 'other ice'), shape = 24, size = 2, color = "white") +
                geom_point(aes(x = lons[which(proxies == 'physical')], 
                               y = lats[which(proxies == 'physical')], 
                               fill = 'physical'), shape = 23, size = 2, color = "white")+
                geom_point(aes(x = lons[which(proxies == 'biophysical')], 
                               y = lats[which(proxies == 'biophysical')], 
                               fill = 'biophysical'), shape = 23, size = 2, color = "white")+
                
                geom_point(aes(x = lons[which(proxies == 'diatom')], 
                               y = lats[which(proxies == 'diatom')], 
                               fill = 'diatom'), shape = 25, size = 2, color = "white")+
                
                scale_fill_manual(name = '', values = c('pollen' = col_m[8],
                                                        'isotope' = col_m[6],
                                              'other microfossil' = col_m[3],
                                                        'physical' = col_m[8],
                                                        'biophysical' = col_m[5],
                                                        'diatom' = col_m[2],
                                                        'Mg/Ca' = col_m[2], 
                                  'chironomid' = col_m[4],
                                  'other ice' = col_m[9],
                                  'alkenone' = col_m[1]),
                                  guide = guide_legend(override.aes = list(shape = c(21,23,21,25,22,22,24,23,23,24), fill = col_m[c(1,5,4,2,6,2,9,3,8,8)],size = rep(3,10), color = rep('white',10)))) +
                theme_bw()+
                theme(plot.title = element_text(hjust = 0.5),
                      axis.ticks = element_blank(), 
                      axis.text.y = element_blank(), 
                      axis.text.x = element_blank(), 
                      axis.title.x = element_blank(), 
                      axis.title.y = element_blank(), 
                      panel.border = element_blank()) +
                ggtitle('Moisture proxies')
              
             pdf(file.path(figDir, 'moistureProxyGeneral.pdf'))
              print(tt)
              dev.off()
              
                          
    }, error = function(e){cat("Error:", conditionMessage(e), " no plot was produced ")})
    
    Time.plot2 <- plotTimeAvailabilityTs(TS = filterTs(TS_all, expression = 'plot == 1'), group.var = "paleoData_proxyGeneral")+
      scale_fill_manual(name = '', values = c('pollen' = col_m[8],
                                              'isotope' = col_m[6],
                                              'other microfossil' = col_m[3],
                                              'physical' = col_m[8],
                                              'biophysical' = col_m[5],
                                              'diatom' = col_m[2],
                                              'Mg/Ca' = col_m[2], 
                                              'chironomid' = col_m[4],
                                              'other ice' = col_m[9],
                                              'alkenone' = col_m[1],
                                              'NULL' = ""),
                        guide = guide_legend(override.aes = list(fill = col_m[c(1,5,4,2,6,2,9,3,8,8)])))+
      theme_bw()
    
    pdf(file.path(figDir, 'MoistureProxyTime.pdf'))
    print(Time.plot2)
    dev.off()
    
  }
  
  if (climateVar == 'T') {
    ## Temperature
    sss = c(seq(-90,90), seq(90,-90, by = -1))
    boundcirc = data.frame(y = sss, x = c(rep(-180, 
                                              length.out = length(sss)/2), 
                                          rep(180, length.out=length(sss)/2)))
    
    tryCatch({tt  <- baseMAP(c(-180,180), c(-90,90),map.type='line',global=T,
                            projection='mollweide',restrict.map.range=F, country.boundaries = F) + 
                geom_point(aes(x = lons[which(proxies == 'other microfossil')], 
                               y = lats[which(proxies == 'other microfossil')], 
                               fill = 'other microfossil'), shape = 23, size = 2, color = "white") +
                geom_point(aes(x = lons[which(proxies == 'alkenone')], 
                               y = lats[which(proxies == 'alkenone')], 
                               fill = 'alkenone'),  shape = 21, size = 2, 
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'Mg/Ca')],
                               y = lats[which(proxies == 'Mg/Ca')], 
                               fill = 'Mg/Ca'), shape = 22, size = 2, 
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'pollen')], 
                               y = lats[which(proxies == 'pollen')], 
                               fill = 'pollen'), shape = 24, size = 2, 
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'chironomid')], 
                               y = lats[which(proxies == 'chironomid')], 
                               fill = 'chironomid'), shape = 21, size = 2,
                           color = "white") +
                geom_point(aes(x = lons[which(proxies == 'isotope')], 
                               y = lats[which(proxies == 'isotope')], 
                               fill = 'isotope'), shape =22, size = 2, color = "white") +
                geom_point(aes(x = lons[which(proxies == 'other ice')], 
                               y = lats[which(proxies == 'other ice')], 
                               fill = 'other ice'), shape = 24, size = 2, color = "white") +
                geom_point(aes(x = lons[which(proxies == 'other biomarker')], 
                               y = lats[which(proxies == 'other biomarker')], 
                               fill = 'other biomarker'), shape = 24, size = 2, color = "white")+
                geom_point(aes(x = lons[which(proxies == 'biophysical')], 
                               y = lats[which(proxies == 'biophysical')], 
                               fill = 'biophysical'), shape = 23, size = 2, color = "white")+
                geom_point(aes(x = lons[which(proxies == 'GDGT')], 
                               y = lats[which(proxies == 'GDGT')], 
                               fill = 'GDGT'), shape = 25, size = 2, color = "white")+
                scale_fill_manual(name = '', values = c('other microfossil' = col_m[3],'alkenone' = col_m[1],'Mg/Ca' = col_m[2],'pollen' = col_m[8],'chironomid' = col_m[4],'isotope' = col_m[6],'other ice'=col_m[9],'other biomarker' = col_m[7], 'biophysical' = col_m[5], 'GDGT' = col_m[3]), 
                                  guide = guide_legend(override.aes = list(shape = c(21,23,21,25,22,22,24,24,23,24), fill = col_m[c(1,5,4,3,6,2,7,9,3,8)], size = rep(3,10), color = rep('white',10)))) +
                theme_bw()+
                theme(plot.title = element_text(hjust = 0.5),
                      axis.ticks = element_blank(), 
                      axis.text.y = element_blank(), 
                      axis.text.x = element_blank(), 
                      axis.title.x = element_blank(), 
                      axis.title.y = element_blank(), 
                      panel.border = element_blank()) +
                ggtitle('Temperature proxies')
              
              pdf(file.path(figDir, 'temperatureProxyGeneral.pdf'))
              print(tt)
              dev.off()
              
                           
    }, error = function(e){cat("Error:", conditionMessage(e), " no plot was produced ")})
    
    Time.plot2 <- plotTimeAvailabilityTs(TS = filterTs(TS_all, expression = 'plot == 1'), group.var = "paleoData_proxyGeneral")+
      scale_fill_manual(name = '', values = c('other microfossil' = col_m[3],
                                              'alkenone' = col_m[1],
                                              'Mg/Ca' = col_m[2],
                                              'pollen' = col_m[8],
                                              'chironomid' = col_m[4],
                                              'isotope' = col_m[6],
                                              'other ice'=col_m[9],
                                              'other biomarker' = col_m[7],
                                              'biophysical' = col_m[5],
                                              'GDGT' = col_m[3]), 
                        guide = guide_legend(override.aes = list(fill = col_m[c(1,5,4,3,6,2,7,9,3,8)]))) +
      theme_bw()
    
    pdf(file.path(figDir, 'TemperatureProxyTime.pdf'))
    print(Time.plot2)
    dev.off()
    
  }
  
  # Sum these categories for temperature and moisture
  print(paste('other microfossil:', length(proxies[which(proxies == 'other microfossil')])))
  print(paste('alkenone:', 
              length(proxies[which(proxies == 'alkenone')])))
  print(paste('Mg/Ca:', 
              length(proxies[which(proxies == 'Mg/Ca')])))
  print(paste('pollen:', length(proxies[which(proxies == 'pollen')])))
  print(paste('chironomid:', length(proxies[which(proxies == 'chironomid')])))
  print(paste('isotope:', length(proxies[which(proxies == 'isotope')])))
  print(paste('other ice:', length(proxies[which(proxies == 'other ice')])))
  print(paste('other biomarker:', length(proxies[which(proxies == 'other biomarker')])))
  print(paste('biophysical:', length(proxies[which(proxies == 'biophysical')])))
  print(paste('GDGT:', length(proxies[which(proxies == 'GDGT')])))
  print(paste('physical:', length(proxies[which(proxies == 'physical')])))
  print(paste('diatom:', length(proxies[which(proxies == 'diatom')])))
  print(paste('Total proxies:', length(proxies)))
  print(paste('Total sites:', length(unique(names))))
  
  prox = c('other microfossil','alkenone', 'Mg/Ca','pollen','chironomid','isotope', 'other ice ','other biomarker','biophysical','GDGT', 'physical' , 'diatom','records_proxies','sites')
  num_prox = c(length(proxies[which(proxies == 'other microfossil')]),
               length(proxies[which(proxies == 'alkenone')]),
               length(proxies[which(proxies == 'Mg/Ca')]),
               length(proxies[which(proxies == 'pollen')]),
               length(proxies[which(proxies == 'chironomid')]),
               length(proxies[which(proxies == 'isotope')]),
               length(proxies[which(proxies == 'other ice')]),
               length(proxies[which(proxies == 'other biomarker')]),
               length(proxies[which(proxies == 'biophysical')]),
               length(proxies[which(proxies == 'GDGT')]),
               length(proxies[which(proxies == 'physical')]),
               length(proxies[which(proxies == 'diatom')]),
                length(proxies),length(unique(names)))
  rec_df = data.frame(prox, num_prox)
  
  # ONLY UNIQUE SITES
  names_unique_inds = order(names)[!duplicated(sort(names))]
  names_unique = names[names_unique_inds]
  proxies2 = proxies[names_unique_inds]
  num_prox2 = c(length(proxies2[which(proxies2 == 'other microfossil')]),
                length(proxies2[which(proxies2 == 'alkenone')]),
                 length(proxies2[which(proxies2 == 'Mg/Ca')]),
                 length(proxies2[which(proxies2 == 'pollen')]),
                 length(proxies2[which(proxies2 == 'chironomid')]),
                 length(proxies2[which(proxies2 == 'isotope')]),
                length(proxies2[which(proxies2 == 'other ice')]),
                 length(proxies2[which(proxies2 == 'other biomarker')]),
                 length(proxies2[which(proxies2 == 'biophysical')]),
                 length(proxies2[which(proxies2 == 'GDGT')]),
                length(proxies2[which(proxies2 == 'physical')]),
                length(proxies2[which(proxies2 == 'diatom')]),
                 length(proxies2),length(unique(names)))
  rec_df2 = data.frame(prox, num_prox2)
  
  if (climateVar == 'T') {
    write.csv(rec_df2, file.path(figDir, 'temperature_records_proxies_unique.csv'))
  } else if (climateVar == "M") {
    write.csv(rec_df2, file.path(figDir, 'moisture_records_proxies_unique.csv'))
  } else {
    write.csv(rec_df, file.path(figDir, "all_records_proxies_unique.csv"))
  }
  
  
  print("THE END")
  plots  <- list(a, t, tt)
  return(plots)
  
}

