hist_MS_plot <- function(data_in,param, climateVar){
  
  eventDetector = 'MS'
  
  figDir = file.path(createPaths(), 'histograms')
  dataDir = file.path(createPaths(), 'RData')
  
  eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))
  eventDetectorStr = ifelse(eventDetector == 'MS', 'Mean shift', 'Broken stick')
  
  posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
  negCol = ifelse(climateVar == 'M', '#543005', '#053061')
  posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
  negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
  quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')
  
  eventYrs = param$eventYrs
  allEvents = data_in$allEvents
  allNullQuants = data_in$allNullQuants
  negEvents = data_in$negEvents
  posEvents = data_in$posEvents
  diffEvents = data_in$diffEvents
  recordCounts = data_in$recordCounts
  posNullQuants = data_in$posNullQuants
  negNullQuants = data_in$negNullQuants
  diffNullQuants = data_in$diffNullQuants
  allNullEvents = data_in$allNullEvents
  posNullEvents = data_in$posNullEvents
  negNullEvents = data_in$negNullEvents
  
  
  if (plotOpt) {
    
    if (climateVar == 'All') {
      
      s  <- ggplot() + geom_col(aes(x = eventYrs, y = allEvents),
                          fill = 'grey60', color = 'grey10') +
        geom_line(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        geom_point(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        scale_x_continuous(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom", 
                           legend.text = element_text(size = 12)) +
        ggtitle(paste0(eventDetectorStr, ':\nAll events'))
      pdf(file.path(figDir, paste0(eventDetector, '_', 
                                   climateVar,'.pdf')))
          print(s)
          dev.off()
          
          plots = s
      
    } else {
      
      p1 = ggplot() + geom_col(aes(x = eventYrs, y = allEvents), 
                               fill = 'grey60', color = 'grey10') +
        geom_line(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        geom_point(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        scale_x_continuous(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0(eventDetectorStr, ':\nAll ', eventTypeStr,' events'))
      
      p2 = ggplot() + geom_col(aes(x = eventYrs, y = posEvents), fill = posCol) +
        geom_line(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        geom_point(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        scale_x_continuous(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0('Positive ', eventTypeStr,' events'))
            
      p3 = ggplot() + geom_col(aes(x = eventYrs, y = negEvents), fill = negCol) +
        geom_line(aes(x = eventYrs, y = negNullQuants[1,], color = '0.9')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[2,], color = '0.95')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[3,], color = '0.99')) +
        geom_point(aes(x = eventYrs, y = negNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = negNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = negNullQuants[3,]), color = quantCol[3]) +
        scale_color_manual(name = 'Quantile', 
                           values = c('0.9' = quantCol[1], 
                                      '0.95' = quantCol[2], '0.99' = quantCol[3])) +
        scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom", 
                           legend.text = element_text(size = 12)) +
        ggtitle(paste0('Negative ', eventTypeStr,' events'))
            
      g1 <- ggplotGrob(p1)
      g2 <- ggplotGrob(p2)
      g3 <- ggplotGrob(p3)
      
      g <- rbind(g1, g2, g3, size = "first")
      g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
      pdf(file.path(figDir, paste0(eventDetector, '_', 
                                   climateVar, '_v2.pdf')))
      grid.newpage()
      grid.draw(g)
      dev.off()
      
      ## NET HISTOGRAM
      posDiff = diffEvents
      negDiff = diffEvents
      posDiff[diffEvents < 0] = 0
      negDiff[diffEvents > 0] = 0
      
      p1 = ggplot() +
        geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
        geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
        geom_point(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3]) +
        geom_point(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3]) +
        scale_x_continuous(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0(eventDetectorStr, ':\nNet ', eventTypeStr,' events'))
      
      p2 = ggplot() + geom_col(aes(x = eventYrs, y = posEvents), fill = posCol) +
        geom_line(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        geom_point(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        scale_x_continuous(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0('Positive ', eventTypeStr,' events'))
            
      p3 = ggplot() + geom_col(aes(x = eventYrs, y = negEvents), fill = negCol) +
        geom_line(aes(x = eventYrs, y = negNullQuants[1,], color = '0.9')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[2,], color = '0.95')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[3,], color = '0.99')) +
        geom_point(aes(x = eventYrs, y = negNullQuants[1,]), color = quantCol[1]) +
        geom_point(aes(x = eventYrs, y = negNullQuants[2,]), color = quantCol[2]) +
        geom_point(aes(x = eventYrs, y = negNullQuants[3,]), color = quantCol[3]) +
        scale_color_manual(name = 'Quantile',
                           values = c('0.9' = quantCol[1], 
                                      '0.95' = quantCol[2], '0.99' = quantCol[3])) +
        scale_x_continuous(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom", 
                           legend.text = element_text(size = 12)) +
        ggtitle(paste0('Negative ', eventTypeStr,' events'))
      
      
      d1 <- ggplotGrob(p1)
      d2 <- ggplotGrob(p2)
      d3 <- ggplotGrob(p3)
      
      d <- rbind(d1, d2, d3, size = "first")
      d$widths <- unit.pmax(d1$widths, d2$widths, d3$widths)
      pdf(file.path(figDir, paste0(eventDetector, '_', climateVar, '.pdf')))
      grid.newpage()
      grid.draw(d)
      dev.off()
      
      plots = list(g,d)
      
    }
    
  }
  
  if (climateVar == 'T') {
    
      posDiff_T_MS = posDiff
      negDiff_T_MS = negDiff
      posEvents_T_MS = posEvents
      negEvents_T_MS = negEvents
      quants_T_MS = diffNullQuants
      
      save(posDiff_T_MS, negDiff_T_MS, quants_T_MS,
           file = file.path(dataDir, 'histogram_T_MS.RData'))
      save(posEvents_T_MS, negEvents_T_MS, 
           file = file.path(dataDir, 'histogram_ALL_T_MS.RData'))
      
      hist_plot <- list(posDiff_T_MS, negDiff_T_MS, quants_T_MS,
                        posEvents_T_MS, negEvents_T_MS)
      
  } 
  
  if (climateVar == 'M'){
    
      posDiff_M_MS = posDiff
      negDiff_M_MS = negDiff
      posEvents_M_MS = posEvents
      negEvents_M_MS = negEvents
      quants_M_MS = diffNullQuants
      
      save(posDiff_M_MS, negDiff_M_MS, quants_M_MS,
           file = file.path(dataDir, 'histogram_M_MS.RData'))
      save(posEvents_M_MS, negEvents_M_MS, file = 
             file.path(dataDir, 'histogram_ALL_M_MS.RData'))
      
      hist_plot  <- list(posDiff_M_MS, negDiff_M_MS, quants_M_MS,
                         posEvents_M_MS, negEvents_M_MS)    
      
    }
    
  save(recordCounts, 
       file = file.path(dataDir, 
                        paste0('recordCountStats_', 
                               eventDetector, '_', 
                               climateVar, '.RData')))
  
  output  <- list(plots,hist_plot)
  
  return(output)
    
  }
  
  
