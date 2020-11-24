histEX_plot <- function(data_in,param, climateVar){
  
  figDir = file.path(createPaths(), 'histograms')
  dataDir = file.path(createPaths(), 'RData')
  
  eventTypeStr = ifelse(climateVar == 'M', 'moisture', 
                        ifelse(climateVar == 'T', 'temperature', ''))
  posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
  negCol = ifelse(climateVar == 'M', '#543005', '#053061')
  posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
  negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
  quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')
  
  eventYrs = param$eventYrs[1:25]
  allEvents = data_in$allEvents
  allNullQuants = data_in$allNullQuants
  negEvents = data_in$negEvents
  posEvents = data_in$posEvents
  diffEvents = data_in$diffEvents
  recordCounts = data_in$recordCounts
  posNullQuants = data_in$posNullQuants
  negNullQuants = data_in$negNullQuants
  diffNullQuants = data_in$diffNullQuants
  
  
  if (plotOpt) {
    
    if (climateVar == 'All') {
      s <- ggplot() + 
        geom_col(aes(x = eventYrs, y = allEvents), 
                 fill = 'grey60', color = 'grey10')+
        geom_line(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        #geom_point(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom", 
              legend.text = element_text(size = 12)) +
        ggtitle(paste0('Excursion:\nAll events'))
      
      d  <- ggplot()+ 
        geom_line(aes(x = eventYrs, y = recordCounts[,2]),
                  color = 'grey10') +
        scale_x_reverse(name = 'ky BP', 
                        breaks = eventYrs[seq(1,25,by=2)], 
                        labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('# records') +
        theme_bw() + 
        theme(legend.position = "none")
      
      p1 <- ggplotGrob(s)
      p2 <- ggplotGrob(d)
      
      p <- rbind(p1, p2, size = "first")
      p$widths <- unit.pmax(p1$widths, p2$widths)
      pdf(file.path(figDir, paste0('EX_', climateVar, '.pdf')))
      grid.newpage()
      grid.draw(p)
      dev.off()
      
      plots = list(p)
      
        
    } else {
      
      b  <- ggplot()+ 
        geom_line(aes(x = eventYrs, y = recordCounts[,2]),
                 color = 'grey10') +
        scale_x_reverse(name = 'ky BP', 
                           breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('# records') +
        theme_bw() + 
        theme(legend.position = "none")
      
      p1 = ggplot() + geom_col(aes(x = eventYrs, y = allEvents), 
                               fill = 'grey60', color = 'grey10') +
        geom_line(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        #geom_point(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0('Excursion:\nAll ', eventTypeStr,' events'))
      
      p2 = ggplot() + geom_col(aes(x = eventYrs, y = posEvents), fill = posCol) +
        geom_line(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
       # geom_point(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0('Positive ', eventTypeStr,' events'))
      #ggtitle(paste0(eventDetectorStr, ': positive ', eventTypeStr,' events'))
      
      p3 = ggplot() + geom_col(aes(x = eventYrs, y = negEvents), fill = negCol) +
        geom_line(aes(x = eventYrs, y = negNullQuants[1,], color = '0.9')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[2,], color = '0.95')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[3,], color = '0.99')) +
        #geom_point(aes(x = eventYrs, y = negNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = negNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = negNullQuants[3,]), color = quantCol[3]) +
        scale_color_manual(name = 'Quantile', values = c('0.9' = quantCol[1], 
                                                         '0.95' = quantCol[2], 
                                                         '0.99' = quantCol[3])) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom", 
                           legend.text = element_text(size = 12)) +
        ggtitle(paste0('Negative ', eventTypeStr,' events'))
            
      g0 <- ggplotGrob(b)
      g1 <- ggplotGrob(p1)
      g2 <- ggplotGrob(p2)
      g3 <- ggplotGrob(p3)
      
      g <- rbind(g1, g2, g3, g0, size = "first")
      g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g0$widths)
      pdf(file.path(figDir, paste0('EX_', 
                                   climateVar, '_v2.pdf')))
      grid.newpage()
      grid.draw(g)
      dev.off()
      
      ## NET HISTOGRAM
      
      b  <- ggplot()+ 
        geom_line(aes(x = eventYrs, y = recordCounts[,2]),
                  color = 'grey10') +
        scale_x_reverse(name = 'ky BP', 
                        breaks = eventYrs[seq(1,25,by=2)], 
                        labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('# records') +
        theme_bw() + 
        theme(legend.position = "none")
      
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
        #geom_point(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3]) +
        #geom_point(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
       # geom_point(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
       # geom_point(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3]) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0('Excursion:\nNet ', eventTypeStr,' events'))
      
      p2 = ggplot() + geom_col(aes(x = eventYrs, y = posEvents), fill = posCol) +
        geom_line(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        geom_line(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        #geom_point(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0('Positive ', eventTypeStr,' events'))
            
      p3 = ggplot() + geom_col(aes(x = eventYrs, y = negEvents), fill = negCol) +
        geom_line(aes(x = eventYrs, y = negNullQuants[1,], color = '0.1')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[2,], color = '0.05')) +
        geom_line(aes(x = eventYrs, y = negNullQuants[3,], color = '0.01')) +
        #geom_point(aes(x = eventYrs, y = negNullQuants[1,]), color = quantCol[1]) +
        #geom_point(aes(x = eventYrs, y = negNullQuants[2,]), color = quantCol[2]) +
        #geom_point(aes(x = eventYrs, y = negNullQuants[3,]), color = quantCol[3]) +
        scale_color_manual(name = 'Signficance level', 
                           values = c('0.1' = quantCol[1], 
                                      '0.05' = quantCol[2], '0.01' = quantCol[3]),
                           breaks = c('0.1', '0.05', '0.01')) +
        scale_x_reverse(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                           labels = eventYrs[seq(1,25,by=2)]/1000) +
        ylab('Fraction of events') +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom", legend.text = element_text(size = 12)) +
        ggtitle(paste0('Negative ', eventTypeStr,' events'))
      
      d0 <- ggplotGrob(b)
      d1 <- ggplotGrob(p1)
      d2 <- ggplotGrob(p2)
      d3 <- ggplotGrob(p3)
      
      d <- rbind(d1, d2, d3, d0, size = "first")
      d$widths <- unit.pmax(d1$widths, d2$widths, d3$widths, d0$widths)
      pdf(file.path(figDir, paste0('EX_', climateVar, '.pdf')))
      grid.newpage()
      grid.draw(d)
      dev.off()
      
      plots = list(g,d)
      
      
      
    }
    
  }
  
  if (climateVar == 'M') {
    
    posDiff = diffEvents
    negDiff = diffEvents
    posDiff[diffEvents < 0] = 0
    negDiff[diffEvents > 0] = 0
    
    posDiff_M_EX = posDiff
    negDiff_M_EX = negDiff
    posEvents_M_EX = posEvents
    negEvents_M_EX = negEvents
    quants_M_EX = diffNullQuants
    
    save(posDiff_M_EX, negDiff_M_EX, quants_M_EX, 
         file = file.path(dataDir, 'histogram_M_EX.RData'))
    save(posEvents_M_EX, negEvents_M_EX, 
         file = file.path(dataDir, 'histogram_ALL_M_EX.RData'))
    
    histogram  <- list(posDiff_M_EX, negDiff_M_EX, quants_M_EX,
                       posEvents_M_EX, negEvents_M_EX)
    
  } 
  
  if (climateVar == 'T') {
    
    posDiff = diffEvents
    negDiff = diffEvents
    posDiff[diffEvents < 0] = 0
    negDiff[diffEvents > 0] = 0
    
    posDiff_T_EX = posDiff
    negDiff_T_EX = negDiff
    posEvents_T_EX = posEvents
    negEvents_T_EX = negEvents
    quants_T_EX = diffNullQuants
    
    save(posDiff_T_EX, negDiff_T_EX, quants_T_EX, 
         file = file.path(dataDir, 'histogram_T_EX.RData'))
    save(posEvents_T_EX, negEvents_T_EX, 
         file = file.path(dataDir, 'histogram_ALL_T_EX.RData'))
    
    histogram  <- list(posDiff_T_EX, negDiff_T_EX, quants_T_EX,
                       posEvents_T_EX, negEvents_T_EX)
    
  }
  
  if (climateVar == 'All') {

    histogram <- NA

  }  
     
  save(recordCounts, file = file.path(dataDir, paste0('recordCountStats_EX_', climateVar, '.RData')))
  
  output  <- list(plots, histogram)
  
  return(output)
  
}
