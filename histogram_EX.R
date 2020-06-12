library(lipdR)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('createPaths.R')

figDir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/histograms/'
figDir = file.path(createPaths(), 'histograms')
datDir = file.path(createPaths(), 'RData')

# parameters
eventYrs = seq(1000,10600,by = 400)
numIt = 1000           # number of iterations in null model
climateVar = 'T'       # M for moisture, T for temperature, All for all
doPlot = T

allNullEvents = matrix(NA, nrow = length(eventYrs), ncol = numIt)
posNullEvents = matrix(NA, nrow = length(eventYrs), ncol = numIt)
negNullEvents = matrix(NA, nrow = length(eventYrs), ncol = numIt)
diffNullEvents = matrix(NA, nrow = length(eventYrs), ncol = numIt)
allEvents = rep(NA, length(eventYrs))
posEvents = rep(NA, length(eventYrs))
negEvents = rep(NA, length(eventYrs))
diffEvents = rep(NA, length(eventYrs))

recordCounts = matrix(NA, nrow = length(eventYrs), ncol = 4) # years, all records, + events, - events
recordCounts[,1] = eventYrs

for (y in 1:length(eventYrs)) {
  
  eventYr = eventYrs[y]
  print(paste('Year:', eventYr))
  
  load(file.path(datDir, paste0('EX_results_plusNull_complete_', eventYr/1000, '.RData')))
  
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
  TS = filterTs(TS_EX, 'useEX == 1')
    
  # isolate only the records corresponding to the chosen climate interpretation
  interps = unlist(sapply(TS,"[[","interpretation1_variable"))
  if (climateVar == 'M') {
    inds = which(interps == 'M' | interps == 'P')
  } else if (climateVar == 'T') {
    inds = which(interps == 'T' | interps == 'TM')
  } else {
    inds = 1:length(interps)
  }
  
  events = as.numeric(sapply(TS,"[[","eventEX"))[inds]
  interps = interps[inds]
  
  # add these for printing table of results by site
  names = unlist(sapply(TS, '[[', 'dataSetName'))[inds]
  lat = unlist(sapply(TS, '[[', 'geo_latitude'))[inds]
  lon = unlist(sapply(TS, '[[', 'geo_longitude'))[inds]
  for (i in 1:length(TS)) {
    if (length(TS[[i]]$pub1_doi) == 0) {
      TS[[i]]$pub1_doi = NA
    }
  }
  doi = unlist(sapply(TS, '[[', 'pub1_doi'))[inds]
  
  # calculate the climate event direction based on the proxy climate dir and event dir
  dirs = unlist(sapply(TS,"[[","interpretation1_interpDirection"))[inds]
  dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
  dirs[dirs == 'negative'] = -1
  dirs[dirs == 'NA' | is.na(dirs)] = 0
  dirs = as.numeric(dirs)
  dirEvents = unlist(sapply(TS,"[[","dirEx"))[inds]  # (0, 1, -1): (no, positive, negative) event
  dirChange = dirs * dirEvents                       # (0, 1, -1): (no, positive, negative) climate event
  
  # save table of results by site
  event_by_site_df = data.frame(Site = names, Lat = lat, Lon = lon, Event = dirChange,
                                Year = rep(eventYr, length(names)), Var = rep(climateVar, length(names)),
                                DOI = doi)
  write.table(event_by_site_df, file = file.path(datDir, 'event_each_site_EX.csv'), append = T,
              row.names = F, col.names = F, sep = ',')
  
  # store real events summary for the year
  allEvents[y] = sum(events == 1) / length(events)
  posEvents[y] = sum(dirChange == 1) / length(events)
  negEvents[y] = sum(dirChange == -1) / length(events)
  diffEvents[y] = posEvents[y] - negEvents[y]
  recordCounts[y,2:4] = c(length(events), posEvents[y], negEvents[y])
  
  # Assign event occurrence
  totNullEvents = matrix(0, nrow = length(inds), ncol = numIt) 
  for (i in 1:length(inds)) {
    
    totNullEvents[i,] = TS[[inds[i]]]$null_events_dir * dirs[i]
    
  }
  
  allNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x != 0)) / length(inds)
  posNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == 1)) / length(inds)
  negNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == -1)) / length(inds)
  
} # end event year loop

allNullQuants = apply(allNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
posNullQuants = apply(posNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
negNullQuants = apply(negNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
diffNullQuants = apply(posNullEvents - negNullEvents, 1, function(x) quantile(x, probs = c(0.1, 0.05, 0.01, 0.9, 0.95, 0.99)))

## HISTOGRAMS
eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))

posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
negCol = ifelse(climateVar == 'M', '#543005', '#053061')
posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')

if (doPlot) {
  
  if (climateVar == 'All') {
    
    ggplot() + geom_col(aes(x = eventYrs, y = allEvents), fill = 'grey60', color = 'grey10') +
      geom_line(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
      geom_line(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
      geom_line(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
      geom_point(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
      geom_point(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
      geom_point(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom", legend.text = element_text(size = 12)) +
      ggtitle(paste0('Excursion:\nAll events'))
    ggsave(paste0(figDir, 'EX_', climateVar, '.pdf'))
    
  } else {
    
    p1 = ggplot() + geom_col(aes(x = eventYrs, y = allEvents), fill = 'grey60', color = 'grey10') +
      geom_line(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
      geom_line(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
      geom_line(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
      geom_point(aes(x = eventYrs, y = allNullQuants[1,]), color = quantCol[1]) +
      geom_point(aes(x = eventYrs, y = allNullQuants[2,]), color = quantCol[2]) +
      geom_point(aes(x = eventYrs, y = allNullQuants[3,]), color = quantCol[3]) +
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0('Excursion:\nAll ', eventTypeStr,' events'))
    
    p2 = ggplot() + geom_col(aes(x = eventYrs, y = posEvents), fill = posCol) +
      geom_line(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
      geom_line(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
      geom_line(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
      geom_point(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
      geom_point(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
      geom_point(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0('Positive ', eventTypeStr,' events'))
    #ggtitle(paste0(eventDetectorStr, ': positive ', eventTypeStr,' events'))
    
    p3 = ggplot() + geom_col(aes(x = eventYrs, y = negEvents), fill = negCol) +
      geom_line(aes(x = eventYrs, y = negNullQuants[1,], color = '0.9')) +
      geom_line(aes(x = eventYrs, y = negNullQuants[2,], color = '0.95')) +
      geom_line(aes(x = eventYrs, y = negNullQuants[3,], color = '0.99')) +
      geom_point(aes(x = eventYrs, y = negNullQuants[1,]), color = quantCol[1]) +
      geom_point(aes(x = eventYrs, y = negNullQuants[2,]), color = quantCol[2]) +
      geom_point(aes(x = eventYrs, y = negNullQuants[3,]), color = quantCol[3]) +
      scale_color_manual(name = 'Quantile', values = c('0.9' = quantCol[1], '0.95' = quantCol[2], '0.99' = quantCol[3])) +
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom", legend.text = element_text(size = 12)) +
      ggtitle(paste0('Negative ', eventTypeStr,' events'))
    #ggtitle(paste0(eventDetectorStr, ': negative ', eventTypeStr,' events'))
    
    g1 <- ggplotGrob(p1)
    g2 <- ggplotGrob(p2)
    g3 <- ggplotGrob(p3)
    #colnames(g1) <- paste0(seq_len(ncol(g1)))
    #colnames(g2) <- paste0(seq_len(ncol(g2)))
    #colnames(g3) <- paste0(seq_len(ncol(g3)))
    #grid.draw(combine(g1, g2, along=2))
    g <- rbind(g1, g2, g3, size = "first")
    g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
    pdf(file.path(figDir, paste0('EX_', climateVar, '_v2.pdf')))
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
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0('Excursion:\nNet ', eventTypeStr,' events'))
    
    p2 = ggplot() + geom_col(aes(x = eventYrs, y = posEvents), fill = posCol) +
      geom_line(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
      geom_line(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
      geom_line(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
      geom_point(aes(x = eventYrs, y = posNullQuants[1,]), color = quantCol[1]) +
      geom_point(aes(x = eventYrs, y = posNullQuants[2,]), color = quantCol[2]) +
      geom_point(aes(x = eventYrs, y = posNullQuants[3,]), color = quantCol[3]) +
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0('Positive ', eventTypeStr,' events'))
    #ggtitle(paste0(eventDetectorStr, ': positive ', eventTypeStr,' events'))
    
    p3 = ggplot() + geom_col(aes(x = eventYrs, y = negEvents), fill = negCol) +
      geom_line(aes(x = eventYrs, y = negNullQuants[1,], color = '0.1')) +
      geom_line(aes(x = eventYrs, y = negNullQuants[2,], color = '0.05')) +
      geom_line(aes(x = eventYrs, y = negNullQuants[3,], color = '0.01')) +
      geom_point(aes(x = eventYrs, y = negNullQuants[1,]), color = quantCol[1]) +
      geom_point(aes(x = eventYrs, y = negNullQuants[2,]), color = quantCol[2]) +
      geom_point(aes(x = eventYrs, y = negNullQuants[3,]), color = quantCol[3]) +
      scale_color_manual(name = 'Signficance level', values = c('0.1' = quantCol[1], '0.05' = quantCol[2], '0.01' = quantCol[3]),
                         breaks = c('0.1', '0.05', '0.01')) +
      scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], labels = eventYrs[seq(1,25,by=2)]/1000) +
      ylab('Fraction of events') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                         legend.position = "bottom", legend.text = element_text(size = 12)) +
      ggtitle(paste0('Negative ', eventTypeStr,' events'))
    
    g1 <- ggplotGrob(p1)
    g2 <- ggplotGrob(p2)
    g3 <- ggplotGrob(p3)
    #colnames(g1) <- paste0(seq_len(ncol(g1)))
    #colnames(g2) <- paste0(seq_len(ncol(g2)))
    #colnames(g3) <- paste0(seq_len(ncol(g3)))
    #grid.draw(combine(g1, g2, along=2))
    g <- rbind(g1, g2, g3, size = "first")
    g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
    pdf(file.path(figDir, paste0('EX_', climateVar, '.pdf')))
    grid.newpage()
    grid.draw(g)
    dev.off()
    
  }
  
} # end plotting if statement

if (climateVar == 'M') {
  posDiff_M_EX = posDiff
  negDiff_M_EX = negDiff
  posEvents_M_EX = posEvents
  negEvents_M_EX = negEvents
  quants_M_EX = diffNullQuants
  
  save(posDiff_M_EX, negDiff_M_EX, quants_M_EX, file = file.path(datDir, 'histogram_M_EX.RData'))
  save(posEvents_M_EX, negEvents_M_EX, file = file.path(datDir, 'histogram_ALL_M_EX.RData'))
  
} else {
  posDiff_T_EX = posDiff
  negDiff_T_EX = negDiff
  posEvents_T_EX = posEvents
  negEvents_T_EX = negEvents
  quants_T_EX = diffNullQuants
  
  save(posDiff_T_EX, negDiff_T_EX, quants_T_EX, file = file.path(datDir, 'histogram_T_EX.RData'))
  save(posEvents_T_EX, negEvents_T_EX, file = file.path(datDir, 'histogram_ALL_T_EX.RData'))
}

save(recordCounts, file = file.path(datDir, paste0('recordCountStats_EX_', climateVar, '.RData')))
