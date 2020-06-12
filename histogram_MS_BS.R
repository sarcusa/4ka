library(lipdR)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
setwd('/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/4ka')
source('createPaths.R')

figDir = '/Users/hannah/Documents/Arctic Group/Proxy analysis/forGitHub/histograms_alt/'
figDir = file.path(createPaths(), 'histograms')
datDir = file.path(createPaths(), 'RData')

# parameters
eventYrs = seq(1000,10600,by = 400)
#eventYrs = seq(900,10600,by = 400)
eventWindow = 199.9      # total event window = eventYr +/- eventWindow
numIt = 1000           # number of iterations in null model
climateVar = 'M'       # M for moisture, T for temperature, All for all
eventDetector = 'BS'   # MS for mean shift, BS for broken stick
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

# Load in relevant results
if (eventDetector == 'MS') {
  load(file.path(datDir, 'MS_results_plusNull_complete.RData'))
  TS_MS_orig = TS_MS
} else {
  # broken stick
  load(file.path(datDir, 'BS_results_plusNull_complete.RData'))
  TS_BS_orig = TS_BS
}

for (y in 1:length(eventYrs)) {
  
  eventYr = eventYrs[y]
  print(paste('Year:', eventYr))
  
  ## MEAN SHIFT PROCESSING
  if (eventDetector == 'MS') {
    
    TS_MS = TS_MS_orig
    
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
    TS = TS_MS
  } else { ## BROKEN STICK PROCESSING
    TS_BS = TS_BS_orig
    
    for (i in 1:length(TS_BS)) {
      # Filter records that don't contain event year in age range
      if (min(TS_BS[[i]]$age) > eventYr || max(TS_BS[[i]]$age) < eventYr) {
        TS_BS[[i]]$useBS = 0
        print('Out of range')
      }
      
      # Only include temp12k records for temperature analysis
      if (climateVar == 'T') {
        # if the flag is defined, but not temp12k, exclude
        if (length(TS_BS[[i]]$paleoData_inCompilation) > 0) {
          if (tolower(TS_BS[[i]]$paleoData_inCompilation) != 'temp12k') {
            TS_BS[[i]]$useBS = -1
          }
        } else {
          # if the flag isn't defined, exclude
          TS_BS[[i]]$useBS = -1
        }
      }
      
      # Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
      if (length(TS_BS[[i]]$interpretation1_seasonalityGeneral) > 0) {
        if (tolower(TS_BS[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
            tolower(TS_BS[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
          TS_BS[[i]]$useBS = -1
          print(TS_BS[[i]]$interpretation1_seasonalityGeneral)
        }
      }
    }
    TS_BS = filterTs(TS_BS, 'useBS == 1')
    
    # Assign event occurrence
    for (i in 1:length(TS_BS)) {
      
      TS_BS[[i]]$eventBS = 0
      TS_BS[[i]]$dirBS = 0
      if (!is.na(TS_BS[[i]]$brk_pts) && any(TS_BS[[i]]$brk_pts >= eventYr - eventWindow & TS_BS[[i]]$brk_pts <= eventYr + eventWindow)) {
        eve_i = which(TS_BS[[i]]$brk_pts >= eventYr - eventWindow & TS_BS[[i]]$brk_pts <= eventYr + eventWindow)
        
        # Metric combining events occurring in event window (event takes the sign of the sum)
        dir = sum(TS_BS[[i]]$brk_dirs[eve_i]) / abs(sum(TS_BS[[i]]$brk_dirs[eve_i]))
        TS_BS[[i]]$eventBS = ifelse(is.na(dir), 0, 1)
        TS_BS[[i]]$dirBS = ifelse(is.na(dir), 0, dir)
      }
      
    }
    TS = TS_BS
  }
  
  # isolate only the records corresponding to the chosen climate interpretation
  interps = unlist(sapply(TS,"[[","interpretation1_variable"))
  if (climateVar == 'M') {
    inds = which(interps == 'M' | interps == 'P')
  } else if (climateVar == 'T') {
    inds = which(interps == 'T' | interps == 'TM')
  } else {
    inds = 1:length(interps)
  }
  
  
  if (eventDetector == 'MS') {
    events = as.numeric(sapply(TS,"[[","eventMS"))[inds]
  } else {
    events = as.numeric(sapply(TS,"[[","eventBS"))[inds]
  }
  interps = interps[inds]
  
  # calculate the climate event direction based on the proxy climate dir and event dir
  dirs = unlist(sapply(TS,"[[","interpretation1_interpDirection"))[inds]
  dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
  dirs[dirs == 'negative'] = -1
  dirs[dirs == 'NA' | is.na(dirs)] = 0
  dirs = as.numeric(dirs)
  if (eventDetector == 'MS') {
    dirEvents = unlist(sapply(TS,"[[","dirMS"))[inds] # (0, 1, -1): (no, positive, negative) event
  } else {
    dirEvents = unlist(sapply(TS,"[[","dirBS"))[inds]
  }
  dirChange = dirs * dirEvents                          # (0, 1, -1): (no, positive, negative) climate event
  
  # store real events summary for the year
  allEvents[y] = sum(events == 1) / length(events)
  posEvents[y] = sum(dirChange == 1) / length(events)
  negEvents[y] = sum(dirChange == -1) / length(events)
  diffEvents[y] = posEvents[y] - negEvents[y]
  recordCounts[y,2:4] = c(length(events), posEvents[y], negEvents[y])
  
  # Assign event occurrence
  totNullEvents = matrix(0, nrow = length(inds), ncol = numIt) 
  for (i in 1:length(inds)) {
    
    if (eventDetector == 'MS') {
      nullBreaks = TS[[inds[i]]]$null_sig_brks
    } else {
      nullBreaks = TS[[inds[i]]]$null_brk_pts
    }
    nullDirs = TS[[inds[i]]]$null_brk_dirs
    
    for (j in 1:numIt) {
      
      eventInd = which(nullBreaks[[j]] >= eventYr - eventWindow & nullBreaks[[j]] <= eventYr + eventWindow)
      
      if (length(eventInd) > 0) {
        
        if (eventDetector == 'MS') {
          totNullEvents[i,j] = nullDirs[[j]][eventInd[1]] * dirs[i]
        } else {
          # Metric combining events occurring in event window (event takes the sign of the sum)
          dir = sum(nullDirs[[j]][eventInd]) / abs(sum(nullDirs[[j]][eventInd]))
          totNullEvents[i,j] = ifelse(is.na(dir), 0, dir) * dirs[i]
        }
        
      }
      
    } # end it loop
    
  } # end record loop
  
  allNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x != 0)) / length(inds)
  posNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == 1)) / length(inds)
  negNullEvents[y,] = apply(totNullEvents, 2, function(x) sum(x == -1)) / length(inds)
  
} # end event year loop

allNullQuants = apply(allNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
posNullQuants = apply(posNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
negNullQuants = apply(negNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
diffNullQuants = apply(posNullEvents - negNullEvents, 1, function(x) quantile(x, probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99)))

## HISTOGRAMS
eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))
eventDetectorStr = ifelse(eventDetector == 'MS', 'Mean shift', 'Broken stick')

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
      ggtitle(paste0(eventDetectorStr, ':\nAll events'))
    ggsave(paste0(figDir, eventDetector, '_', climateVar, '.pdf'))
    
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
       ggtitle(paste0(eventDetectorStr, ':\nAll ', eventTypeStr,' events'))
    
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
     pdf(file.path(figDir, paste0(eventDetector, '_', climateVar, '_v2.pdf')))
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
      ggtitle(paste0(eventDetectorStr, ':\nNet ', eventTypeStr,' events'))

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
    pdf(file.path(figDir, paste0(eventDetector, '_', climateVar, '.pdf')))
    grid.newpage()
    grid.draw(g)
    dev.off()
    
  }
  
} # end plotting if statement

if (climateVar == 'T') {
  
  if (eventDetector == 'MS') {
    posDiff_T_MS = posDiff
    negDiff_T_MS = negDiff
    posEvents_T_MS = posEvents
    negEvents_T_MS = negEvents
    quants_T_MS = diffNullQuants
    
    save(posDiff_T_MS, negDiff_T_MS, quants_T_MS, file = file.path(datDir, 'histogram_T_MS.RData'))
    save(posEvents_T_MS, negEvents_T_MS, file = file.path(datDir, 'histogram_ALL_T_MS.RData'))
    
  } else {
    posDiff_T_BS = posDiff
    negDiff_T_BS = negDiff
    posEvents_T_BS = posEvents
    negEvents_T_BS = negEvents
    quants_T_BS = diffNullQuants
    
    save(posDiff_T_BS, negDiff_T_BS, quants_T_BS, file = file.path(datDir, 'histogram_T_BS.RData'))
    save(posEvents_T_BS, negEvents_T_BS, file = file.path(datDir, 'histogram_ALL_T_BS.RData'))
  }

} else {
  
  if (eventDetector == 'MS') {
    posDiff_M_MS = posDiff
    negDiff_M_MS = negDiff
    posEvents_M_MS = posEvents
    negEvents_M_MS = negEvents
    quants_M_MS = diffNullQuants
    
    save(posDiff_M_MS, negDiff_M_MS, quants_M_MS, file = file.path(datDir, 'histogram_M_MS.RData'))
    save(posEvents_M_MS, negEvents_M_MS, file = file.path(datDir, 'histogram_ALL_M_MS.RData'))
    
  } else {
    posDiff_M_BS = posDiff
    negDiff_M_BS = negDiff
    posEvents_M_BS = posEvents
    negEvents_M_BS = negEvents
    quants_M_BS = diffNullQuants
    
    save(posDiff_M_BS, negDiff_M_BS, quants_M_BS, file = file.path(datDir, 'histogram_M_BS.RData'))
    save(posEvents_M_BS, negEvents_M_BS, file = file.path(datDir, 'histogram_ALL_M_BS.RData'))
  }

}
save(recordCounts, file = file.path(datDir, paste0('recordCountStats_', eventDetector, '_', climateVar, '.RData')))

