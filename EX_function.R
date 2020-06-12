EX_fun = function(age, vals, event_yr, event_window, ref_window, plotOpt = F, figDir = '',
                  resCriteria = 50, sigNum = 2, datNam = '', varNam = '', lat = NA, lon = NA, 
                  units = '', proxy = '') {
  ## Written by Hannah Kolus, 09/04/2018 
  ## Determines whether an excursion event has occurred within the specified event window.
  ## Excursion events are defined as two consecutive values within the event window that
  ## are more extreme than the avg +/- X std of the reference windows.
  
  # yr_start:yr_end defines boundaries of analysis (i.e. both reference windows and the event window)
  yr_start = event_yr - event_window / 2 - ref_window
  yr_end = event_yr + event_window / 2 + ref_window
  
  # event_start:event_end defines the boundaries of the event
  event_start = event_yr - event_window / 2
  event_end = event_yr + event_window / 2
  
  analysis_i = which(age >= yr_start & age <= yr_end) # define analysis window indices
  
  age = age[analysis_i]
  vals = vals[analysis_i]
  
  statusEX = 1      # store whether record was used in analysis or filtered
  eventEX = FALSE   # store event occurrence
  dirEX = 0         # store direction of excursion
  onsetEX = NA      # store onset age of excursion
  
  if (class(vals) != 'numeric' & class(vals) != 'integer') {
    print('NON-NUMERIC DATA')
    print(paste('DATA IS OF CLASS:', class(vals)))
    statusEX = 0
    return(list(statusEX, eventEX, dirEX, onsetEX))
  }
  
  # Calculate median resolution if there are data points in analysis window
  if (length(analysis_i) > 1) {
    medRes = median(diff(age))
  } else {
    statusEX = 0
    return(list(statusEX, eventEX, dirEX, onsetEX))
  }
  
  # Skip record if resolution doesn't meet criteria
  if (medRes > resCriteria) {
    statusEX = 0
    return(list(statusEX, eventEX, dirEX, onsetEX))
  }
  
  pre_i = which(age < event_start)                        # define pre-event (ref) window indices
  event_i = which(age >= event_start & age <= event_end)  # define event window indices
  post_i = which(age > event_end)                         # define post-event (ref) window indices
  
  # Skip records that don't contain enough data points to perform statistics
  if (length(event_i) < 5 | length(post_i) <= length(event_i)/2 | length(pre_i) <= length(event_i)/2) {
    statusEX = 0
    return(list(statusEX, eventEX, dirEX, onsetEX))
  }
  ## ------------------------------ END RECORD FILTERING ------------------------------ ##
  
  # Detrend over analysis window
  a = predict(lm(vals ~ age))
  values = as.vector(vals - a)
  
  # Calculate the avg and sd for pre-event ref window, excluding the most extreme data point
  preAVG = mean(values[pre_i])
  extremeInd = which(max(abs(preAVG - values[pre_i])) == abs(preAVG - values[pre_i]))
  preAVG = mean(values[pre_i[-extremeInd[1]]])
  preSD = sd(values[pre_i[-extremeInd[1]]]) 
  
  # Calculate the avg and sd for pre-event ref window, excluding the most extreme data point
  postAVG = mean(values[post_i])
  extremeInd = which(max(abs(postAVG - values[post_i])) == abs(postAVG - values[post_i]))
  postAVG = mean(values[post_i[-extremeInd[1]]])
  postSD = sd(values[post_i[-extremeInd[1]]]) 
  
  # Identify mean and sd used to test the high/positive anomaly threshold 
  if (preAVG + sigNum * preSD > postAVG + sigNum * postSD) {
    sd_hi = preSD
    avg_hi = preAVG
  } else {
    sd_hi = postSD
    avg_hi = postAVG
  }
  
  # Identify mean and sd used to test the low/negative anomaly threshold 
  if (preAVG - sigNum * preSD < postAVG - sigNum * postSD) {
    sd_lo = preSD
    avg_lo = preAVG
  } else {
    sd_lo = postSD
    avg_lo = postAVG 
  }
  
  # Identify points in event window that exceed the thresholds defined above
  abovePts = which(values[event_i] > avg_hi + sigNum * sd_hi)
  belowPts = which(values[event_i] < avg_lo - sigNum * sd_lo)
  
  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (medRes >= 10) {
    if (any(diff(abovePts) == 1) | any(diff(belowPts) == 1)) {
      eventEX = TRUE
      
      # Check direction
      indsA = which(diff(abovePts) == 1)
      indsB = which(diff(belowPts) == 1)
      
      if (length(indsA) > 0 && length(indsB) > 0) {
        diffA = max(values[event_i[abovePts[c(indsA, indsA+1)]]] - (avg_hi + sigNum * sd_hi))
        diffB = max((avg_lo - sigNum * sd_lo) - values[event_i[belowPts[c(indsB, indsB+1)]]])
        dirEX = ifelse(diffA > diffB, 1, -1)
      } else {
        dirEX = ifelse(length(indsA) > 0, 1, -1)
      }
      
      # find event onset
      onsetEX = ifelse(dirEX == 1, age[event_i[abovePts[indsA[length(indsA)]+1]]], age[event_i[belowPts[indsB[length(indsB)]+1]]])
      
    }
  } else {
    if (any(diff(diff(abovePts)) == 0) | any(diff(diff(belowPts)) == 0)) {
      eventEX = TRUE
      
      # Check direction
      indsA = which(diff(diff(abovePts)) == 0)
      indsB = which(diff(diff(belowPts)) == 0)
      
      if (length(indsA) > 0 && length(indsB) > 0) {
        diffA = max(values[event_i[abovePts[c(indsA, indsA+1, indsA+2)]]] - (avg_hi + sigNum * sd_hi))
        diffB = max((avg_lo - sigNum * sd_lo) - values[event_i[belowPts[c(indsB, indsB+1, indsB+2)]]])
        dirEX = ifelse(diffA > diffB, 1, -1)
      } else {
        dirEX = ifelse(length(indsA) > 0, 1, -1)
      }
      
      # find event onset
      onsetEX = ifelse(dirEX == 1, age[event_i[abovePts[indsA[length(indsA)]+2]]], age[event_i[belowPts[indsB[length(indsB)]+2]]])
    }
  }
  
  if (plotOpt) {
    
    varNam = gsub('/','_',varNam)
    varNam = gsub(',','_',varNam)
    datNam = gsub('/','_',datNam)
    datNam = gsub(',','_',datNam)
    dirStr = ifelse(dirEX == 0, 'NA', ifelse(dirEX == 1, 'positive', 'negative'))
    lat = as.numeric(lat)
    lon = as.numeric(lon)
    pdf(file.path(figDir, paste0(datNam, '_', varNam, '_', eventEX, '_', dirStr, '.pdf')))
    
    p1 = ggplot() + geom_point(aes(x = age, y = values)) + geom_line(aes(x = age, y = values)) +
      geom_line(aes(x = age[pre_i], y = rep(preAVG, length(pre_i))), linetype = 1) +
      geom_line(aes(x = age[post_i], y = rep(postAVG, length(post_i))), linetype = 1) +
      geom_line(aes(x = age[pre_i], y = rep(preAVG + sigNum * preSD, length(pre_i))), color = 'red', linetype = 4) +
      geom_line(aes(x = age[pre_i], y = rep(preAVG - sigNum * preSD, length(pre_i))), color = 'red', linetype = 4) +
      geom_line(aes(x = age[post_i], y = rep(postAVG + sigNum * postSD, length(post_i))), color = 'red', linetype = 4) +
      geom_line(aes(x = age[post_i], y = rep(postAVG - sigNum * postSD, length(post_i))), color = 'red', linetype = 4) +
      ylab(paste0(proxy, ':\n', varNam, ' (', units, ')')) +
      xlab('yr BP') + ggtitle(paste0(datNam))
    
    latLims = ifelse(c(lat >= 0,lat >= 0), c(0,90), c(-90,0))
    
    colStr = ifelse(eventEX == T, ifelse(dirEX == 1, 'red', 'blue'), 'gray50')
    p2 = ggplot() + borders("world", colour = "black") + 
      geom_point(aes(x = lon, y = lat), color = 'lightgoldenrod1', size = 5) +
      geom_point(aes(x = lon, y = lat), color=colStr, size = 4) +
      theme_bw() + xlab('Longitude') + ylab('Latitude') +
      xlim(-180, 180) + ylim(latLims)
    
    grid.arrange(p1, p2, nrow = 2)
    
    dev.off()
    #ggsave(file.path(figDir, paste0(datNam, '_', varNam, '_', eventEX, '_', dirStr, '.pdf')))
    
  }
  
  return(list(statusEX, eventEX, dirEX, onsetEX))
  
}