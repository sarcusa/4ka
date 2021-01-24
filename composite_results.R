#############
# Combine all parameter choices for each detector (rather than the central parameters)
# November 2020 Arcusa

dir = "/projects/pd_lab/sha59/4ka" #set the directory of the input scripts
setwd(dir)

outDir = paste0(dir,"/Composite/")

source('packages.R')
source('functions.R')
source('set_parameters.R')
source('plan_prep.R')
source('plan_prep_excursion.R')
source('plan_prep_meanshift.R')
source('plan_prep_trendchanges.R')
source('plan_excursion.R')
source('plan_meanshift.R')
source('plan_trendchanges.R')
#source('plan_plotting.R')
source('plan_var.R')

parameters = list(eventYrs = eventYrs, event_window = event_window, 
                  ref_window = ref_window, resCriteria = resCriteria, 
                  sigNum = sigNum, plotOpt = plotOpt, 
                  mainDir = mainDir, numIt = numIt, 
                  res = res, radius = radius, 
                  eventWindow = eventWindow, 
                  CName = CName, CVers = CVers, 
                  eventDetector = eventDetector, OutDat = OutDat,
                  ncores = ncores, maxDiff = maxDiff)

param = parameters

# lists folders in directory, and separates for each detector
listdirs <- list.dirs(path = dir, full.names = T,recursive = F)
seldirs <- listdirs[grep(pattern = "Signif", x = listdirs)]
seldirs <- seldirs[-grep(pattern = "random", seldirs)]
seldirs <- seldirs[-grep(pattern = "seq", seldirs)]
seldirsEX <- c(seldirs[grep(pattern = "EX", seldirs)])
seldirsMS <- c(seldirs[grep(pattern = "MS", seldirs)])

#For the histogram
EX_T  <- list()
MS_T <- list()
EX_M  <- list()
MS_M <- list()
#For the spatial maps
EX_T_s  <- list()
MS_T_s <- list()
EX_M_s  <- list()
MS_M_s <- list()

#' Composites the results from the histogram script and plots the composite.
#'@description This function takes the output of the MS and EX analysis, creates a composite, and analyzes this composite as a net histogram and as tables.
#'@import cowplot
#'@param dat a list of results
#'@param param a list of parameters used in the analysis
#'@param climateVar climate variables (T or M)
#'@param detector the detector to analyze (MS and EX)
#'@return a composite plot (a) net histogram and (b) records available through time
#'@return a list of results

compositeHist <- function(dat, param, climateVar, detector){
  
  figDir = file.path(dir, 'histograms')
  eventYrs = param$eventYrs[1:25] #removes event yr 11ka (not enough data)
  
  allEvents = rep(NA, length(eventYrs))
  posEvents = rep(NA, length(eventYrs))
  negEvents = rep(NA, length(eventYrs))
  diffEvents = rep(NA, length(eventYrs))
  
  allNullEvents = matrix(NA, nrow = length(eventYrs), ncol = param$numIt)
  posNullEvents = matrix(NA, nrow = length(eventYrs), ncol = param$numIt)
  negNullEvents = matrix(NA, nrow = length(eventYrs), ncol = param$numIt)
  
  recordCounts = matrix(NA, nrow = length(eventYrs), ncol = length(dat)+3) # years, all records, + events, - events
  recordCounts[,1] = eventYrs
  
  for(j in 1:length(eventYrs)){
    
    print(paste0("eventYr: ", eventYrs[j]))
    
    #Composite
    events <- unlist(sapply(dat, function(x) x$everyEvent[[j]]))
    dirChange <- unlist(sapply(dat, function(x) x$everyDirChange[[j]]))
    totNullEventsRec <- lapply(dat, function(x) x$everyTotNullEvent[[j]])
    
    if(class(totNullEventsRec) == "list"){
      
      totNullEvents <- do.call(rbind,totNullEventsRec)
      
    }else{
      
      totNullEvents <- totNullEventsRec
    }
    
    Rec <- sapply(dat, function(x) x$recordCounts[j,2])
    
    # store real events summary for the year
    allEvents[j] = sum(events != 0) / length(events)
    posEvents[j] = sum(dirChange == 1) / length(events)
    negEvents[j] = sum(dirChange == -1) / length(events)
    diffEvents[j] = posEvents[j] - negEvents[j]
    
    #Store real event summary for the year (composite and individual)
    recordCounts[j,2:ncol(recordCounts)] = c(Rec, posEvents[j], negEvents[j])
    
    allNullEvents[j,] = apply(totNullEvents, 2, function(x) sum(x != 0)) / nrow(totNullEvents)
    posNullEvents[j,] = apply(totNullEvents, 2, function(x) sum(x == 1))  / nrow(totNullEvents)
    #sum(sapply(totNullEvents, function(x) nrow(x))) used to be this
    negNullEvents[j,] = apply(totNullEvents, 2, function(x) sum(x == -1)) / nrow(totNullEvents)
    #sum(sapply(totNullEvents, function(x) nrow(x)))
    
  } #end of year loop
  
  allNullQuants = apply(allNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  posNullQuants = apply(posNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  negNullQuants = apply(negNullEvents, 1, function(x) quantile(x, probs = c(0.9, 0.95, 0.99)))
  diffNullQuants = apply(posNullEvents - negNullEvents, 1, function(x) quantile(x, probs = c(0.1, 0.05, 0.01, 0.9, 0.95, 0.99)))
  
  eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))
  
  posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
  negCol = ifelse(climateVar == 'M', '#543005', '#053061')
  posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
  negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
  quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')
  
  posDiff = diffEvents
  negDiff = diffEvents
  posDiff[diffEvents < 0] = 0
  negDiff[diffEvents > 0] = 0
  
  recordCounts <- as.data.frame(recordCounts)
  names(recordCounts)[1] <- "eventYr"
  m <- reshape2::melt(recordCounts[,1:(length(dat)+1)], "eventYr")
  
  record_num_plot  <- ggplot(m, aes(x = eventYr, y = value, group = variable))+ 
    geom_line() +
    scale_x_reverse(name = 'ky BP', 
                    breaks = eventYrs[seq(1,25,by=2)], 
                    labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('# records') +
    theme_bw() + 
    theme(legend.position = "none")
  
  netHist = ggplot() +
    geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
    geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
    geom_line(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3]) +
    scale_x_reverse(name = 'ky BP', 
                    breaks = eventYrs[seq(1,25,by=2)], 
                    labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('Fraction of events') +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(detector,':\nNet ', eventTypeStr,' events'))
  
  together <- cowplot::plot_grid(netHist, record_num_plot,
                                 align = "hv", nrow = 2,rel_heights = c(1,0.25))
  
  pdf(file.path(figDir, paste0('Composite_',detector ,'_', climateVar,'.pdf')))
  grid.draw(together)
  dev.off() 
  
  return(recordCounts)
  
}

OriginalCompositeHist <- function(dat, param, climateVar, detector, cal){
  
  figDir = file.path(dir, 'histograms')
  
  DE_all <- sapply(dat,"[[","diffEvents")
  PNE_all <- sapply(dat, "[[", "posNullEvents")
  NNE_all <- sapply(dat, "[[", "negNullEvents")
  recordCounts_list <- lapply(dat, "[[", "recordCounts")
  recordCounts_all <- sapply(recordCounts_list, function(x) x[,2])
  recordCounts <- recordCounts_all[,1]
  
  eventYrs = param$eventYrs[1:25] #removes event yr 11ka (not enough data)
  
  sig <- matrix(NA,nrow = length(eventYrs), ncol = length(dat)+1)
  
  #For the composite
  
  n <- seq(1,length(eventYrs)*1000,1000)
  m <- seq(1000,length(eventYrs)*1000,1000)
  
  # Prepapre the data
  posDiff = DE_all
  negDiff = DE_all
  
  # Temporarily set the values above and below zero to NA to composite
  posDiff[DE_all < 0] = NA
  negDiff[DE_all > 0] = NA
  
  # How should the real event composite be calculated?
  if(cal == "mean"){
    
    posDiff  <- apply(posDiff, MARGIN = 1, FUN = mean, na.rm = T)
    negDiff  <- apply(negDiff, MARGIN = 1, FUN = mean, na.rm = T)
    
  }else{
    
    posDiff  <- apply(posDiff, MARGIN = 1, FUN = median, na.rm = T)
    negDiff  <- apply(negDiff, MARGIN = 1, FUN = median, na.rm = T)
  }
  
  # Reset NA values to 0
  posDiff[posDiff == "NaN"] <- 0
  negDiff[negDiff == "NaN"] <- 0
  posDiff[is.na(posDiff)] <- 0
  negDiff[is.na(negDiff)] <- 0
  
  DNQ <- matrix(data = NA, nrow = 6, ncol = length(eventYrs))
  
  for(j in 1:length(eventYrs)){
    
    # Takes the difference of the positive and negative null events of all
    DN <- as.vector(PNE_all[n[j]:m[j],]) - as.vector(NNE_all[n[j]:m[j],])
        
    # calculates the quantiles of the composite
    DNQ[,j] <- unname(quantile(DN, 
                               probs = c(0.1,0.05,0.01,0.9,0.95,0.99)))
    
    # Objective identification of significant events from composite
    if(posDiff[j] > DNQ[6,j] | negDiff[j] < DNQ[3,j]){
      
      sig[j,1] = 1
      
    }else{
      
      sig[j,1] = 0
      
    }
    
  }  
  
  
  # For each individual analysis
  for(i in 1:length(dat)){
    
    n <- seq(1,length(eventYrs)*1000,1000)
    m <- seq(1000,length(eventYrs)*1000,1000)
    
    # Prepapre the data
    DE <- DE_all[,i]
    posDiffind = DE
    negDiffind = DE
    posDiffind[posDiffind < 0] = 0
    negDiffind[posDiffind > 0] = 0
    
    DNQind <- matrix(data = NA, nrow = 6, ncol = length(eventYrs))
    
    for(k in 1:length(eventYrs)){
      
      # Takes the difference of the positive and negative null events of all
      DN <- as.vector(PNE_all[n[k]:m[k],i]) - as.vector(NNE_all[n[k]:m[k],i])
      # calculates the quantiles 
      DNQind[,k] <- unname(quantile(DN, 
                                    probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99)))
      
      # Objective identification of significant events from composite
      if(posDiffind[k] > DNQind[6,k] | negDiffind[k] < DNQind[3,k]){
        
        sig[k,i+1] = 1
        
      }else{
        
        sig[k,i+1] = 0
        
      }
      
    }  
    
  }
  
  tot = colSums(sig) # columns are 1= composite, 2-10 = indiv analysis
  
  eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))
  
  posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
  negCol = ifelse(climateVar == 'M', '#543005', '#053061')
  posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
  negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
  quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')
  
  recordCounts_all <- as.data.frame(recordCounts_all)
  recordCounts_all <- cbind(eventYrs,recordCounts_all)
  m <- reshape2::melt(recordCounts_all[,1:(length(dat)+1)], "eventYrs")
  
  record_num_plot  <- ggplot(m, aes(x = eventYrs, y = value, group = variable))+ 
    geom_line() +
    scale_x_reverse(name = 'ky BP', 
                    breaks = eventYrs[seq(1,25,by=2)], 
                    labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('# records') +
    theme_bw() + 
    theme(legend.position = "none")  
  
  #record_num_plot  <- ggplot()+ 
  #  geom_line(aes(x = eventYrs, y = recordCounts),
  #            color = 'grey10') +
  #  scale_x_reverse(name = 'ky BP', 
  #                  breaks = eventYrs[seq(1,25,by=2)], 
  #                  labels = eventYrs[seq(1,25,by=2)]/1000) +
  #  ylab('# records') +
  #  theme_bw() + 
  #  theme(legend.position = "none")
  
  netHist = ggplot() +
    geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
    geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
    geom_line(aes(x = eventYrs, y = DNQ[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = DNQ[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = DNQ[3,]), color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = DNQ[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = DNQ[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = DNQ[6,]), color = quantCol[3]) +
    scale_x_reverse(name = 'ky BP', 
                    breaks = eventYrs[seq(1,25,by=2)], 
                    labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('Fraction of events') +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(detector,':\nNet ', eventTypeStr,' events'))
  
  together <- cowplot::plot_grid(netHist, record_num_plot,
                                 align = "hv", nrow = 2,rel_heights = c(1,0.25))
  
  pdf(file.path(figDir, paste0('Composite_',detector ,'_', climateVar,"_",cal, '.pdf')))
  grid.draw(together)
  dev.off()    
    
  out <- list(PosDiff = posDiff, NegDiff = negDiff, DNQ = DNQ, Sig = sig, Tot = tot, recordCounts_all)
  return(out)
}

#' Spatial composite maps
#'@description This function takes the results from the GridCreate function and plots the composite spatially for each event yr and each climate variable
#'@import cowplot
#'@param D a list of directories
#'@param param a list of parameters used in the analysis
#'@param climateVar climate variables (T or M)
#'@param cal the calculation type to use for the composite of the real events (mean or median only)
#'@param detector the detector to analyze (MS and EX)
#'@return a composite plot (a) net histogram and (b) records available through time
#'@return a list of results

compositeSpatial <- function(D, param, climateVar, cal, detector){
  
  msDir = file.path(createPaths(), 'mean_shift')
  
  # Make a grid of lat/lon center points
  latitude = seq(-90 + param$res/2, 90 - param$res/2, param$res)
  longitude = seq(-180 + param$res/2, 180 - param$res/2, param$res)
  
  for(k in 1:length(param$eventYrs[1:25])){
    
    eventYr = param$eventYrs[k]
    
    res <- list()
    
    n <- seq(1,length(D)*1000,1000)
    m <- seq(1000,length(D)*1000,1000)
    
    gridNullEvents = array(data = NA, dim = c(72,36,9000))
    gridNumRecs = array(data = NA, dim = c(72,36,9))
    gridEvents = array(data = NA, dim = c(72,36,9))
    totNullEvents = list()
    gridPercentEvents_null = matrix(NA, nrow = length(longitude), 
                                    ncol = length(latitude))
    allTsLon <-list()
    allTsLat  <- list()
    dirChange <- list()
    
    # collects results from the null analysis processed by GridCreate
    for(i in 1:length(D)){
      
      list2env(loading(paste0(D[i], '/RData/results_spatial_MS_', climateVar, "_", eventYr/1000, ".RData")),envir=.GlobalEnv)
      
      names(OP) <- c("gridNumRecs","gridEvents","gridNullEvents",
                     "totNullEvents","allTsLon", "allTsLat", "dirChange")
      
      gridNumRecs[,,i] <- OP$gridNumRecs
      gridNullEvents[,,n[i]:m[i]] <- OP$gridNullEvents
      gridEvents[,,i] <- OP$gridEvents
      totNullEvents[[i]] <- OP$totNullEvents
      allTsLon[[i]]  <- OP$allTsLon
      allTsLat[[i]] <- OP$allTsLat
      dirChange[[i]] <- OP$dirChange
    
    } # end of D loop
    
    # Compositing the real events
          
      gridEvents  <- apply(gridEvents, MARGIN = c(1,2), function(x) sum(x ), na.rm = T)
            
         
    # Compositing the event locations
    
    
    
    # Calculating probabilities
    for (i in 1:length(longitude)) {
      
      print(paste('iteration i = ', i))
      
      for (j in 1:length(latitude)) {
        
        if (is.na(gridEvents[i, j])) {
          next
        } else if (gridEvents[i, j] < 0) {
          gridPercentEvents_null[i,j] = -sum(gridNullEvents[i,j,] > gridEvents[i, j]) / param$numIt
        } else if (gridEvents[i, j] > 0) {
          gridPercentEvents_null[i,j] = sum(gridNullEvents[i,j,] < gridEvents[i, j]) / param$numIt
        } else {
          gridPercentEvents_null[i,j] = 0
        }
        
      }
    } #end of longitude loop
    
    # Preparing to plot
    locs = which(!is.na(gridEvents), arr.ind = T)
    locs_na = which(is.na(gridEvents), arr.ind = T)
    percentEvents_NULL = gridPercentEvents_null[locs]
    
    bins = c('<= 0.05', '0.05 - 0.1', '0.1 - 0.25', '0.25 - 0.5', '> 0.5', '0.5 - 0.25', '0.25 - 0.1', '0.1 - 0.05', ' <= 0.05')
    locs_binned = matrix(nrow = nrow(locs) + length(bins), ncol = 2)
    locs_binned[1:nrow(locs),] = locs
    percentEvents_NULL_binned = percentEvents_NULL
    percentEvents_NULL_binned[percentEvents_NULL >= 0.95] = bins[1]
    percentEvents_NULL_binned[percentEvents_NULL < 0.95 & percentEvents_NULL >= 0.9] = bins[2]
    percentEvents_NULL_binned[percentEvents_NULL < 0.9 & percentEvents_NULL >= 0.75] = bins[3]
    percentEvents_NULL_binned[percentEvents_NULL < 0.75 & percentEvents_NULL >= 0.5] = bins[4]
    percentEvents_NULL_binned[percentEvents_NULL < 0.5 & percentEvents_NULL > -0.5] = bins[5]
    percentEvents_NULL_binned[percentEvents_NULL > -0.75 & percentEvents_NULL <= -0.50] = bins[6]
    percentEvents_NULL_binned[percentEvents_NULL > -0.9 & percentEvents_NULL <= -0.75] = bins[7]
    percentEvents_NULL_binned[percentEvents_NULL > -0.95 & percentEvents_NULL <= -0.9] = bins[8]
    percentEvents_NULL_binned[percentEvents_NULL <= -0.95] = bins[9]
    
    dummy = length(percentEvents_NULL_binned) + 1
    dummy_i = dummy:(dummy+length(bins)-1)
    for (i in 1:length(bins)) {
      percentEvents_NULL_binned[dummy] = bins[i]
      locs_binned[dummy,] = locs_na[10*i,]
      dummy = dummy + 1
    }
    
    percentEvents_NULL_binned = factor(percentEvents_NULL_binned)
    percentEvents_NULL_binned = factor(percentEvents_NULL_binned,levels(percentEvents_NULL_binned)[c(2,4,6,8,3,9,7,5,1)])
    
    if (climateVar == 'M') {
      
      myCol = c('#543005','#bf812d','#dfc27d','#f6e8c3','snow2','#c7eae5','#80cdc1','#35978f','#003c30')
      
      tryCatch({s <-baseMAP(lon=0,lat = 0,projection = "mollweide",global = TRUE,map.type = "line",restrict.map.range = F, country.boundaries = F) + 
                  geom_tile(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]],width = 5.5,height = 5.5, fill = as.factor(percentEvents_NULL_binned))) + 
                  borders("world", colour="grey70") + 
                  geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
                  geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white', size = 3) +
                  geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white', size = 3) +
                  geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 2, shape = 21, fill = "grey50") +
                  geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'), size = 3, shape = 24, fill = "blue") +
                  geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'), size = 3, shape = 25, fill = "tomato4") +
                  scale_color_manual(name = '', values = c('no event' = 'black', '+ event' = 'white', '- event' = 'white'),breaks = c('+ event', '- event', 'no event'),guide = guide_legend(override.aes = list(shape = c(24, 25, 21), fill = c('blue','tomato4','grey50'),color = c('black','black','black')))) +
                  scale_fill_manual(name = '', values = rev(myCol))+
                  ggtitle(paste0('MS: Fraction of null events < real event #\n', eventYr/1000,'+/-',param$eventWindow/2/1000, 'ka events'))
                  #geom_rect(aes(xmax=180.1,xmin=-180.1,ymax=90.1,ymin=-90.1),fill=NA, colour="black")     
                  
                  pdf(file.path(msDir, 'spatial_MS_M', 
                                paste0(y,'-', eventYr/1000, '.pdf')))
                print(s)
                dev.off() }, error = function(e){cat("Error:", conditionMessage(e), " may not have plotted")})   
      
    }
    
    if (climateVar == 'T') {
      
      myCol = c('#67001f','#d6604d','#f4a582','#fddbc7','snow2','#d1e5f0','#92c5de','#4393c3','#053061')
      
      tryCatch({ s <- baseMAP(lon=0,lat = 0,projection = "mollweide",global = TRUE,map.type = "line",restrict.map.range = F, country.boundaries = F) + 
                   geom_tile(aes(x = longitude[locs_binned[,1]], y = latitude[locs_binned[,2]],width = 5.5,height = 5.5, fill = as.factor(percentEvents_NULL_binned))) + 
                   borders("world", colour="grey70") + 
                   geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)]), color='white', size = 3) +
                   geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)]), color='white', size = 3) +
                   geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)]), color='white', size = 3) +
                   geom_point(aes(x = allTsLon[which(dirChange == 0)], y = allTsLat[which(dirChange == 0)], color='no event'), size = 2, fill = "grey50", shape = 21) +
                   geom_point(aes(x = allTsLon[which(dirChange == 1)], y = allTsLat[which(dirChange == 1)], color='+ event'), size = 3, fill = "red", shape = 24) +
                   geom_point(aes(x = allTsLon[which(dirChange == -1)], y = allTsLat[which(dirChange == -1)], color='- event'), size = 3, fill = "royalblue", shape = 25) +
                   scale_color_manual(name = '', values = c('no event' = 'black', '+ event' = 'white', '- event' = 'white'),
                                      breaks = c('+ event', '- event', 'no event'),
                                      guide = guide_legend(override.aes = list(shape = c(24, 25, 21), fill = c('red','royalblue','grey50'),color = c('black','black','black')))) +
                   scale_fill_manual(name = '', values = myCol) +
                   ggtitle(paste0('MS: Fraction of null events < real event #\n', eventYr/1000,'+/-',param$eventWindow/2/1000, 'ka events'))
                   #geom_rect(aes(xmax=180.1,xmin=-180.1,ymax=90.1,ymin=-90.1),fill=NA, colour="black")
                   
                   pdf(file.path(msDir, 'spatial_MS_T', 
                                 paste0(y,'-', eventYr/1000, '.pdf')))
                 print(s)
                 dev.off()}, error = function(e){cat("Error:", conditionMessage(e), " may not have plotted")})
      
      
      
    }
    
    
  } #end of eventYr loop
  
  
}#end of function


### Mean shift

# collects results for each detector for the histogram 
for(i in 1:length(seldirsMS)){
  
  load(paste0(seldirsMS[i],"/RData/results_M_MS.RData"))
  MS_M[[i]] <- output
  
  load(paste0(seldirsMS[i],"/RData/results_T_MS.RData"))
  MS_T[[i]] <- output
  
}

MS_results_T <- compositeHist(dat = MS_T, param = param, climateVar = "T", detector = "MS")
MS_results_M <- compositeHist(dat = MS_M, param = param, climateVar = "M", detector = "MS")

MS_results_T_ori <- OriginalCompositeHist(dat = MS_T, param = param, climateVar = "T", detector = "MS", cal = "median")
MS_results_M_ori <- OriginalCompositeHist(dat = MS_M, param = param, climateVar = "M", detector = "MS", cal = "median")

# Spatial maps 

GridCreate(D = seldirsMS[1], param = param, climateVar = "M")

##### Excursion

for(i in 1:length(seldirsEX)){
  
  load(paste0(seldirsEX[i],"/RData/results_M_EX.RData"))
  EX_M[[i]] <- output
  
  load(paste0(seldirsEX[i],"/RData/results_T_EX.RData"))
  EX_T[[i]] <- output
  
}

EX_results_T <- compositeHist(dat = EX_T, param = param, climateVar = "T", detector = "EX")
EX_results_M <- compositeHist(dat = EX_M, param = param, climateVar = "M",  detector = "EX")

EX_results_T_ori <- OriginalCompositeHist(dat = EX_T, param = param, climateVar = "T", detector = "EX", cal = "median")
EX_results_M_ori <- OriginalCompositeHist(dat = EX_M, param = param, climateVar = "M", detector = "EX", cal = "median")

# Write results out

names_EX <- vector()
names_MS <- vector()

for(i in 1:length(seldirsEX)){
  
  a <- substr(seldirsEX[i], 43, nchar(seldirsEX[i]))
  names_EX[i] <- a
  
  b <- substr(seldirsMS[i], 43, nchar(seldirsMS[i]))
  names_MS[i] <- b
  
}
names_EX <- c("comp",names_EX)
names_MS <- c("comp",names_MS)

library(xlsx)
write.xlsx(MS_results_T[[4]], file=paste0(dir,"/sensitivity_results.xlsx"), sheetName="MS_T", row.names=FALSE)
write.xlsx(MS_results_M[[4]], file=paste0(dir,"/sensitivity_results.xlsx"), sheetName="MS_M", row.names=FALSE, append = T)
write.xlsx(EX_results_T[[4]], file=paste0(dir,"/sensitivity_results.xlsx"), sheetName="EX_T", row.names=FALSE, append = T)
write.xlsx(EX_results_M[[4]], file=paste0(dir,"/sensitivity_results.xlsx"), sheetName="EX_M", row.names=FALSE, append = T)
write.xlsx(c(names_EX,names_MS), file=paste0(dir,"/sensitivity_results.xlsx"), sheetName="names", row.names=FALSE, append = T)


