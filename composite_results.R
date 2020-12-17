#############
# Combine all parameter choices for each detector (rather than the central parameters)
# November 2020 Arcusa

dir = "/projects/pd_lab/sha59/4ka" #set the directory of the input scripts
setwd(dir)

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
EX_T  <- list()
MS_T <- list()
EX_M  <- list()
MS_M <- list()

# Mean shift

# collects results for each detector for the histogram 
for(i in 1:length(seldirsMS)){
  
  load(paste0(seldirsMS[i],"/RData/results_M_MS.RData"))
  MS_M[[i]] <- output
  
  load(paste0(seldirsMS[i],"/RData/results_T_MS.RData"))
  MS_T[[i]] <- output
  
}

#' Composites the results from the histogram script and plots the composite.
#'@description This function takes the output of the MS and EX analysis, creates a composite, and analyzes this composite as a net histogram and as tables.
#'@import cowplot
#'@param dat a list of results
#'@param param a list of parameters used in the analysis
#'@param climateVar climate variables (T or M)
#'@param cal the calculation type to use for the composite of the real events (mean or median only)
#'@param detector the detector to analyze (MS and EX)
#'@return a composite plot (a) net histogram and (b) records available through time
#'@return a list of results

compositeAndPlot <- function(dat, param, climateVar, cal, detector){
  
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
      if(posDiffind[k] > DNQ[6,k] | negDiffind[k] < DNQ[3,k]){
        
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
  
  record_num_plot  <- ggplot()+ 
    geom_line(aes(x = eventYrs, y = recordCounts),
              color = 'grey10') +
    scale_x_reverse(name = 'ky BP', 
                    breaks = eventYrs[seq(1,25,by=2)], 
                    labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('# records') +
    theme_bw() + 
    theme(legend.position = "none")
  
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
    
  out <- list(posDiff, negDiff, DNQ, sig, tot)
  return(out)
}

MS_results_T <- compositeAndPlot(dat = MS_T, param = param, climateVar = "T", cal = "mean", detector = "MS")
MS_results_M <- compositeAndPlot(dat = MS_M, param = param, climateVar = "M", cal = "mean", detector = "MS")


##### Excursion

for(i in 1:length(seldirsEX)){
  
  load(paste0(seldirsEX[i],"/RData/results_M_EX.RData"))
  EX_M[[i]] <- output
  
  load(paste0(seldirsEX[i],"/RData/results_T_EX.RData"))
  EX_T[[i]] <- output
  
}

EX_results_T <- compositeAndPlot(dat = EX_T, param = param, climateVar = "T", cal = "mean", detector = "EX")
EX_results_M <- compositeAndPlot(dat = EX_M, param = param, climateVar = "M", cal = "mean", detector = "EX")