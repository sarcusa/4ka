########
# Stability test
# 23-11-20
# This needs to be run after running 10 small changes (A to J folders)

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



# Data are stable when run this way, but not when run in full mode -- why?
# Because Drake freeze randomness. Need to change the set_parameters enough for Drake to register that something needs to be updated, i.e. can't just rely on random numbers (also seed will not work)

name <- c("A","B","C","D","E","F","G","H","I","J")
r <- sample(name, 10)

stab <- list()

# upload the null results. Not useful for the following analysis.
for(i in 1:length(name)){
  
  n <- name[i]
  
  load(paste0("/projects/pd_lab/sha59/4ka/Stability_test/",n,"/RData/MS_results_plusNull_complete.RData"))
  
  stab[[i]]  <- data_MS
  
}



# How many iterations until null results are stable?


eventDetector == 'MS'
stab_res <- list()

# Upload the results from the analyzes.
for(k in 1:length(name)){
  
  n <- r[k]
  
  load(paste0("/projects/pd_lab/sha59/4ka/Stability_test/",n,"/RData/results_T_MS.RData"))
  
  stab_res[[k]]  <- output 
  
  
}

save(stab_res, file = "/projects/pd_lab/sha59/4ka/Stability_test/stability_test.RData")

###### Once the above code is run

load("/projects/pd_lab/sha59/4ka/Stability_test/stability_test.RData")
load("/projects/pd_lab/sha59/4ka/Stability_test/B/RData/results_T_MS.RData")
stab_res[[2]] <- output
load("/projects/pd_lab/sha59/4ka/Stability_test/C/RData/results_T_MS.RData")
stab_res[[3]] <- output

eventTypeStr = ifelse(climateVar == 'M', 'moisture', ifelse(climateVar == 'T', 'temperature', ''))
eventDetectorStr = ifelse(eventDetector == 'MS', 'Mean shift', 'Broken stick')

posCol = ifelse(climateVar == 'M', '#003c30', '#67001f')
negCol = ifelse(climateVar == 'M', '#543005', '#053061')
posFill = ifelse(climateVar == 'M', '#35978f', '#d6604d')
negFill = ifelse(climateVar == 'M', '#bf812d', '#4393c3')
quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')

eventYrs = param$eventYrs[1:25]

m  <- 2
posDiff = stab_res[[m]]$diffEvents
negDiff = stab_res[[m]]$diffEvents
posDiff[stab_res[[m]]$diffEvents < 0] = 0
negDiff[stab_res[[m]]$diffEvents > 0] = 0
diffNullQuants <- stab_res[[m]]$diffNullQuants

ggplot() +
  geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = eventYrs, y = negDiff), fill = negCol)

ggplot() +
  geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[1,]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[2,]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[3,]), color = quantCol[3]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[4,]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[5,]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[6,]), color = quantCol[3])   

# Checks for stability of the data results first, then check for the stability of number of iterations

DE_all <- sapply(stab_res,"[[","diffEvents")
PNE_all <- sapply(stab_res, "[[", "posNullEvents")
NNE_all <- sapply(stab_res, "[[", "negNullEvents")

# Measure of significant events in each run

sig <- matrix(NA,nrow = length(eventYrs), ncol = length(name))

for(o in 1:length(name)){
  
  n <- seq(1,length(eventYrs)*1000,1000)
  m <- seq(1000,length(eventYrs)*1000,1000)
    
  # Prepapre the data
  DE <- DE_all[,o]
  posDiff = DE
  negDiff = DE
  posDiff[DE < 0] = 0
  negDiff[DE > 0] = 0
  
  #Take difference of the null pos and neg events
  DN <- PNE_all[,o] - NNE_all[,o]
  DNQ <- matrix(data = NA, nrow = 6, ncol = length(eventYrs))
  
  for(i in 1:length(eventYrs)){
    
    
  #Calculate the quantiles per event year
  DNQ[, i ] <- quantile(DN[n[i]:m[i]], 
                        probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99))
    
    
    #DNQ[i, o ]  <- apply(DN[n[o]:m[o]], 2, function(x) quantile(x, probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99)))
      
    if(DE[i] > DNQ[6,i] | DE[i] < DNQ[3,i]){
      
      sig[i,o] = 1
      
    }else{
      
      sig[i,o] = 0
      
    }
    
  }
    
}

tot <- colSums(sig)

plot(1:10, tot, ylab = "Total # events (1-10.6 ka)", xlab = "# iterations")

# Increasing iterations

inc <- seq(1,10,1)
siginc <- matrix(NA,nrow = length(eventYrs), ncol = length(name))

for(o in 1:length(name)){
  
  print(paste0("run name: ",o))
  n <- seq(1,length(eventYrs)*1000,1000)
  m <- seq(1000,length(eventYrs)*1000,1000)
  
  # Prepapre the data
  DE <- DE_all[,o]
  posDiff = DE
  negDiff = DE
  posDiff[DE < 0] = 0
  negDiff[DE > 0] = 0
  
  for(i in 1:length(eventYrs)){
    for(s in inc){
      print(paste0("it num: ",s))
      
      #Take difference of the null pos and neg events for x iterations
      DN <- as.vector(PNE_all[n[i]:m[i],1:s]) - as.vector(NNE_all[n[i]:m[i],1:s])
      
      DNQ <- matrix(data = NA, nrow = 6, ncol = length(eventYrs))
      
      print(paste0("event year: ",i))
      #Calclate the quantiles per event year
      DNQ[, i ] <- quantile(DN,probs = c(0.1,0.05,0.01,0.9, 0.95, 0.99))
      
     
      if(DE[i] > DNQ[6,i] | DE[i] < DNQ[3,i]){
        
        siginc[i,o] = 1
        
      }else{
        
        siginc[i,o] = 0
        
      }
      
    }
    
  }
  
}
totinc <- colSums(siginc)

plot(1:10, totinc, ylab = "Total # events (1-10.6 ka)", xlab = "cumulative # iterations x1000")








ggplot() +
  geom_col(aes(x = eventYrs, y = posDiff), fill = posCol) +
  geom_col(aes(x = eventYrs, y = negDiff), fill = negCol) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[1]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[2]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[3]), color = quantCol[3]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[4]), color = quantCol[1]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[5]), color = quantCol[2]) +
  geom_line(aes(x = eventYrs, y = diffNullQuants[6]), color = quantCol[3]) 

