#######
# Abrupt climate change 
# Analysis of 4.2 ka event
# Filter and plot records that do show an event

dir = "/projects/pd_lab/sha59/4ka" #set the directory of the input scripts
setwd(dir)

source('packages.R')
source('functions.R')
source('set_parameters.R')

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
  
eventYr = 4200
 
# Excursion 
list2env(loading('/projects/pd_lab/sha59/4ka/RData/EX_results_plusNull_complete_4.2.RData'),envir=.GlobalEnv)
  
for (i in 1:length(data_EX)) {
  
  # Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
  if (length(data_EX[[i]]$interpretation1_seasonalityGeneral) > 0) {
    if (tolower(data_EX[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
          tolower(data_EX[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
      data_EX[[i]]$useEX = -1
      print(data_EX[[i]]$interpretation1_seasonalityGeneral)
    }
  }
  
  
}



# Records with and without events
# Moisture
interps <- pullTsVariable(TS = data_EX,variable = "interpretation1_variable")
index <- which(interps == 'M' | interps == 'P' | interps == 'P-E' | interps ==  'P/E')
TS_EX_M <- data_EX[index]
dirs <- pullTsVariable(TS = TS_EX_M,variable = "interpretation1_interpDirection")
dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
dirs[dirs == 'negative'] = -1
dirs[dirs == 'NA' | is.na(dirs)] = 0
dirs = as.numeric(dirs)
dirEvents = pullTsVariable(TS = TS_EX_M, variable = "dirEx")
dirChange = dirs * dirEvents 
for(i in 1:length(TS_EX_M)){
  if(dirChange[i] == 0) next
  #TS_EX_M[[i]]$paleoData_values = TS_EX_M[[i]]$paleoData_values*dirChange[i]
  TS_EX_M[[i]]$dirEx = dirChange[i]
  
}

tidy.M <- tidyTs(TS_EX_M, age.var = "age")

df.M <- tidy.M %>% 
  filter(between(age,3000,5000)) %>% 
  group_by(eventEX) %>% 
  arrange(eventEX) %>%
  arrange(dirEx)

all.M  <- plotTimeseriesStack(df.M, color.var = "eventEX", time.var = "age")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.M$dataSetName))+1,alpha = 0.2)+
  ggtitle("Moisture: Excursion")
  
ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.M_EX.pdf", plot = all.M, limitsize = F, height = 40)
dev.off()

map.all.M <- mapTs(TS_EX_M, color = 'eventEX')+
  scale_color_manual(values = c("FALSE" = "#1B9E77", "TRUE" = "#7570B3"))+
  ggtitle("Moisture: Excursion")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/map.all.M_EX.pdf", plot = map.all.M)
dev.off()

# Temperature

index <- which(interps == 'T' | interps == 'TM')
TS_EX_T <- data_EX[index]
dirs <- pullTsVariable(TS = TS_EX_T,variable = "interpretation1_interpDirection")
dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
dirs[dirs == 'negative'] = -1
dirs[dirs == 'NA' | is.na(dirs)] = 0
dirs = as.numeric(dirs)
dirEvents = pullTsVariable(TS = TS_EX_T, variable = "dirEx")
dirChange = dirs * dirEvents 
for(i in 1:length(TS_EX_T)){
  if(dirChange[i] == 0) next
  #TS_EX_T[[i]]$paleoData_values = TS_EX_T[[i]]$paleoData_values*dirChange[i]
  TS_EX_T[[i]]$dirEx = dirChange[i]
  
}

tidy.T <- tidyTs(TS_EX_T, age.var = "age")

df.T <- tidy.T %>% 
  filter(between(age,3000,5000)) %>% 
  group_by(eventEX) %>% 
  arrange(eventEX) %>%
  arrange(dirEx)

map.all.T <- mapTs(TS_EX_T, color = 'eventEX', global = T)+
  scale_color_manual(values = c("FALSE" = "#1B9E77", "TRUE" = "#7570B3"))+
  ggtitle("Temperature: Excursion")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/map.all.T_EX.pdf", plot = map.all.T)
dev.off()

all.T  <- plotTimeseriesStack(df.T, color.var = "eventEX", time.var = "age")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.T$dataSetName))+28,alpha = 0.2)+
  ggtitle("Temperature: Excursion")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.T_EX.pdf", plot = all.T, limitsize = F, height = 56)
dev.off()

# European records

sapply(TS_EX_M, "[[", "dirEx")
df.T <- tidy.T %>% 
  filter(between(age,3000,5000)) %>% 
  filter(between(geo_longitude,0,50)) %>%  
  filter(between(geo_latitude,20,60)) %>% 
  group_by(paleoData_TSid) %>% 
  arrange(archiveType) 

EU.T <- plotTimeseriesStack(df.T, color.var = "archiveType", color.ramp = "black", time.var = "age")
ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/EU.T.pdf", plot = EU.T, limitsize = F)

#####################################
# Mean shift

load(file.path( '/projects/pd_lab/sha59/4ka/RData/MS_results_plusNull_complete.RData'))

TS_MS  <- analysis_2b

for (i in 1:length(TS_MS)) {
  # Filter records that don't contain event year in age range
  if (min(TS_MS[[i]]$age) > eventYr || max(TS_MS[[i]]$age) < eventYr) {
    TS_MS[[i]]$useMS = 0
    #print('Out of range')
  }
  
  #Only include annual, winterOnly, and summerOnly (exclude winter+ and summer+)
  if (length(TS_MS[[i]]$interpretation1_seasonalityGeneral) > 0) {
    if (tolower(TS_MS[[i]]$interpretation1_seasonalityGeneral) == 'summer+' | 
          tolower(TS_MS[[i]]$interpretation1_seasonalityGeneral) == 'winter+') {
      TS_MS[[i]]$useMS = -1
      print(TS_MS[[i]]$interpretation1_seasonalityGeneral)
    }
  }
}

for (i in 1:length(TS_MS)) {
  
  TS_MS[[i]]$eventMS = 0
  TS_MS[[i]]$dirMS = 0
  if (!is.na(TS_MS[[i]]$sig_brks) && any(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow)) {
    TS_MS[[i]]$eventMS = 1
    eve_i = which(TS_MS[[i]]$sig_brks >= eventYr - param$eventWindow & TS_MS[[i]]$sig_brks <= eventYr + param$eventWindow)
    TS_MS[[i]]$dirMS = TS_MS[[i]]$brk_dirs[eve_i]
  }
  
}

for(i in 1:length(TS_MS)){
  
  if(TS_MS[[i]]$eventMS == 0){
    TS_MS[[i]]$eventMs = "FALSE"
  }else{
    TS_MS[[i]]$eventMs = "TRUE"
  } 
  
}

interps <- pullTsVariable(TS = TS_MS,variable = "interpretation1_variable")

# Moisture
index <- which(interps == 'M' | interps == 'P' | interps == 'P-E' | interps ==  'P/E')
TS_MS_M <- TS_MS[index]
dirs <- pullTsVariable(TS = TS_MS_M,variable = "interpretation1_interpDirection")
dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
dirs[dirs == 'negative'] = -1
dirs[dirs == 'NA' | is.na(dirs)] = 0
dirs = as.numeric(dirs)
dirEvents = pullTsVariable(TS = TS_MS_M, variable = "dirMS")
dirChange = dirs * dirEvents 
for(i in 1:length(TS_MS_M)){
  
  if(dirChange[i] == 0) next
  TS_MS_M[[i]]$dirMS = dirChange[i]
  
}

tidy.M <- tidyTs(TS_MS_M, age.var = "age")

df.M <- tidy.M %>% 
  filter(between(age,3000,5000)) %>% 
  filter(eventMs == "FALSE") 
  
all.M  <- plotTimeseriesStack(df.M, color.var = "eventMs", time.var = "age",
                              color.ramp = "#1B9E77")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.M$dataSetName))+1,alpha = 0.2)+
  ggtitle("Moisture: Mean Shift (FALSE)")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.M_MS_FALSE.pdf", plot = all.M, limitsize = F, height = 50)

map.all.M <- mapTs(TS_MS_M, color = 'eventMs', global = T, size = 2)+
  scale_color_manual(values = c("FALSE" = "#1B9E77", "TRUE" = "#7570B3"))+
  ggtitle("Moisture: Mean Shift")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/map.all.M_MS.pdf", plot = map.all.M)

df.M <- tidy.M %>% 
  filter(between(age,3000,5000)) %>% 
  filter(eventMs == "TRUE") %>%
  filter(dirMS == 1) 

all.M  <- plotTimeseriesStack(df.M, color.var = "eventMs", time.var = "age", color.ramp = "#7570B3")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.M$dataSetName))+3,alpha = 0.2)+
  ggtitle("Moisture: Mean Shift (TRUE & POSITIVE)")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.M_MS_TRUE_POSITIVE.pdf", plot = all.M, limitsize = F, height = 50)

df.M <- tidy.M %>% 
  filter(between(age,3000,5000)) %>% 
  filter(eventMs == "TRUE") %>%
  filter(dirMS == -1) 

all.M  <- plotTimeseriesStack(df.M, color.var = "eventMs", time.var = "age", color.ramp = "#7570B3")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.M$dataSetName))+4,alpha = 0.2)+
  ggtitle("Moisture: Mean Shift (TRUE & NEGATIVE)")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.M_MS_TRUE_NEGATIVE.pdf", plot = all.M, limitsize = F, height = 50)

# Temperature
index <- which(interps == 'T' | interps == 'TM')
TS_MS_T <- TS_MS[index]
dirs <- pullTsVariable(TS = TS_MS_T,variable = "interpretation1_interpDirection")
dirs[dirs == 'positive' | dirs == 'positve' | dirs == 'postitive'] = 1
dirs[dirs == 'negative'] = -1
dirs[dirs == 'NA' | is.na(dirs)] = 0
dirs = as.numeric(dirs)
dirEvents = pullTsVariable(TS = TS_MS_T, variable = "dirMS")
dirChange = dirs * dirEvents 
for(i in 1:length(TS_MS_T)){
  if(dirChange[i] == 0) next
  
  TS_MS_T[[i]]$dirMS = dirChange[i]
  
}

tidy.T <- tidyTs(TS_MS_T, age.var = "age")

map.all.T <- mapTs(TS_MS_T, color = 'eventMs', global = T, size = 2)+
  scale_color_manual(values = c("FALSE" = "#1B9E77", "TRUE" = "#7570B3"))+
  ggtitle("Temperature: Mean Shift")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/map.all.T_MS.pdf", plot = map.all.T)

df.T <- tidy.T %>% 
  filter(between(age,3000,5000)) %>% 
  filter(eventMs == "FALSE") 

all.T  <- plotTimeseriesStack(df.T, color.var = "eventMs", 
                              time.var = "age", color.ramp = "#1B9E77")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.T$dataSetName))+1,alpha = 0.2)+
  ggtitle("Temperature: Mean Shift (FALSE)")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.T_MS_FALSE.pdf", plot = all.T, limitsize = F, height = 200)

df.T <- tidy.T %>% 
  filter(between(age,3000,5000)) %>% 
  filter(eventMs == "TRUE") %>%
  filter(dirMS == 1) 

all.T  <- plotTimeseriesStack(df.T, color.var = "eventMs", time.var = "age", color.ramp = "#7570B3")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.T$dataSetName))+17,alpha = 0.2)+
  ggtitle("Temperature: Mean Shift (TRUE & POSITIVE)")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.T_MS_TRUE_POSITIVE.pdf", plot = all.T, limitsize = F, height = 50)

df.T <- tidy.T %>% 
  filter(between(age,3000,5000)) %>% 
  filter(eventMs == "TRUE") %>%
  filter(dirMS == -1) 

all.T  <- plotTimeseriesStack(df.T, color.var = "eventMs", time.var = "age", color.ramp = "#7570B3")+
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 3900, xmax = 4500, ymin = 0, ymax = length(unique(df.T$dataSetName))+17,alpha = 0.2)+
  ggtitle("Temperature: Mean Shift (TRUE & NEGATIVE)")

ggsave(filename = "/projects/pd_lab/sha59/4ka/4.2k_event/all.T_MS_TRUE_NEGATIVE.pdf", plot = all.T, limitsize = F, height = 50)

# Test for well known 4.2 ka records

which(sapply(data_EX, "[[", "dataSetName") == "Corchia.Regattieri.2014")

Corchia  <- data_EX[[46]]
plot(x = Corchia$age, y = Corchia$paleoData_values, type = "l")
Corchia$interpretation1_interpDirection
