# Reference table

dir = "/projects/pd_lab/sha59/4ka" #set the directory of the input scripts
setwd(dir)

source('packages.R')
source('functions.R')
source('set_parameters.R')
source('plan_prep.R')

parameters =list(eventYrs = eventYrs, event_window = event_window, 
                         ref_window = ref_window, resCriteria = resCriteria, 
                         sigNum = sigNum, plotOpt = plotOpt, 
                         mainDir = mainDir, numIt = numIt, 
                         res = res, radius = radius, 
                         eventWindow = eventWindow, 
                         CName = CName, CVers = CVers, 
                         eventDetector = eventDetector, OutDat = OutDat,
                         ncores = ncores, maxDiff = maxDiff)
param = parameters

external_data  = readLipd(file_in("http://lipdverse.org/HoloceneAbruptChange/0_10_0/HoloceneAbruptChange0_10_0.zip"))
data_all = extractTs(external_data)
data = plan_prep(param = parameters)
Name <- unlist(sapply(data, "[[", "paleoData_TSid"))
Values <- lipdR::pullTsVariable(data, "paleoData_values")
Age <- lipdR::pullTsVariable(data, "age")
Site <- unique(lipdR::pullTsVariable(data, "geo_siteName"))

header.names <- c("TVerseID", "Dataset", "Site", "Lat", "Lon", "Event","Year", "Var", "Archive", "Proxy", "OriginalURL", "DOI")

EX_M <- read.csv(file = "/projects/pd_lab/sha59/4ka/RData/event_each_site_EX_M.csv", stringsAsFactors = FALSE, header = F)
colnames(EX_M) <- header.names
EX_M <- subset(EX_M, select = -c(Event,Var,Year))
EX_M <- dplyr::distinct(EX_M)

EX_T <- read.csv(file = "/projects/pd_lab/sha59/4ka/RData/event_each_site_EX_T.csv", stringsAsFactors = FALSE, header = F)
colnames(EX_T) <- header.names
EX_T <- subset(EX_T, select = -c(Event,Var,Year))
EX_T <- dplyr::distinct(EX_T)

EX <- dplyr::union(EX_M,EX_T)

indx <- which(Name %in% EX$TVerseID)
EX_proxydata <- list(SiteName = Name[indx], Values = Values[indx], Age = Age[indx])

index <- vector()
for(i in 1:nrow(EX)){
  ind <- which(EX$TVerseID[i] == Name)
  index[i] <- ind
}

EX$OriginalURL <- Name[index]
EX <- subset(EX, select = -c(Dataset))

# Mean shift
MS_M <- read.csv(file = "/projects/pd_lab/sha59/4ka/RData/event_each_site_MS_M.csv", stringsAsFactors = FALSE, header = F)
colnames(MS_M) <- header.names
MS_M <- subset(MS_M, select = -c(Event,Var,Year))
MS_M <- dplyr::distinct(MS_M)

MS_T <- read.csv(file = "/projects/pd_lab/sha59/4ka/RData/event_each_site_MS_T.csv", stringsAsFactors = FALSE, header = F )
colnames(MS_T) <- header.names
MS_T <- subset(MS_T, select = -c(Event,Var,Year))
MS_T <- dplyr::distinct(MS_T)

MS <- dplyr::union(MS_M,MS_T)

indx <- which(Name %in% MS$TVerseID)
MS_proxydata <- list(SiteName = Name[indx],Values = Values[indx], Age = Age[indx])

index <- vector()
for(i in 1:nrow(MS)){
  ind <- which(MS$TVerseID[i] == Name)
  index[i] <- ind
}
MS$OriginalURL <- Name[index]
MS <- subset(MS, select = -c(Dataset))

# Joining tables

All <- dplyr::full_join(EX,MS)
All$Excursion  <- rep(NA,nrow(All))
All$MeanShift <- rep(NA,nrow(All))
All$TrendChange <- rep(NA,nrow(All))

for(i in 1:nrow(All)){
  if(All$Site[i] %in% EX$Site){
    All$Excursion[i]  <- "X"
  }
  if(All$Site[i] %in% MS$Site){
    All$MeanShift[i]  <- "X"
    All$TrendChange[i]  <- "X"
  }
}
write.table(x = All, file = paste0(dir,"/RData/references.csv"))

unique_all <- unique(All$Site)

# How many records of temperature?

T_all <- dplyr::distinct(dplyr::union(MS_T, EX_T))

# Types of archives and proxy

archives <- T_all$Archive
proxies  <- T_all$Proxy

prox = c(length(proxies[which(proxies == 'other microfossil')]),
             length(proxies[which(proxies == 'alkenone')]),
             length(proxies[which(proxies == 'Mg/Ca')]),
             length(proxies[which(proxies == 'pollen')]),
             length(proxies[which(proxies == 'chironomid')]),
             length(proxies[which(proxies == 'isotope')]),
             length(proxies[which(proxies == 'other ice')]),
             length(proxies[which(proxies == 'other biomarker')]),
             length(proxies[which(proxies == 'biophysical')]),
             length(proxies),length(unique(T_all$Site)))
prox_df_T = data.frame(Proxy = c(unique(proxies),"records", "site"), Number = prox)
prox_df_T = dplyr::mutate(prox_df_T, Percentage = round((Number/Number[which(Proxy == "records")])*100,digits = 0))
prox_df_T$Percentage[which(prox_df_T$Proxy == "site")]  <- NA

arch = c(length(archives[which(archives == 'MarineSediment')]),
           length(archives[which(archives== 'LakeSediment')]),
           length(archives[which(archives == 'GlacierIce')]),
           length(archives[which(archives == 'Peat')]),
           length(archives[which(archives == 'Speleothem')]),
           length(archives[which(archives == 'Ice-other')]),
           length(archives[which(archives == 'Wood')]),
           length(archives[which(archives== 'Midden')]),NA,
           length(archives),length(unique(T_all$Site)))
arch_df_T = data.frame(Archive = c(unique(archives),NA,"records", "site"), Number = arch)
arch_df_T = dplyr::mutate(arch_df_T, Percentage = round((Number/Number[which(Archive == "records")])*100,digits = 0))
arch_df_T$Percentage[which(arch_df_T$Archive == "site")]  <- NA

T_ProxyArchive <- cbind(arch_df_T,prox_df_T)

write.table(T_ProxyArchive,file = paste0(dir,"/RData/T_numbers.csv"),sep = ",")

#How many records of moisture?

M_all <- dplyr::distinct(dplyr::union(MS_M, EX_M))

archives <- M_all$Archive
proxies  <- M_all$Proxy

prox = c(length(proxies[which(proxies == 'pollen')]),
         length(proxies[which(proxies == 'isotope')]),
         length(proxies[which(proxies == 'biophysical')]),
         length(proxies[which(proxies == 'other microfossil')]),
         length(proxies[which(proxies == 'chironomid')]),
         length(proxies[which(proxies == 'hybrid')]),
         length(proxies[which(proxies == 'other ice')]),
         length(proxies),length(unique(M_all$Site)))
prox_df_M = data.frame(Proxy = c(unique(proxies),"records", "site"), Number = prox)
prox_df_M = dplyr::mutate(prox_df_M, Percentage = round((Number/Number[which(Proxy == "records")])*100,digits = 0))
prox_df_M$Percentage[which(prox_df_M$Proxy == "site")]  <- NA

arch = c(length(archives[which(archives == 'LakeSediment')]),
         length(archives[which(archives== 'MarineSediment')]),
         length(archives[which(archives == 'Peat')]),
         length(archives[which(archives == 'Speleothem')]),
         length(archives[which(archives == 'Midden')]),
         length(archives[which(archives == 'Wood')]),
         length(archives[which(archives == 'GlacierIce')]),
         length(archives),length(unique(M_all$Site)))
arch_df_M = data.frame(Archive = c(unique(archives),"records", "site"), Number = arch)
arch_df_M = dplyr::mutate(arch_df_M, Percentage = round((Number/Number[which(Archive == "records")])*100,digits = 0))
arch_df_M$Percentage[which(arch_df_M$Archive == "site")]  <- NA

M_ProxyArchive <- cbind(arch_df_M,prox_df_M)

write.table(M_ProxyArchive,file = paste0(dir,"/RData/M_numbers.csv"),sep = ",")

# All records counted in plotProxyMap script

load("/projects/pd_lab/sha59/4ka/RData/final_TS_T.RData")
TS_T <- TS_all
load("/projects/pd_lab/sha59/4ka/RData/final_TS_M.RData")
TS_M <- TS_all

proxiesT <- pullTsVariable(TS_T, variable = "paleoData_proxyGeneral")
proxiesM <- pullTsVariable(TS_M, variable = "paleoData_proxyGeneral")

archivesT  <- pullTsVariable(TS_T, variable = "archiveType")
archivesM  <- pullTsVariable(TS_M, variable = "archiveType")

IDT <- pullTsVariable(TS_T, "paleoData_TSid")
IDM <- pullTsVariable(TS_M, "paleoData_TSid")

IDall <- c(IDT,IDM)

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
