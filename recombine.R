### Recombining the datasets into an ordered list

library(gtools)


files  <- list.files(path = "/projects/pd_lab/sha59/4ka/RData/", pattern = "BS_results_plusNull_")
print(files)
files  <-  mixedsort(sort(files))
print(files)

wanted  <-  files[1:26]

dir = "/projects/pd_lab/sha59/4ka/RData/"
load("/projects/pd_lab/sha59/4ka/RData/BS_results_complete.RData")
TS_BS_ori = TS_BS

ind  <- list()

for(i in 1:length(wanted)){
  
  load(paste0(dir,wanted[i]))
  
  #finds the index of filled lists within each record's null break directions
  filled_dirs <- sapply(TS_BS, "[[", "null_brk_dirs")
  filled_pts  <- sapply(TS_BS, "[[", "null_brk_pts")
  filled_ptsErr <-  sapply(TS_BS, "[[", "null_brk_ptsErr")
  ind[[i]]  <- which(!unlist(lapply(filled_dirs, is.null))) 
 
  file  <- ind[[i]]
  
  # Finds the appropriate record and adds the new values to the original
  for(j in 1:length(file)){
  
    TS_BS_ori[[file[[j]]]]$null_brk_pts = filled_pts[[file[[j]]]]
    TS_BS_ori[[file[[j]]]]$null_brk_ptsErr = filled_ptsErr[[file[[j]]]]
    TS_BS_ori[[file[[j]]]]$null_brk_dirs = filled_dirs[[file[[j]]]]
  
  }
  
  
} # end of loop

# This is a manual check for duplicates
print(ind)  
which(duplicated(ind))
wanted
wanted[which(duplicated(ind))]
  
TS_BS = TS_BS_ori
save(TS_BS, file = "/projects/pd_lab/sha59/4ka/RData/BS_results_plusNull_complete.RData")

# Plot of null break point detection

load("/projects/pd_lab/sha59/4ka/RData/BS_results_plusNull_complete.RData")

for(i in 1:length(TS_BS)){
  data[i] =  length(which(!is.na(TS_BS[[i]]$null_brk_pts)))
}

data = (data/250)*100

barplot(data, xlab = seq(1, length(data), 1))
hist(data)
