dir = "/projects/pd_lab/sha59/4ka"
setwd(dir)

#envir <- new.env(parent = globalenv())

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
source('plan_var.R') #, local = envir

#file.exists("report.Rmd")

#future::plan(future::multisession, workers = 4)
future::plan(batchtools_slurm, template = file.path(dir,"slurm_batchtools.tmpl"))

# Download necessary data
#set.seed(1)
external_data  = readLipd(file_in("http://lipdverse.org/HoloceneAbruptChange/0_9_0/HoloceneAbruptChange0_9_0.zip"))
print("data downloaded")
#data_sub = list.sample(external_data, size = 50, replace = F)
#print("subset done")
#data_all  =  extractTs(data_sub)
#print("TS extracted")
data_all = extractTs(external_data)

#vis_drake_graph(my_plan, file = "dependency_graph.html")
#vis_drake_graph(my_plan)
#predict_runtime(my_plan)

make(my_plan, lock_envir = FALSE, lock_cache = FALSE, 
     parallelism = "future", jobs = 8)
