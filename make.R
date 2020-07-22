dir = '/home/sha59/4ka'
setwd(dir)

#envir <- new.env(parent = globalenv())

source('packages.R')
source('functions.R')
source('set_parameters.R')
source('plan_prep.R')
source('plan_excursion.R')
source('plan_meanshift.R')
source('plan_trendchanges.R')
source('plan_plotting.R')
source('plan_var.R') #, local = envir

file.exists("report.Rmd")
#drake_hpc_template_file("slurm_batchtools.tmpl")
#future::plan(batchtools_slurm, template = "slurm_batchtools.tmpl")

# Download necessary data
set.seed(1)
external_data  = readLipd(file_in("http://lipdverse.org/HoloceneAbruptChange/0_9_0/HoloceneAbruptChange0_9_0.zip"))
print("data downloaded")
data_sub = list.sample(external_data, size = 50, replace = F)
print("subset done")
data_all  =  extractTs(data_sub)
print("TS extracted")

#vis_drake_graph(my_plan)

#make(
#  envir$my_plan,
#  parallelism = "future",
#  jobs = 4,
#  console_log_file = "drake.log",
#  lock_envir = FALSE,
#  envir = envir
#)

#make(plan, max_expand = 2)
#vis_drake_graph(plan)
# plot(plan)

make(my_plan, lock_envir = FALSE, lock_cache = FALSE)
