BrokenStick_null_test <- function(data_in, param){
  
  print("Running brokenStick_null now")
  load(paste0(mainDir,"RData/BS_results_complete.RData"))
      
  chunks  <- length(TS_BS)
  chunk_size  <- c(1,floor(chunks/4), ceiling(chunks/4), ceiling(chunks/4)*2,
                   ceiling(chunks/4)*2+1, ceiling(chunks/4)*3, 
                   ceiling(chunks/4)*3+1, chunks)
  
  future::plan(batchtools_slurm, template = file.path(folder.dir,"slurm_batchtools.tmpl"), resources = list(cores = 16))
  
  a <- future({
    A = BSWrapper(TS = TS_BS,par = param,chunkS = chunk_size[1],
              chunkE = chunk_size[2])
    fileName = file.path(mainDir, 'RData', 'BS_results_plusNull_A.RData')
    save(A, file = fileName)
    
  })
  b <- future({
    B = BSWrapper(TS = TS_BS,par = param,chunkS = chunk_size[3],
               chunkE = chunk_size[4])
    fileName = file.path(mainDir, 'RData', 'BS_results_plusNull_B.RData')
    save(B, file = fileName)
  })
  c <- future({
    C = BSWrapper(TS = TS_BS,par = param,chunkS = chunk_size[5],
               chunkE = chunk_size[6])
    fileName = file.path(mainDir, 'RData', 'BS_results_plusNull_C.RData')
    save(C, file = fileName)
  })
  d <- future({
    D = BSWrapper(TS = TS_BS,par = param,chunkS = chunk_size[7],
               chunkE = chunk_size[8])
    fileName = file.path(mainDir, 'RData', 'BS_results_plusNull_D.RData')
    save(D, file = fileName)
  })
  
  
  chunk_out <- value(c(a[[1]],b[[1]],c[[1]],d[[1]]))
  
  records_error <- value(c(a[[2]],b[[2]],c[[2]],d[[2]]) )
  
  plan(sequential)
  save(chunk_out, file = 'BS_results_plusNull_test.RData')
  write.table(records_error, file = file.path(dir,"skipped_records_BS.txt"))
  
  output  <- list(chunk_out,s)
  return(output)
  
}
