createPaths = function() {
  
  if(identical(dir,output_destination)){
    mainDir = dir
  }else{
    mainDir = output_destination
  }
  
  dir.create(mainDir)
  dir.create(file.path(mainDir, 'excursion'))
  dir.create(file.path(mainDir, 'mean_shift'))
  dir.create(file.path(mainDir, 'broken_stick'))
  dir.create(file.path(mainDir, 'histograms'))
  dir.create(file.path(mainDir, 'RData'))
  
  return(mainDir)
}
