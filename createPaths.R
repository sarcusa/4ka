createPaths = function() {
  
  mainDir = '/Users/hannah/Dropbox/2019'
  
  dir.create(mainDir)
  dir.create(file.path(mainDir, 'excursion'))
  dir.create(file.path(mainDir, 'mean_shift'))
  dir.create(file.path(mainDir, 'broken_stick'))
  dir.create(file.path(mainDir, 'histograms'))
  dir.create(file.path(mainDir, 'RData'))
  
  return(mainDir)
}