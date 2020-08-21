# Set Parameters 

# Parameters common to all analyses

CName = "HoloceneAbruptChange" #the database name
CVers = "0_9_0" #the database version
OutDat = 'TS_climateInterp_2020.RData' #What should the processed dataset be called?
eventYrs = seq(1000,11000,by = 400) #seq(1000,2000, by = 1000) # was seq(1000,11000,by = 400)
numIt = 250 # number of iterations in null model, was 1000
res = 5 # grid resolution, deg
radius = 2000 # search radius, km for the spatial excursion
climateVar = 'T'   # M for moisture, T for temperature, All for all
plotOpt = T # Plotting?
ncores = 128 #how many cores to run in parallel?

# Excursion

event_yr = 3800 # Set the event year (select any year, this will change)
event_window = 600 # (event_yr +/- event_window/2) as the event window
ref_window = 500  # Set the duration of the reference windows 
resCriteria = 50 # Resolution criteria (yr) required to include record
sigNum = 2 # Set how many sigmas below and above the mean define an excursion
mainDir = createPaths() # do not change

# Mean shift & Trend Change

maxDiff = 1000    # min segment length between change points
alpha = 0.05      # set confidence level for mean difference testing
gaussOpt = F      # option to gaussianize values before analysis
eventYr = 8200     # event year, yr BP
eventWindow = 199.9  # total event window = eventYr +/- eventWindow
eventDetector = 'MS'   # MS for mean shift, BS for broken stick

