# Set Parameters 

#output_destination = dir
output_destination = "/projects/pd_lab/sha59/4ka/Signif_test_MS_99_750/" #Use dir for the same folder as the inputs. 

# Parameters common to all analyses

CName = "HoloceneAbruptChange" #the database name
CVers = "0_10_0" #the database version
OutDat = 'TS_climateInterp_2020.RData' #What should the processed dataset be called?
eventYrs = seq(1000,11000,by = 400) #
numIt = 1000 # number of iterations in null model, was 1000
res = 5 # grid resolution, deg
radius = 2000 # search radius, km for the spatial excursion
climateVar = 'T' # Deprecated. M for moisture, T for temperature, All for all. 
plotOpt = T # Plotting?
ncores = 8 #was 128 #how many cores to run in parallel? Deprecated
det = FALSE #detrend records before running the analyses? Deprecated

# Excursion

event_yr = 3800 # Set the event year (select any year, this will change)
event_window = 500 #(event_yr +/- event_window/2) as the event window. Was 600.
ref_window = 600  # Set the duration of the reference windows. Was 500. 
resCriteria = 50 # Resolution criteria (yr) required to include record
sigNum = 2 # Set how many sigmas below and above the mean define an excursion
mainDir = createPaths() # do not change.

# Mean shift & Trend Change

maxDiff = 750    # min segment length between change points. Was 1000. 
alpha = 0.05      # set confidence level for mean difference testing
gaussOpt = F      # option to gaussianize values before analysis
eventYr = 8200     # event year, yr BP. Placeholder.
eventWindow = 99.9  # total event window = eventYr +/- eventWindow. Was 199.9
eventDetector = 'MS'   # Deprecated. MS for mean shift, BS for broken stick

