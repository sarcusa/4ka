# Hannah Kolus, hrk37@nau.edu
# 
# Determines whether an excursion event has occurred within the specified event window.
# Excursion events are defined as two consecutive values within the event window that
# are more extreme than the avg +/- X std of the reference windows.

## --- SET PARAMETERS DEFINING THE ANALYSIS WINDOW AND EVENTS --- ##

event_yr = 4200     # Set the event year
event_window = 600  # (event_yr +/- event_window/2) as the event window
ref_window = 500    # Set the duration of the reference windows 

resCriteria = 50    # Resolution criteria (yr) required to include record
sigNum = 2          # Set how many sigmas below and above the mean define an excursion

# yr_start:yr_end defines boundaries of analysis (i.e. both reference windows and the event window)
yr_start = event_yr - event_window / 2 - ref_window
yr_end = event_yr + event_window / 2 + ref_window

# event_start:event_end defines the boundaries of the event
event_start = event_yr - event_window / 2
event_end = event_yr + event_window / 2

## --- FILTER RECORDS --- ##

# UPDATE THIS SECTION AFTER FIGURING OUT HOW TO INTEGRATE LiPD FILES - THIS MAY GO INSIDE THE FOR LOOP

# 1. Check that there is a climate interpretation

# 2. Check that there are data points in event window and both reference windows

# 3. Calculate median resolution, check that this meets resolution criteria

# 4. Save all data and fields: 
#         loc, site name, proxy type, archive type, resolution, units, variable name, ages, values

# FOR NOW ASSUME:
# ages = list of length(records that match crtieria) containing ages
# values = list of length(records that match crtieria) containing paleo data values

status = rep(0, length(ages)) # store whether the record is used (1) or not (0) in analysis 
for (i in 1:length(ages)) {
	
	age = ages[[i]]
	vals = values[[i]]
	
	# MIGHT DO ALL THIS PROCESSING IN DIFFERENT FILE, SAME STEP AS EXTRACTING FROM ALL LiPD FILES
	# Setting all missing values to NA
  age[age == "NaN"] = NA
  vals[vals == "NaN"] = NA
  age[!is.finite(age)] = NA
  vals[!is.finite(vals)] = NA
  vals[vals == -9999] = NA
  vals[vals == -999] = NA
  vals[vals == 999999] = NA
  vals[vals == 99999] = NA
  vals[vals == 9999] = NA
  
  # ---------------------- MORE DATA FILTERS ----------------------- #
  if (all(is.na(vals))) {
    print(paste('NA VECTOR: INDEX ', i, 'SITE: ', siteNames[i]))
    status[i] = -1
    next()
  }
  
  if (length(age) != length(vals)) {
    print(paste('DIFFERENT VECTOR LENGTHS: INDEX ', i, 'SITE: ', siteNames[i]))
    status[i] = -1
    next()
  }
  
  if (sum(diff(vals) == 0, na.rm = T) > length(vals) / 2) {
    print(paste('STEP FUNCTION/QUANTIZED VALUES: INDEX ', i))
    status[i] = -1
    next()
  }
  
  if (all.equal(diff(age), rep(50, length(age) - 1)) == TRUE) {
    paste('50-YR INTERPOLATED DATA SET: INDEX ', i)
    status[i] = -1
    next()
  }
  
  if (!all(sort(age[!is.na(age)]) == age[!is.na(age)])) {
    print(paste('OUT OF ORDER AGES: INDEX ', i, 'SITE: ', siteNames[i]))
    ind_order = order(age)
    age = age[ind_order]
    vals = vals[ind_order]
  }
	
	# get rid of NAs
  inds = which(!is.na(vals))
  age = age[inds]
  vals = vals[inds]
  inds = which(!is.na(age))
  age = age[inds]
  vals = vals[inds]

	
}	
