##################
# input parameters for binning and plotting AAVSO data
llightcurve_name <- "data/aavso_latest_data.csv"
# the name of the asassn file
asassn.csv.file <- "data/ASASSN_latest.csv"
merge.asassn <- TRUE # merge asasn data into light curve
asassn.code <- "ASASSN"

maxairmass <- 2.8 # air mass values above this will be filtered out, as well as missing air masses. Set >= 100 to turn this off
maxuncertainty <- 0.2  # maximum AAVSO uncertainty estimate
maxBinUncertainty <- 0.2 # worst standard deviation to accept for a binned set of observations
wildsd <- 10.0 # worst number of standard deviations from mean allowed

earliestJD = 2457290 # only data on or after this JD will be used
#earliestJD = 2457800 # only data on or after this JD will be used
#earliestJD <- 2457700
startPlot <- earliestJD
#startPlot <- 2457690
#startPlot <- 2457980
#startPlot <- 2457880
#startPlot <- 2457580
startPlot <- 2457990
stop.plot <- NA
#stop.plot <- 2458020
add.predict <- 0 # additional prediction days to add to plot
plotRelTimes <- TRUE
##########
includeExclude <- TRUE # TRUE if your list of observer codes is to to be included, FALSE if excluded or not used
ExclCodes <- "None"
ExclCodes <- c("LDJ","DUBF","PXR","DKS","OJJ","HBB","SDB","VBPA","OAS","MJB","MATA","JSJA",
				"WROC","MAND","VBPA","NOT","PALE","GKA","AMID","SGEA","ELYA","GCJ","LBG","HJW","OAR","ASASSN") # V & B
#ExclCodes <- c("LDJ")
#ExclCodes <- c("ASASSN")
#ExclCodes <- c("DUBF")
#ExclCodes <- c("SWIA","CPP","DRZ","LPAC","HDHA") 	# B band exclusions
#ExclCodes <- c("DUBF","GKA","BPAD","LPB","SJAR","LBG","LDJ","LWHA") # R ensemble
#ExclCodes <- c("OAR","OJJ","GKA","MJB","SJAR","LWHA","LBG","LPB","LDJ","CMP","JM") # I ensemble
ExclCodes <- c("DUBF","MJB","LDJ","GKA","ELYA","HJW","JSJA","VBPA","DKS","OAR","JM","HBB") # B ensemble
plotMee <- NA # do not highlight any particular observer code
#plotMee <- "ASASSN"
meeColor <- "darkviolet"
weightless <- NA
#weightless <- c("JM","LDJ") # observers to plot, but not use in fit.
#weightless <- c("ASASSN") # observers to plot, but not use in fit.
########
allBands <- data.frame(bandinQ=c("I","R","V","B"),plotColor=c("darkviolet","red","green","blue"), stringsAsFactors=FALSE)
allBands <- data.frame(bandinQ=c("V"),plotColor=c("darkgreen"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V","B"),plotColor=c("green","blue"), stringsAsFactors=FALSE)
allBands <- data.frame(bandinQ=c("B"),plotColor=c("blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I"),plotColor=c("darkviolet"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("R"),plotColor=c("red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R"),plotColor=c("darkviolet","red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R","B"),plotColor=c("darkviolet","red","blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("R","V","B"),plotColor=c("red","green","blue"), stringsAsFactors=FALSE)

########################
deltaJD <- 1.0 # bin width in days
########################
plotExcluded <- FALSE # set to TRUE to plot the points in the lightcurve not used in the fit.
plotQuadratic <- FALSE # set to TRUE to plot a quadratic fit
generateTS <- FALSE # set to TRUE to creat a time series from the data
tsBinWidth <- 5.0 # time series bin width in days. Important if generateTS is TRUE
smoothTS <- TRUE # set to TRUE to smooth the times series
tsSmoothOrder <- 6 # the order for the moving average to smooth the time series
tryLQS <- FALSE # set to TRUE is you want to try resistant regression.
userlm <- FALSE # set to TRUE to use robust lm, or rlm() if not using MARS.
plotResiduals <- TRUE # set to true to plot the residuals vs. time
res.plot.type = "p" # should be a legit plot() type
plot2Lines <-  FALSE  # two line feature doesn't work well
lqsColor <- "darkgreen"
ebar.color <- "grey"	 # color to use for error bars


####### MARS
marsOrder <- 11 # maximum number of knots
marsPenalty <- 3 # set to 0 to avoid penalizing knots in pruning pass
#marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "exhaustive" # set to "none" to avoid pruning
splineRaw <-  FALSE # do the spline on the raw lightcurve, not binned.

############################## Comparison Stars
#okComparison <- "(000-?BLS-?556)|(000-?BLS-?551)|(000-?BLS-?553)|(000-?BLS-?552)|(000-?BLS-?554)|(000-?BLS-?549)|(BLS-549)|(BLS-555)|(000-?BLS-?555)|(108)|(113)|(116)|(118)|(121)|(124)|(128)|(ENSEMBLE)|(APASS20062365[+-]442738)|(3162-1853-1)" # regular expression from AAVSO photometry table
okComparison <- "."
bad.comp.star <- "(000-?BLS-549)|(000-?BLS-?549)|(BLS-549)"
#bad.comp.star <- NA

##########################################################################
pretty.interval <- 100 # controls rounding of plot limits
pretty.JD.interval <- 100
