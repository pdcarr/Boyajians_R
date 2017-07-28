##################
# input parameters for binning and plotting AAVSO data
llightcurve_name <- "data/aavso_latest_data.csv"
maxairmass <-1.25 # air mass values above this will be filtered out, as well as missing air masses. Set >= 100 to turn this off
maxuncertainty <- 0.1  # maximum AAVSO uncertainty estimate
maxBinUncertainty <- 0.1 # worst standard deviation to accept for a binned set of observations
wildsd <- 10.0 # worst number of standard deviations from mean allowed

earliestJD = 2457294 # only data on or after this JD will be used
#earliestJD <- 2457700
startPlot <- earliestJD
#startPlot <- 2457880
#startPlot <- 2457580
plotRelTimes <- TRUE
##########
includeExclude <- TRUE # TRUE if your list of observer codes is to to be included, FALSE if excluded or not used
ExclCodes <- "None"
#ExclCodes <- c("JM","LDJ","ELYA","DKS","OJJ","OAR","ATE","BPAD","HJW")
#ExclCodes <- c("JM","LDJ","DUBF","HJW","PXR","DKS","OJJ","HBB","SDB","VBPA","OAS","MJB","MATA","JSJA","UJHA","WROC","MAND","HDHA","ELYA","VBPA","NOT","PALE") # V & B
#ExclCodes <- c("LDJ","DUBF","HJW","PXR","DKS","OJJ","HBB","SDB","VBPA","OAS","MJB","MATA","JSJA","WROC","MAND","HDHA","ELYA","VBPA","NOT","PALE") # V & B
ExclCodes <- c("LDJ","DUBF","HJW","PXR","DKS","OJJ","HBB","SDB","VBPA","OAS","MJB","MATA","JSJA","WROC","MAND","HDHA","ELYA","VBPA","NOT","PALE","KTHC") # V
#ExclCodes <- c("JM","LDJ","DUBF","HJW","PXR","DKS","OJJ","HBB","SDB","VBPA","OAS","MJB","MATA","JSJA","UJHA","WROC","MAND","HDHA","ELYA","VBPA","NOT") # V & B
#ExclCodes <- c("DUBF","DKS","ELYA","OAR","ATE","HJW","BPAD","OJJ","LBG","LDJ","UJHA","OYE","GFRB","OAS","MJB","EEY") # V ensemble
#ExclCodes <- c("DUBF","GKA","BPAD","LPB","SJAR","LBG","LDJ","LWHA") # R ensemble
#ExclCodes <- c("OAR","OJJ","GKA","MJB","SJAR","LWHA","LBG","LPB","LDJ","CMP") # I ensemble
#ExclCodes <- c("DUBF","MJB","LDJ","GKA","ELYA","HJW","DUBF","JSJA","VBPA","DKS") # candidate B ensemble
#ExclCodes <- c("LDJ","DUBF","HBB")
#ExclCodes <- c("SWIA","CPP")	# B band exclusions
#ExclCodes <- c("BJFB","COO","LBG","BMAK") # I band exclusion
#ExclCodes <- c("LDJ","DUBF","UJHA","DKS","HBB")
#ExclCodes <- "JM"
########
plotMee <- NA # do not highlight any particular observer code
#plotMee <- "ELYA"
meeColor <- "darkviolet"
########
allBands <- data.frame(bandinQ=c("I","R","V","B"),plotColor=c("darkviolet","red","green","blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V"),plotColor=c("darkgreen"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V","B"),plotColor=c("green","blue"), stringsAsFactors=FALSE)
allBands <- data.frame(bandinQ=c("B"),plotColor=c("blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I"),plotColor=c("darkviolet"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("R"),plotColor=c("red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R"),plotColor=c("darkviolet","red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R","B"),plotColor=c("darkviolet","red","blue"), stringsAsFactors=FALSE)

########################
deltaJD <- 1.0 # bin width in days
########################
plotExcluded <- FALSE # set to TRUE to plot the points in the lightcurve not used in the fit.
plotQuadratic <- FALSE # set to TRUE to plot a quadratic fit
generateTS <- TRUE # set to TRUE to creat a time series from the data
tsBinWidth <- 8.0 # time series bin width in days. Important if generateTS is TRUE
smoothTS <- TRUE # set to TRUE to smooth the times series
tsSmoothOrder <- 10 # the order for the moving average to smooth the time series
tryLQS <- FALSE # set to TRUE is you want to try resistant regression.
userlm <- TRUE # set to TRUE to use robust lm, or rlm() if not using MARS.
plotMARS <-  TRUE # set to TRUE to try a MARS fit instead of lm() or rlm()
plotResiduals <- TRUE # set to true to plot the residuals vs. time
plot2Lines <-  FALSE  # two line feature doesn't work well
lqsColor <- "darkgreen"
weightedBins <- FALSE # set to TRUE to weight lower uncertainty bins more.

####### MARS
marsOrder <- 11 # maximum number of knots
marsPenalty <- 4 # set to 0 to avoid penalizing knots in pruning pass
marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "backward" # set to "none" to avoid pruning
splineRaw <-  FALSE # do the spline on the raw lightcurve, not binned.
############################## Comparison Stars
okComparison <- "(000-?BLS-?556)|(000-?BLS-?551)|(000-?BLS-?553)|(000-?BLS-?552)|(000-?BLS-?554)|(000-?BLS-?549)|(BLS-549)|(BLS-555)|(000-?BLS-?555)|(108)|(113)|(116)|(118)|(121)|(124)|(128)|(ENSEMBLE)|(APASS20062365[+-]442738)|(3162-1853-1)" # regular expression from AAVSO photometry table

##########################################################################
