##################
# input parameters for binning and plotting AAVSO data
llightcurve_name <- "aavso_latest_data.csv"
maxairmass <- 1.9 # air mass values above this will be filtered out, as well as missing air masses. Set >= 100 to turn this off
maxuncertainty <- 0.029  # maximum AAVSO uncertainty estimate
maxBinUncertainty <- 0.1 # worst standard deviation to accept for a binned set of observations
wildsd <- 10.0 # worst number of standard deviations from mean allowed

earliestJD = 2457294 # only data on or after this JD will be used
#earliestJD <- 2457700
startPlot <- earliestJD
startPlot <- 2457880
plotRelTimes <- TRUE
##########
includeExclude <- FALSE # TRUE if your list of observer codes is to to be included, FALSE if excluded or not used
ExclCodes <- "None"
#ExclCodes <- c("LDJ","DUBF","ELYA","HJW","JM")
#ExclCodes <- c("DUBF","LDJ","OAR","SGEA","DKS","OJJ","LPB","BPAD")
#ExclCodes <- c("DUBF","LDJ","OAR","SGEA","DKS")
#ExclCodes <- c("DUBF","ELYA","LPB","OJJ","HJW","OAR")
#ExclCodes <- c("ATE","OJJ") # observers to be used/not used in the fit. set to an invalid code (.e.g "None") if not interested
#ExclCodes <- c("ATE") # observers to be used/not used in the fit. set to an invalid code (.e.g "None") if not interested
#ExclCodes <- c("LDJ","UJHA","DKS","OJJ","JM","DUBF","ELYA","HJW")
#ExclCodes <- c("JM","UJHA","LDJ","OAR","LPB","DUBF","ELYA","DKS","OJJ","BSM","SDB","SWIA","VBPA","OAS","MJB","PXR")	 # B and V ensemble
#ExclCodes <- c("JM","LDJ","OAR","LPB","DUBF","ELYA","DKS","OJJ","BSM","SDB","SWIA","VBPA","OAS","MJB")	 # B and V ensemble
#ExclCodes <- c("DUBF","OAR","LDJ","LPB","ELYA","LBG","OJJ","JM","SGEA")
#ExclCodes <- c("DUBF","DKS","ELYA","OAR","ATE","HJW","BPAD","OJJ","LBG","LDJ","UJHA","OYE","GFRB","OAS","MJB","EEY") # V ensemble
#ExclCodes <- c("DUBF","GKA","BPAD","LPB","SJAR","LBG","LDJ","LWHA") # R ensemble
#ExclCodes <- c("OAR","OJJ","GKA","MJB","SJAR","LWHA","LBG","LPB","LDJ","CMP") # I ensemble
#ExclCodes <- c("ELYA","DUBF","OAR","HJW","DKS","LBG","UJHA","JM")
#ExclCodes <- "DUBF"
ExclCodes <- c("BJFB","COO")
########
plotMee <- NA # do not highlight any particular observer code
#plotMee <- "LPAC" # observer code to plot with special character
#plotMee <- "PXR"
#plotMee <- "HJW"
#plotMee <- "DUBF"
#plotMee <- "ELYA"
#plotMee <- "CPP"
#plotMee <- "OAR"
#plotMee <- "MJB"
#plotMee <- "WROC"
#plotMee <- "GKA"
#plotMee <- "JM"
#meeColor <- "darkviolet"
########
allBands <- data.frame(bandinQ=c("I","R","V","B"),plotColor=c("darkviolet","red","green","blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V"),plotColor=c("darkgreen"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V","B"),plotColor=c("green","blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("B"),plotColor=c("blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I"),plotColor=c("darkviolet"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("R"),plotColor=c("red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R","B"),plotColor=c("darkviolet","red","blue"), stringsAsFactors=FALSE)

########################
deltaJD <- 1.0 # bin width in days
########################
plotExcluded <- FALSE # set to TRUE to plot the points in the lightcurve not used in the fit.
plotQuadratic <- FALSE # set to TRUE to plot a quadratic fit
generateTS <- TRUE # set to TRUE to creat a time series from the data
tsBinWidth <- 10.0 # time series bin width in days. Important if generateTS is TRUE
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
marsOrder <- 45
marsPenalty <- 2 # set to 0 to avoid penalizing knots in pruning pass
marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "backward" # set to "none" to avoid pruning
splineRaw <-  FALSE # do the spline on the raw lightcurve, not binned.
##############################
okComparison <- "(000-?BLS-?556)|(000-?BLS-?551)|(000-?BLS-?553)|(000-?BLS-?552)|(000-?BLS-?554)|(000-?BLS-?549)|(000-?BLS-?555)|(108)|(113)|(116)|(118)|(121)|(124)|(128)|(ENSEMBLE)|(APASS20062365[+-]442738)" # regular expression from AAVSO photometry table

editUser <- data.frame(obsCode= "UJHA",startJD=2457650,endJD=2457760,band="V",stringsAsFactors=FALSE)
editUser <- rbind(editUser,c("PXR",2457512,2457513,"V"),c("JSJA",2457646,2457647,"V"),c("SGEA", 2457879,2457880,"V"))
# fill in some missing airmass values
lasCruces <- c(32.31994,-106.763654) # center of Las Cruces, NM in decimal degrees latitude, longitude.
Leominster <- c(52.226529,-2.741) # approximate location of PXR
tabbysLoc <- c("+44d 27m 24.61s","20h 06m 15.457s") # right ascension and declination of the star.
missingAirmass <- c("JM","PXR")	# observer code
missingAMLocs = rbind(lasCruces,Leominster)
### option to draw a vertical date line
drawDateLine <-  TRUE
jdLine <- c(2457892.0, 2457917.5)
jdLineColor <- "red"
jdLineText <- c("18May17","11Jun17")
##########################################################################
