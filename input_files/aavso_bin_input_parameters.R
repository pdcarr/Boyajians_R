##################
# input parameters for binning and plotting AAVSO data
llightcurve_name <- "data/aavso_latest_data.csv"
################# the name of the asassn data file
asassn.csv.file <- "data/ASASSN_latest.csv"
asassn.cameras <- c("bd","bt","bD")
# asassn.cameras <- c("bt")
merge.asassn <- TRUE # merge asasn data into V light curve
convert.asassn <- TRUE # TRUE if converting ASASSN g to V using Jordi, et. al. forumlas
asassn.bands <- c("g","SG","V")
our.BminusV <- 0.52 # B-V for our star in question
converted.V.bias <- 0.031 # additional bias to apply to subtract from converted V observations.
asassn.code <- "ASASSN"
######### filters #########################
maxairmass <- 2.5 # air mass values above this will be filtered out, as well as missing air masses. Set >= 100 to turn this off
maxuncertainty <- 0.1 # maximum AAVSO uncertainty estimate
maxBinUncertainty <- 0.025 # worst standard deviation to accept for a binned set of observations
wildsd <- 100.0 # worst number of standard deviations from mean allowed
##########################################
earliestJD = 2457290 # only data on or after this JD will be used
#earliestJD = 2458000 # only data on or after this JD will be used
#earliestJD <- 2457700
startPlot <- earliestJD
#startPlot <- 2457690
#startPlot <- 2457980
#startPlot <- 2457880
#startPlot <- 2457580
#startPlot <- 2457990
startPlot <- 2458150
#startPlot <- 2458000
#startPlot <- 2457620
# startPlot <- 2458200
# startPlot <- 2458700
#stop.plot <- 2457680
#stop.plot <- 2458000
add.predict <- 0 # additional prediction days to add to plot (experimental feature)
plotRelTimes <- TRUE
##########
includeExclude <- TRUE # TRUE if your list of observer codes is to to be included, FALSE if excluded or not used
ExclCodes <- "None"
weightless <- NA
# ExclCodes <- c("LDJ","DUBF","PXR","DKS","OJJ","SDB","VBPA","OAS","MJB","MATA",
# 				"WROC","MAND","VBPA","NOT","PALE","GKA","AMID","SGEA","GCJ","LBG","HJW","OAR","ASASSN","MMAO","HBB","EEY","MNIC","KHAB",
# 				"FRGA","BJFB","PTFA","TRE","ATE","DFS","FJAA","CIVA") # V & B
# ExclCodes <- c("LDJ","ASASSN")
# ExclCodes <- c("LDJ","DUBF","ASASSN","OAR")
# # ExclCodes <- c("DUBF","LDJ","GKA","CIVA","NOT")
# # ExclCodes <- c("DUBF")
# #ExclCodes <- c("VMT")
# ExclCodes <- c("ASASSN")
#ExclCodes <- c("LDJ","OAR")
# ExclCodes <- c("LDJ","OAR","ASASSN")
ExclCodes <- c("LDJ","ASASSN","OAR","HBB","DUBF","EEY","DJED","VMT","STFB",
				"TRE","BJFB","NOT","ATE","DFS","TIA","FJAA","CIVA","FRGA",
				"DJED","EEY","BSM","GKA","ODEA","LPAC","DKS","SFLB","RZD","PTFA","SGEA") # good small V band ensemble
weightless <- c("STFB","JM","FRGA","KHAB","TIA","FJAA","DKS","SFLB","RZD","LPAC","SGEA")  # weightless for V band
#ExclCodes <- c("LDJ","OAR","DKS","HBB","SGEA","HJW","TIA","DFS","FJAA") # new B ensemble under development
#ExclCodes <- c("LDJ","OAR","DKS","HBB","SGEA","ASASSN","OAR","EEY","DUBF") # merged B and V
#ExclCodes <- c("ASASSN")
# ExclCodes <- c("DUBF","GKA","BPAD","SJAR","LBG","LDJ","LWHA","SGEA","SDB","RNL","NOT","DFS","TIA","DFS","CIVA","DJED") # R ensemble
# weightless <- c("RNL","TIA","DJED","LWHA","VMT","JM","NOT","DFS") # weightless for R band
#ExclCodes <- c("DUBF","GKA","BPAD","SJAR","LBG","LDJ","LWHA") # R ensemble without JM
# ExclCodes <- c("OAR","OJJ","GKA","MJB","SJAR","LBG","LPB","LDJ","CMP","JM","VMT","SDB","SGEA","LPAC","TIA","TRE","BSM","FJQ") # I ensemble
# weightless <- c("TIA","FJQ") # I band weightless
#ExclCodes <- c("OAR","OJJ","GKA","MJB","SJAR","LWHA","LBG","LPB","LDJ","CMP","VMT","SDB","SGEA") # I ensemble without JM
# ExclCodes <- c("DUBF","MJB","LDJ","GKA","ELYA","HJW","JSJA","VBPA","DKS","OAR","HBB",
#                "SGEA","SDB","NOT","DFS","BJFB","DJED","GJP","SFLB","LPAC") # B ensemble
# weightless <- c("JM","TIA","FJAA","DJED","GJP","LPAC") # B band weightless
#ExclCodes <- c("JM","LDJ","DUBF")
#ExclCodes <- c("DUBF","LDJ","SGEA")
#ExclCodes <- "JM"

#ExclCodes <- "None"
plotMee <- NA # do not highlight any particular observer code
# plotMee <- "CIVA"
#plotMee <- "JM"
#plotMee <- "VMT"
#plotMee <- "MMAO"
#plotMee <- "FRGA"
meeColor <- "darkviolet"
########
allBands <- data.frame(bandinQ=c("I","R","V","B","SG"),plotColor=c("darkviolet","red","green","blue","aquamarine2"), stringsAsFactors=FALSE)
allBands <- data.frame(bandinQ=c("V"),plotColor=c("darkgreen"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V","B","SG"),plotColor=c("green","blue","aquamarine"), stringsAsFactors=FALSE)
# allBands <- data.frame(bandinQ=c("B"),plotColor=c("blue"), stringsAsFactors=FALSE)
# allBands <- data.frame(bandinQ=c("R"),plotColor=c("red"), stringsAsFactors=FALSE)
# allBands <- data.frame(bandinQ=c("I"),plotColor=c("darkviolet"), stringsAsFactors=FALSE)
# allBands <- data.frame(bandinQ=c("R"),plotColor=c("red"), stringsAsFactors=FALSE)
# allBands <- data.frame(bandinQ=c("SG"),plotColor=c("aquamarine2"),stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V","SG"),plotColor=c("darkgreen","aquamarine2"),stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R"),plotColor=c("darkviolet","red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R","B"),plotColor=c("darkviolet","red","blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("R","V","B"),plotColor=c("red","green","blue"), stringsAsFactors=FALSE)

########################
#deltaJD <- 2.0 # bin width in days
trial.bin <- 900 # trial number of bins. The acutal number will be less than this.
min.population <- 1 # smallest population allowed for 1 bin.
########################
plotExcluded <- FALSE # set to TRUE to plot the points in the lightcurve not used in the fit.
plotQuadratic <- FALSE # set to TRUE to plot a quadratic fit
generateTS <- FALSE # set to TRUE to create a time series from the data
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

############################## Comparison Stars
#okComparison <- "(000-?BLS-?556)|(000-?BLS-?551)|(000-?BLS-?553)|(000-?BLS-?552)|(000-?BLS-?554)|(000-?BLS-?549)|(BLS-549)|(BLS-555)|(000-?BLS-?555)|(108)|(113)|(116)|(118)|(121)|(124)|(128)|(ENSEMBLE)|(APASS20062365[+-]442738)|(3162-1853-1)" # regular expression from AAVSO photometry table
okComparison <- "."
bad.comp.star <- "(000-?BLS-549)|(000-?BLS-?549)|(BLS-549)"
#bad.comp.star <- NA

##########################################################################
pretty.interval <- 100 # controls rounding of plot limits
pretty.JD.interval <- 100

############## ensemble average plots?
plot.ensemble = TRUE

########## option to draw a vertical date line
drawDateLine <-  TRUE

###################### set to FALSE to turn off relative static biases
use.static.biases <- TRUE

######################### data from defined dips in dip_mask.R is weightless if this is TRUE
mask.Dips <- TRUE