cleanLC <- function(lightcurve,band,ExclCode,maxair,maxuncertainty,wildsigma,earliestJD) {
	# this function returns a binary vector specfiying which records (rows) in a lightcurve fram to use. 
	# lightcurve is the AAVSO data frame 
	# band is the string for the Band of interest: "I","R","V", and "B" are the most common.
	# you can optionally exclude one observer by specifying his or her observer code as  ExclCode
	# maxair is the maximum acceptable Airmass set to >= 100 to not test for airmass
	# max uncertainty is maximum acceptable reported uncertainty in magnitudes. Set to >=10 to not test for this
	# wildsigma is for filtering wild points. Set to the number of standard deviations you want to edit out.

# lightcurve is a data frame with AAVSO data	
	runningClean <- !is.na(lightcurve$JD) & !is.na(lightcurve$Magnitude)
	
	#exclude the observer in question
	exclObs <- lightcurve$Observer_Code != ExclCode
	runningClean <- runningClean & exclObs
	
	# which observations are in the specified band?
	Iband <- (lightcurve$Band == band)
	runningClean <- runningClean & Iband 
#	print(unique(runningClean))
	
	# which observations have airmass data and acceptable airmass?
	if (maxair < 100) {
		okair <- !is.na(lightcurve$Airmass)
		okair <- okair & (lightcurve$Airmass <= maxair)
		runningClean <- runningClean & okair
	}
#	print(unique(runningClean))
	
	# which observations have uncertainty data and acceptable uncertainty?
	okuncertainty <- !is.na(lightcurve$Uncertainty) & (lightcurve$Uncertainty > 0)
	if (maxuncertainty < 10) {
		okuncertainty <- okuncertainty  & lightcurve$Uncertainty <= maxuncertainty
		
	} 
	runningClean <- runningClean & okuncertainty
#	print(unique(runningClean))
		
	
	# calculate mean and sd over the cleaned set of magnitudes
	Imean = mean(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)
	Istdev = sd(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)

	# Remove anything bigger than max sd from mean
	# set maxsd very large if you don't want to do this
	notwild  <- (abs(lightcurve$Magnitude - Imean) <= wildsigma*Istdev) 
	runningClean <- runningClean & notwild
#	print(unique(runningClean))

# which Julian dates are on or after them minimum date?
	notEarly <- lightcurve$JD >= earliestJD
	runningClean <- runningClean & notEarly
	return(runningClean)
}

###########################################################################################


cleanAAVSO <- function(lightcurve,band,ExclCodes,maxair,maxuncertainty,wildsigma,earliestJD) {
	# this is an improved version of cleanLC
	# this function returns a binary vector specfiying which records (rows) in a lightcurve fram to use. 
	# lightcurve is the AAVSO data frame 
	# band is the string for the Band of interest: "I","R","V", and "B" are the most common.
	# you can optionally exclude one observer by specifying his or her observer code as  ExclCode
	# maxair is the maximum acceptable Airmass set to >= 100 to not test for airmass
	# max uncertainty is maximum acceptable reported uncertainty in magnitudes. Set to >=10 to not test for this
	# wildsigma is for filtering wild points. Set to the number of standard deviations you want to edit out.

# lightcurve is a data frame with AAVSO data	
	# clean any missing data
	runningClean <- !is.na(lightcurve$JD) & !is.na(lightcurve$Magnitude)
	
#	print(unique(runningClean))
	
	#exclude the observers in question
	for (ocode in ExclCodes) {
		exclObs <- lightcurve$Observer_Code != ocode
		runningClean <- runningClean & exclObs
	}
#	print(unique(runningClean))
	
	# which observations are in the specified band?
	Iband <- (lightcurve$Band == band)
	runningClean <- runningClean & Iband 
#	print(unique(runningClean))
	
	# which observations have airmass data and acceptable airmass?
	if (maxair < 100) {
		okair <- !is.na(lightcurve$Airmass)
		okair <- okair & (lightcurve$Airmass <= maxair)
		runningClean <- runningClean & okair
	}
#	print(unique(runningClean))
	
	# which observations have uncertainty data and acceptable uncertainty?
	okuncertainty <- !is.na(lightcurve$Uncertainty) & (lightcurve$Uncertainty > 0)
	if (maxuncertainty < 10) {
		okuncertainty <- okuncertainty  & lightcurve$Uncertainty <= maxuncertainty
		
	} 
	runningClean <- runningClean & okuncertainty
#	print(unique(runningClean))
		
	
	# calculate mean and sd over the cleaned set of magnitudes
	Imean = mean(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)
	Istdev = sd(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)

	# Remove anything bigger than max sd from mean
	# set maxsd very large if you don't want to do this
	notwild  <- (abs(lightcurve$Magnitude - Imean) <= wildsigma*Istdev) 
	runningClean <- runningClean & notwild
#	print(unique(runningClean))

# which Julian dates are on or after them minimum date?
	notEarly <- lightcurve$JD >= earliestJD
	runningClean <- runningClean & notEarly
	
#	print(unique(runningClean))
	
	return(runningClean)
}

###########################################################################################


cleanAAVSO2 <- function(lightcurve,band,ExclCodes,inclexcl,maxair,maxuncertainty,wildsigma,earliestJD) {
	# this is an improved version of cleanAAVSO that allows the calling script to determine whether or include or exclude observers.
	# this function returns a binary vector specfiying which records (rows) in a lightcurve fram to use. 
	# lightcurve is the AAVSO data frame 
	# band is the string for the Band of interest: "I","R","V", and "B" are the most common.
	# you can optionally exclude OR include a subset of observer codes by specifying ExclCodes
	# if inclexcl is TRUE, then include only observer codes listed in ExclCodes, otherwise exclude those observers
	# maxair is the maximum acceptable Airmass set to >= 100 to not test for airmass
	# max uncertainty is maximum acceptable reported uncertainty in magnitudes. Set to >=10 to not test for this
	# wildsigma is for filtering wild points. Set to the number of standard deviations you want to edit out.

# lightcurve is a data frame with AAVSO data	
	# filter out any missing data
	runningClean <- !is.na(lightcurve$JD) & !is.na(lightcurve$Magnitude) &!is.na(lightcurve$Observer_Code)
#	print(unique(runningClean))
	
	# include exclude the observers in question
	firstCode <- TRUE
	for (ocode in ExclCodes) {
		if (inclexcl == FALSE) {
			# exclude the listed observers
			if (firstCode) {
				exclObs <- lightcurve$Observer_Code != ocode
			} else {
				exclObs <- exclObs & (lightcurve$Observer_Code != ocode)
			}
		} else {
			# include the listed observer codes
			if(firstCode) {
				exclObs <- lightcurve$Observer_Code == ocode
			} else {
				exclObs <- exclObs | lightcurve$Observer_Code == ocode
			}
		}
		firstCode <- FALSE	
	}
	runningClean <- runningClean & exclObs
#	print(unique(runningClean))
	
	# which observations are in the specified band?
	Iband <- (lightcurve$Band == band)
	runningClean <- runningClean & Iband 
#	print(unique(runningClean))
	
	# which observations have airmass data and acceptable airmass?
	if (maxair < 100) {
		okair <- !is.na(lightcurve$Airmass)
		okair <- okair & (lightcurve$Airmass <= maxair)
		runningClean <- runningClean & okair
	}
#	print(unique(runningClean))
	
	# which observations have uncertainty data and acceptable uncertainty?
	okuncertainty <- !is.na(lightcurve$Uncertainty) & (lightcurve$Uncertainty > 0)
	if (maxuncertainty < 10) {
		okuncertainty <- okuncertainty  & lightcurve$Uncertainty <= maxuncertainty
		
	} 
	runningClean <- runningClean & okuncertainty
#	print(unique(runningClean))
		
	
	# calculate mean and sd over the cleaned set of magnitudes
	Imean = mean(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)
	Istdev = sd(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)

	# Remove anything bigger than max sd from mean
	# set maxsd very large if you don't want to do this
	notwild  <- (abs(lightcurve$Magnitude - Imean) <= wildsigma*Istdev) 
	runningClean <- runningClean & notwild
#	print(unique(runningClean))

# which Julian dates are on or after them minimum date?
	notEarly <- lightcurve$JD >= earliestJD
	runningClean <- runningClean & notEarly
	
#	print(unique(runningClean))
	
	return(runningClean)
}


###########################################################################################

cleanAAVSO3 <- function(lightcurve,band,ExclCodes,inclexcl,maxair,maxuncertainty,wildsigma,earliestJD,okCompStars,bad.stars=NA) {
	# this is an improved version of cleanAAVSO that allows the calling script to determin whether or include or exclude observers.
	# this function returns a binary vector specfiying which records (rows) in a lightcurve fram to use. 
	# lightcurve is the AAVSO data frame 
	# band is the string for the Band of interest: "I","R","V", and "B" are the most common.
	# you can optionally exclude OR include a subset of observer codes by specifying ExclCodes
	# if inclexcl is TRUE, then include only observer codes listed in ExclCodes, otherwise exclude those observers
	# maxair is the maximum acceptable Airmass set to >= 100 to not test for airmass
	# max uncertainty is maximum acceptable reported uncertainty in magnitudes. Set to >=10 to not test for this
	# wildsigma is for filtering wild points. Set to the number of standard deviations you want to edit out.
	# earliestJD is the earliest Julian date to allow
	# okCompstars is a regular expression that should match any valid comparison star
	# bad.stars is a regular expression that should match bad comparison stars
# debug stuff
#	mnic <- lightcurve$Observer_Code == "MNIC"
#	print(sum(mnic))
# lightcurve is a data frame with AAVSO data	
	# filter out any missing data
	runningClean <- !is.na(lightcurve$JD) & !is.nan(lightcurve$JD) & 
					!is.na(lightcurve$Magnitude) & !is.nan(lightcurve$Magnitude)
					!is.na(lightcurve$Observer_Code) & 
					!is.na(lightcurve$Comp_Star_1)
	
#	print(unique(runningClean))
#	print(sum(mnic & runningClean)) #debug print
	# include exclude the observers in question
	firstCode <- TRUE
	for (ocode in ExclCodes) {
		if (inclexcl == FALSE) {
			# exclude the listed observers
			if (firstCode) {
				exclObs <- lightcurve$Observer_Code != ocode
			} else {
				exclObs <- exclObs & (lightcurve$Observer_Code != ocode)
			}
		} else {
			# include the listed observer codes
			if(firstCode) {
				exclObs <- lightcurve$Observer_Code == ocode
			} else {
				exclObs <- exclObs | lightcurve$Observer_Code == ocode
			}
		}
		firstCode <- FALSE	
	}
	runningClean <- runningClean & exclObs
#	print(unique(runningClean))
	
	# which observations are in the specified band?
	Iband <- (lightcurve$Band == band)
	runningClean <- runningClean & Iband 
#	print(unique(runningClean))

	
	# which observations have airmass data and acceptable airmass?
	if (maxair < 100) {
		okair <- !is.na(lightcurve$Airmass)
		okair <- okair & (lightcurve$Airmass <= maxair)
		runningClean <- runningClean & okair
	}
#	print(unique(runningClean))
	
	# which observations have uncertainty data and acceptable uncertainty?
	okuncertainty <- !is.na(lightcurve$Uncertainty) & (lightcurve$Uncertainty > 0)
	if (maxuncertainty < 10) {
		okuncertainty <- okuncertainty  & lightcurve$Uncertainty <= maxuncertainty
		
	} 
	runningClean <- runningClean & okuncertainty
#	print(unique(runningClean))
		
	
	# # calculate mean and sd over the cleaned set of magnitudes
	# Imean = mean(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)
# #	print(Imean)
	# Istdev = sd(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)

	# # Remove anything bigger than max sd from mean
	# # set maxsd very large if you don't want to do this
	# notwild  <- (abs(lightcurve$Magnitude - Imean) <= wildsigma*Istdev) 
	# runningClean <- runningClean & notwild
# #	print(unique(runningClean))

# which Julian dates are on or after them minimum date?
	notEarly <- lightcurve$JD >= earliestJD
	runningClean <- runningClean & notEarly
	
#	print(unique(runningClean))

# which observers are not reporting ok comparison star 1 (what is "ensemble"?)
	starsOK <- !is.na(lightcurve$Comp_Star_1) & !is.na(lightcurve$Comp_Star_2) # check that a comparison star is provided
	blankMask <- starsOK & FALSE 	# creates a matching mask of indices all FALSE
	okIndex  <- grep(okComparison,as.character(lightcurve$Comp_Star_1),ignore.case=TRUE) # check the first comparison star
	if (length(okIndex) != 0) {blankMask[okIndex] <- TRUE} # checks for an unlikely circumstance
	starsOK <- starsOK & blankMask
	blankMask <- starsOK & FALSE
	okIndex  <- grep(okComparison,as.character(lightcurve$Comp_Star_2),ignore.case=TRUE) # check the second comparison star
	if (length(okIndex) != 0) {blankMask[okIndex] <- TRUE} # checks for an unlikely circumstance
	starsOK <- starsOK & blankMask
	runningClean <- runningClean & starsOK
#	print(length(runningClean))
#	print(sum(runningClean))

	
 # are there any bad comparison stars we need to watch out for?
	if (!is.na(bad.stars)) {
		blankMask <- starsOK | TRUE 	# creates a matching mask of indices all TRUE
		these.are.bad <- grep(bad.stars,as.character(lightcurve$Comp_Star_1,ignore.case=TRUE)) # check if the first comperison star is bad
		if (length(these.are.bad) != 0) {blankMask[these.are.bad] <- FALSE} # checks for an unlikely circumstance
		starsOK <- starsOK & blankMask
		blankMask <- starsOK | TRUE 	# creates a matching mask of indices all TRUE
		these.are.bad <- grep(bad.stars,as.character(lightcurve$Comp_Star_2,ignore.case=TRUE)) # check if the second comperison star is bad
		if (length(these.are.bad) != 0) {blankMask[these.are.bad] <- FALSE} # checks for an unlikely circumstance
		starsOK <- starsOK & blankMask

		runningClean <- runningClean & starsOK
	}
#done
#	print(length(runningClean))
#	print(sum(runningClean))

	return(runningClean)
}


###########################################################################################

cleanLC_oneObs <- function(lightcurve,band,ExclCode,maxair,maxuncertainty,wildsigma) {
# this function picks out a single observer and returns a data frame for just that observer
#lightcurve is a data frame with AAVSO data	

	#exclude all other observers except the one in question
	runningClean <- !is.na(lightcurve$JD) & !is.na(lightcurve$Magnitude)
	exclObs <- lightcurve$Observer_Code == ExclCode
	#omitted <- lightcurve[exclObs,]
	runningClean <- runningClean & exclObs
#	print(unique(runningClean))
	
	if (length(lightcurve[runningClean,"JD"]) > 0) {
		
		# which observations are in the correct band?
		Iband <- (lightcurve$Band == band) 
		runningClean <- runningClean & Iband
#	print(unique(runningClean))
		
		# which observations have airmass data and acceptable airmass?
		if (maxair < 100) {
			okair <- !is.na(lightcurve$Airmass)
			okair <- okair & lightcurve$Airmass <= maxair
			runningClean <- runningClean & okair
		}
#	print(unique(runningClean))
		
		# which observations have uncertainty data and acceptabel uncertainty?
		
		if (maxuncertainty < 10) {
			okuncertainty <- !is.na(lightcurve$Uncertainty) & (lightcurve$Uncertainty > 0)
			okuncertainty <- okuncertainty  & lightcurve$Uncertainty <= maxuncertainty
			runningClean <- runningClean & okuncertainty
		}
#	print(unique(runningClean))
			
		
		# calculate meand and sd over the cleaned set of magnitudes
		if(length(lightcurve[runningClean,"Magnitude"]) > 2) {
			Imean = mean(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)
			Istdev = sd(lightcurve[runningClean,"Magnitude"],na.rm=TRUE)
			print(toString(Istdev))
			# Remove anything bigger than max sd from mean and excessive airmass
			notwild  <- (abs(lightcurve$Magnitude - Imean) <= wildsigma*Istdev)
			runningClean <- runningClean & notwild 
		}
#	print(unique(runningClean))
	
	} 
	return(runningClean)
} 
 
##############################################################################################

LClm <- function(cleancurve,lmweights) {
	myIfit <- lm(cleancurve$Magnitude ~ cleancurve$JD, weights= lmweights)
	return(myIfit)
}

##############################################################################################

AAVSOObsStats <- function(lightcurve,cleanBands,allBands) {
# 
# lightcurve is a data frame extracted from an AAVSO data download.
# cleanBands is the cleaning matrix, with each column corresponding to a passband.
# allBands is a data frame with names of all the bands of interest and their plot colors.

myobscounts <- list(length(allBands$bandinQ))


# logical or all theclean bands together
firstBand <- TRUE

for (thisBand in 1:length(cleanBands[1,])) {
	if (firstBand) {
		goodObs <- cleanBands[,thisBand]
		firstBand <- FALSE
	} else {
		goodObs <- goodObs | cleanBands[,thisBand]
	}
}

weare <- unique(lightcurve[goodObs,"Observer_Code"]) # all the observers with good Observations
#print(length(weare))

obscounts <- matrix(nrow=1,ncol=length(allBands$bandinQ)) # set up the matrix we will return
myobscounts <- matrix(nrow=1, ncol=length(allBands$bandinQ)) #set up the temporary matrix
bandCount <- length(allBands$bandinQ)	# number of passband we are using

#loop over all the observer codes
for (myrow in 1:length(weare)) {
		thisO <- lightcurve$Observer_Code == weare[myrow]
		#loop over each passband and get counts for the passband for the observer
		for (thisBand in 1:bandCount) {
			myTest <- (lightcurve$Band == allBands$bandinQ[thisBand] & thisO)
			myobscounts[1,thisBand]  <-  length(lightcurve[myTest & cleanBands[,thisBand],"Band"])
		}
	if (myrow ==1) {
		obscounts[myrow,] <- myobscounts
	} else {
		obscounts <- rbind(obscounts,myobscounts)
	}
	}
#	create a data frame with the observer codes and the counts for each passband
	obsAndCounts <- data.frame(obscode = weare,countsMatrix = obscounts)
	colnames(obsAndCounts)[2:(1+bandCount)] <- allBands$bandinQ
	return (obsAndCounts) # a a data frame with observer codes and counts
}

####################################################### binAAVSO
binAAVSO <-  function(lightcurve,cleanObs,allBand,deltaJD=1,weightless=NA,trial.bins=100,minimum.membership=1) {
# returns a dataframe with JD, Magnitude for each band, and observer code.
# lightcurve is the unprocessed curve from AAVSO
# cleanObs is the matrix of logical vectors of accepted observations for the passbands
# allBands is the dataframe containing the AAVSO codes of the passbandsa, e.g. "B","V"
# deltaJD is the time cluster size in days
# weightless is the list of observer codes who are given zero weight in the fit or ensemble
# trial.bins is an interger number of bins to go for initially. This may be reduced, but won't be increased.


	allClean <- cleanObs[,1]
	if (ncol(cleanObs) > 1) {
		for (n in 2:ncol(cleanObs)) {allClean <- allClean | cleanObs[,n]}
	}

	#calculate the start and stop times from the light curve Julian Dates 
	
#	startJD = floor(min(lightcurve[allClean,"JD"],na.rm=TRUE))
#	stopJD = ceiling(max(lightcurve[allClean,"JD"],na.rm=TRUE))
#	print(length(lightcurve[allClean,"JD"]))
#	print(startJD)
#	print(stopJD)
#	print(deltaJD)
	weare = unique(lightcurve[allClean,"Observer_Code"]) # a list of all the observer codes
	numBands = length(allBand$bandinQ) # the number of passbands
#	right.colors <- lightcurve$Band[allClean] %in% allBand$bandinQ # logical vector every observation that has the right band(s)
	#set up the data frames we'll be populating - one scratch, and one is the output data frame
	allSuperObs <- data.frame(JD=numeric(),Band=character(),Magnitude=numeric(),Uncertainty=numeric(),nobs=numeric(),Observer_Code=character(),stringsAsFactors=FALSE)
	superObs <- data.frame(JD=numeric(),Band=character(),Magnitude=numeric(),Uncertainty=numeric(),nobs=numeric(),Observer_Code=character(),stringsAsFactors=FALSE)
	ensemble.means <- data.frame(JD = numeric(),Band=character(),Magnitude=numeric(),stringsAsFactors=FALSE)
	ensemble.entry <- data.frame(JD = numeric(),Band=character(),Magnitude=numeric(),nobs=numeric(),stringsAsFactors=FALSE)
	
# determine bin times and membership for the allowed observations (cleaned lightcurve with the required colors)
	bin.defs <- kmeans.time.series(lightcurve$JD[allClean],trial.bins,min.population=minimum.membership,delta.mean=0.04,max.iterations=12) 
		#loop over the bin times
	nbins <- length(bin.defs$bins)
	print(paste("kmeans found",nbins,"bins")) #debug
#	koope over the bins we just defines
	for(this.bin in 1:nbins) {
#			print(mean.bin.time)
#			print(stopNow)
		# print(ensemble.mean)
		# print(accumulated.sups)
		
		#loop over the passbands to create the superobservation in the time frame for the Observer code
		for (bandIndex in 1:numBands) {			
#			ensemble.test  <- testTime & cleanObs[,bandIndex]
		  # to get in the ensemble for a bin, the observation has to be in the right band and in the bin.
		  ensemble.test <- bin.defs$membership == this.bin & (lightcurve$Band[allClean] == allBand$bandinQ[bandIndex])
			# loop over observers
			for (thisObs in weare) {
				testObs <- lightcurve$Observer_Code[allClean] == thisObs # test for observer in question
				allTests <- ensemble.test & testObs # all observations eligible for super observation
				sublight <- lightcurve[allClean,]
				# compile the super observation
				n <- sum(allTests) # number of observations in the bin for this observer
				if(n > 0) {
					superObs[1,"JD"] <- mean(sublight$JD[allTests],na.rm=TRUE)
					superObs[1,"Band"] <- allBand$bandinQ[bandIndex]
					superObs[1,"Magnitude"] <- mean(sublight$Magnitude[allTests],na.rm=TRUE)
					superObs[1,"Uncertainty"] <- mean(sublight$Uncertainty[allTests],na.rm=TRUE)/sqrt(n)
					superObs[1,"nobs"] <- n
					if(is.na(superObs[1,"Uncertainty"]) | is.nan(superObs[1,"Uncertainty"])) {
						print("invalid uncertainty - not a number")
						next
					}
					if(superObs[1,"Uncertainty"] < 0) {
						print("invalid uncertainty")
						print(superObs)
						
						next
					} #invalid uncertainty
						
					superObs[1,"Observer_Code"] <- thisObs # add the observer code to the superobservation
										
	#					append the superobservation to the main data frame using rbind
					allSuperObs <- rbind(allSuperObs,superObs)
				} # if(n> 0) if statement
			} #observer loop end
			
			ensemble.test <- ensemble.test &  
			              	!(sublight$Observer_Code %in% weightless) # remove weightless observers from ensemble
			if(sum(ensemble.test) >0) {
					ensemble.mean <- mean(sublight$Magnitude[ensemble.test])  # unweighted average
					mean.time <- mean(sublight$JD[ensemble.test],na.rm=TRUE)
					ensemble.entry[1,"JD"] <- mean.time
					ensemble.entry[1,"Band"] <- allBand$bandinQ[bandIndex]
					ensemble.entry[1,"Magnitude"] <- ensemble.mean
					ensemble.entry[1,"nobs"] <- sum(ensemble.test)
					ensemble.means <- rbind(ensemble.means,ensemble.entry)
			}
		} # band loop end
	} #time bins
	return(list(allSuperObs,ensemble.means))
}

#######################################################################
binAAVSO_ts <-  function(lightcurve,cleanObs,allBand,deltaJD) {
# returns a time_series (ts) with Magnitude for each band, and observer code.
# lightcurve is the unprocessed curve
# cleanObs is the matric of logical vectors of accepted observations for the passbands
# allBands is the dataframe containing the AAVSO codes of the passbandsa, e.g. "B","V"
# deltaJD is the time cluster size in days

	library(stats)

	allClean <- cleanObs[,1]
	if (ncol(cleanObs) > 1) {
		for (n in 2:ncol(cleanObs)) {allClean <- allClean | cleanObs[,n]}
	}

	#calculate the start and stop times from the llight curve Julian Dates 
	startJD <- floor(min(lightcurve[allClean,"JD"],na.rm=TRUE))
	stopJD <- ceiling(max(lightcurve[allClean,"JD"],na.rm=TRUE))
	
	numBands <- length(allBand$bandinQ) # the number of passbands
	uNames <- sapply(as.list(allBand$bandinQ),paste,"Uncertainty")
	myNames <- c(allBand$bandinQ,uNames)	
	
	#set up the matrix we'll be populating
	allSuperObs <- matrix(ncol=2*numBands)
	#loop over the times
	firstObs = TRUE	
	for (startNow in seq(startJD,stopJD,by=deltaJD)) {
		stopNow <- startNow+deltaJD
		testTime <- (lightcurve$JD >= startNow) & (lightcurve$JD < stopNow) 
	#loop over the passbands
		superObs <- vector(length=2*numBands)
		for (bandIndex in 1:numBands) {
			allTests <- testTime & cleanObs[,bandIndex]
			# compile the super observation
#			print(length(lightcurve$JD[allTests]))
			if(length(lightcurve[allTests,"JD"]) > 0) {
				superObs[bandIndex] <- mean(lightcurve[allTests,"Magnitude"],na.rm=TRUE)
				if (length(lightcurve[allTests,"JD"]) >=4 ) {
					superObs[numBands+bandIndex] <- sd(lightcurve[allTests,"Magnitude"],na.rm=TRUE)
				} else {
					superObs[numBands+bandIndex] <- mean(lightcurve[allTests,"Uncertainty"],na.rm=TRUE)
				}
				if(superObs[numBands+bandIndex] <= 0) {superObs[numBands+bandIndex] <- NA} #invalid uncertainty
#					print(superObs)
			} else {
#				print(c(bandIndex," ",startNow))
#				print(length(lightcurve$JD[cleanObs[,bandIndex]]))
#				print(length(lightcurve$JD[testTime]))
#				print(length(lightcurve$JD[allTests]))
				superObs[bandIndex] <- NA
				superObs[numBands + bandIndex] <- NA
			}
		}
#					append the superobservation to the main data frame using rbind
		if(firstObs) {
			allSuperObs[1,] <- superObs
			firstObs <- FALSE
		} else {
			allSuperObs <- rbind(allSuperObs,superObs)
		}
	}
#	print(myNames)
	return(ts(data=allSuperObs,start=startJD+deltaJD/2,deltat=deltaJD,names=myNames))
}


#######################################################################

eDate2JDate <- function(excelDate) {
	# returns Julian Date corresponding to an Excel numerical date.
	return(excelDate + 2415018.5)
}

#######################################################################
JDDate2Date <- function(jDate) {
	# converts a Julian Date to a Date() object.
	unixT <- (jDate - 2440587.5)*86400 # assumes 1970 origin
	mylt <- strptime(as.POSIXct(unixT,origin="1970-01-01 00:00:00",tz="UTC"),format="%Y-%m-%d",tz="UTC")
	return(as.Date(mylt))
}

#######################################################################
ctDate2JD <- function(ctDateTime) {
	# ctDateTime is a POSIXct class time. A Julian date is returned.
	jdTime = as.numeric(ctDateTime)/86400 + 2440587.5
	return(jdTime)
}

#######################################################################
ObserverJDEdit <- function(editFrame,lightcurve) {
# marks with FALSE all entries from an AAVSO light curve or binned curve that correspond to a particular observer, band and time range

	keepThis <- !is.na(lightcurve$Observer_Code)
	killThese <- is.na(lightcurve$JD)
	
	for (index in 1:length(editFrame$obsCode)) {
#		print(index)
#		print(editFrame$obsCode[index])
		killThis <- lightcurve$Observer_Code == editFrame$obsCode[index]
		killThis <- killThis & (lightcurve$JD >= editFrame$startJD[index] & lightcurve$JD <= editFrame$endJD[index])
#		print(editFrame$startJD[index])
		killThis <- killThis & (lightcurve$Band ==  editFrame$band[index])
		killThese <- killThese | killThis
	}
	keepThis <- keepThis & !killThese

	return (keepThis)
}

#######################################################################

aavsoJDbounds <- function(lightcurve,jdbias=0) {
	# takes an AAVSO light curve as an argument. Based upon last Julian Date in the file, recommends a set of time bounds for downloading new data
	# from AAVSO.org
	natest <- is.numeric(lightcurve$JD)
	jdEnd = max(lightcurve$JD[natest]) + 0.00001 + jdbias
	jdNew = ctDate2JD(Sys.time()) 
	return(c(jdEnd,jdNew))
}

#######################################################################
goldenSectionTest <-  function(low,high) {
	# for this to work right, high > low
	phi = (1+sqrt(5))/2 # golden ratio
	return(c(high - (high-low)/phi,low + (high-low)/phi))
}

###########################################################
filter.dips.JD <- function(lightcurve.JDs,dip.mask) {
	# lightcurve.times are the Julian dates for lightcurve observations
	# dip.mask is a data frame containing the dip names and start and stop times
	
	dip.index <- 1
	good.times <- !is.na(lightcurve.JDs)
	for (myname in dip.mask$dipname) {
		good.times <- good.times & (lightcurve.JDs < dip.mask$JD.begins[dip.index] | lightcurve.JDs > dip.mask$JD.ends[dip.index])
		dip.index <- dip.index + 1
	}
	return(good.times)
}

##############################################
# truncates a POSIX time to the day
trunc.pxct <- function(PXct.time) {
	scratch <- as.POSIXlt(PXct.time)
	scratch$hour <-0
	scratch$min <- 0
	scratch$sec <- 0
	return(as.POSIXct(scratch))
}

########################################
asassn.merge <- function(lightcurve,asassn.data,asassn.code="ASASSN",g.to.V=FALSE,star.BminusV=0.52,V.bias=0.0,asassn.cs=c("UNK","UNK")) {
# this function merges in the ASASSN data in both V and SG bands into the AAVSO lightcurve, appending it to the end.
# lightcurve is the AAVSO lightcurve read in from their .csv file
# asassn.data is the ASAS-SN data read in form the .csv file (https://asas-sn.osu.edu/)
# asassn.code is the observer Code you want ASASSN to have
# g.to.V controls converting g band to V
# star.BminusV is the B-V for the star in question
# V.bias is an additional bias to be subtracted from converted V values
# asassn.cs is the comp stars you want entered into the light curve.
	n <- nrow(asassn.data)
	m = ncol(lightcurve)
	scratch <- lightcurve[1,]
	if(g.to.V) {	cat("\nConverting any ASAS-SN g band magnitudes to V band per Jordi, et. al., (2005)\n")}
	for (index in 1:n) {
		scratch$JD <- asassn.data$HJD[index]
		if(is.numeric(asassn.data$mag[index])) {
			scratch$Magnitude <- asassn.data$mag[index]
		} else {
			scratch$Magnitude <- NaN
		}
		scratch$Observer_Code <- asassn.code
		scratch$Uncertainty <- asassn.data$mag_err[index]
		
		# all ASASSN data should be either "V" or "g" (which we think is g' or SG)
		if(asassn.data$Filter[index] == "g") {
			scratch$Band <- "SG"
			if(g.to.V & is.numeric(scratch$Magnitude)){
				scratch$Band <- "V"
				# see http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
				scratch$Magnitude <- scratch$Magnitude - 0.63*star.BminusV +0.124 - V.bias # per Jordi, et. al., (2005)
#				print(scratch$Magnitude)
			}
		} else if(asassn.data$Filter[index]=="V") {
			scratch$Band <- "V" 
		} else { 
			print(asassn.data$Filter) # debug
			next 
		}
			
		scratch$Comp_Star_1 <- asassn.cs[1]
		scratch$Comp_Star_2 <- asassn.cs[2]
		scratch$Comments <- "Merged ASASSN Data"
		scratch$Credit <- "Ohio State"
		scratch$Observer.Affiliation <- "Ohio State"
		scratch$Airmass <- 1.0
		
		scratch$HQuncertainty <- NA
		scratch$Charts <- NA
		scratch$Transfomed <- NA
		scratch$Validation.Flag <- NA
		scratch$Cmag <- NA
		scratch$Kmag <- NA
		scratch$Measurement.Method <- NA
#		print(scratch)
#		print(ncol(scratch))
#		print(ncol(lightcurve))
		# browser()
		lightcurve <- rbind(lightcurve,scratch)
		
	}
	return(lightcurve)	
}

####################### mean.resid ##########################
mean.resid <- function(binCurve,these.resids,btest,Obs.code,dipless) {
# calculates the mean residual for an observer code when not in a dip. Are the arguments are computed as a matter of course.
# binCurve is the binned lightcurve	
# these.resids are the residuals calculated relative to the fit and the smae length as a columnn of binCurve
# btest is a logical vector the same length as a column of binCurve
# Obs.code is a strong designating an observer, e.g. "OAR"
# dipless is the dip mask, and must be the same length as sum(btest)
	my.obs <- (binCurve$Observer_Code[btest] == Obs.code) & dipless
#	print()
	return(mean(these.resids[btest][my.obs]))
}

################### lightcurve.biases #########################
lightcurve.biases <- function(lightcurve,biasObserver) {
# best to use the cleaned lightcurve for efficiency	
	for(index in 1:nrow(lightcurve)) {
		lightcurve[index,] <- entry.bias.apply(lightcurve[index,],biasObserver)
		if(index %% 1000 == 0) {print(index)}
	}
	return(lightcurve)
}

entry.bias.apply <- function(lightcurve.row, biasObserver) {
#lightcurve.row is a row from athe lightcurve data frame
#bias observer is a data frame with the biases for a subset of observers in a subset of the bands
	this.row = as.data.frame(lightcurve.row)
	# print(class(this.row))
	find.bias <- biasObserver$band == this.row $Band & biasObserver$obsCode == this.row $Observer_Code
#	print(sum(find.bias)) # debug
	if(sum(find.bias) > 0) {
		this.bias <- as.numeric(biasObserver$bias[find.bias][1])
		# print(this.bias)
		this.row$Magnitude <- this.row$Magnitude - this.bias
	}
#	print(class(this.row))
	return(this.row)
}


####################################################################################################################
exclude.LCO <- function(LCO.curve,exclude.these,exclude=TRUE) {
## takes a light curve from LCO and retiurns a logical vector of observations to ignore/include
# LCO Curve is data.frame(Num = numeric(),label=character(),MJD=numeric(),rel_flux_T1_n=numeric(),band=character())
# exclude.these are Perl regula expressions to match
# exclude is a logical. TRUE if you wnat to exclude labels matching the regular expression, FALSE if you want to only include them

# return value is logical vector of itmes that can used to filter in/out whatever is in exclude.these
  ok.obs <- !is.na(LCO.curve$MJD) # creates a logical vector as long as the data frame.
  found <- ok.obs & FALSE # a vector of FALSE
  for (code in exclude.these) {
    found <- found | grepl(pattern=code,x=LCO.curve$label,perl=TRUE)
  }
  return(ok.obs & xor(found,exclude))
}

######################################################################

kmeans.time.series <- function(times,initial.clusters.num=16,min.population=1,delta.mean=0.04,max.iterations=12) {
## clusters a time series using a simple algorithm.
  # times is a vector of real numbers. NA or NaNs are ok.
  # initial.clusters is an integer number of clusters to start with
  # min.population is the minimum number of elements in the cluster else the cluster is pruned. St ot zero to avoid pruning.
  # delta.mean: the algorithm is doen when all the means change less than this number
  # max.iterations is the maximum interger number before we are declaring victory.
# returns a vector of bin centers, which may be smaller than initial.clusters
# error checking  
  # browser()
  if(length(times) <= 1 | 
     as.integer(initial.clusters.num) < 1 | 
     delta.mean <= 0 | 
     as.integer(max.iterations) <= 0) { return(NA)}
# set up the initial centers
  t.min <- min(times,na.rm=TRUE)
  t.max <- max(times,na.rm=TRUE)
  t.span <- t.max - t.min
  t.bin <- t.span/initial.clusters.num
  clusters <- seq(from=t.min + t.bin/2,to=t.max - t.bin/2,t.bin) # first guess at the clusters
  N <- length(clusters)
  old.cluster.mean <- rep(NA,times=N)
  for(iter.num in 1:as.integer(max.iterations)) {
    N <- length(clusters)
    pruned <- is.na(clusters) # logical vector same size as clusters
    cluster.sum = vector(mode="numeric",length=N)
    cluster.mean = vector(mode="numeric",length=N)
    cluster.diff = vector(mode="numeric",length=N)
    # browser()
    situation <- t(sapply(times,k.min.distance,bins=clusters))
    
    for (i.bin in 1:N) {
      my.times <- !is.na(situation[,1]) & situation[,1] == i.bin
      # browser()
      cluster.sum[i.bin] <- sum(my.times,na.rm=TRUE)
      if(cluster.sum[i.bin] >= min.population) {
        pruned[i.bin] <- FALSE
        cluster.mean[i.bin] <- mean(times[my.times],na.rm=TRUE)
        clusters[i.bin] <- cluster.mean[i.bin]
        if(!is.na(old.cluster.mean[i.bin])) {
          cluster.diff[i.bin] <- abs(old.cluster.mean[i.bin] - cluster.mean[i.bin] )
        } 
      } else {
#        browser()
       pruned[i.bin] <-  TRUE # probably not necessary, but makes sure.
       next()
      }
    } # bin loop
    old.cluster.mean <- cluster.mean[!pruned]
    clusters <- clusters[!pruned] # pruned clusters are gone and won't come back
    if(sum(pruned) > 0) {next()}  # need another go-around to pick up the orphaned times

#    if(iter.num == as.integer(max.iterations)) {print("exiting k means on max iterations")} # debug
    if(max(cluster.diff[!pruned],na.rm=TRUE) <= delta.mean & iter.num > 1) {
      # print("breaking loop on convergence")
      # print(iter.num)
      break
      } # Converged!

  } # maximum iteration loop
  situation <- t(sapply(times,k.min.distance,bins = clusters)) # bin membership mapping
  return(list(bins=clusters,membership=situation[,1]))
}

#################
k.min.distance <- function(time,bins) {
# for a given time, calculate the min distance from the distances to all the bins.
  if(is.na(time)) {return(c(NA,0))}
  L <- length(bins)
  distance <- vector(mode="numeric",length=L)
  # create a vector of distances
  distance <- abs(bins - time[1])
  
  # find the shortest distance
  bin.number <- 0 # indicates an error if it stays 0
  min.dist <- min(distance,na.rm=TRUE)
  bin.number <- which.min(distance)
  # browser()
  return(c(bin.number,min.dist))
}

bin.fluxes <- function(bin.number,situation,fluxes) {
  membership.element <- situation == bin.number
  bin.mean.flux <- mean(fluxes[membership.element],na.rm=TRUE)

}

#########################
# function to calcuate the standard error of the slope of a linear regression. You need at last three data points.
slope.SE <- function(x,residuals){
# x is the independent variable for the regression
# residuals are the residuals vs. the line
  nx <- sum(!is.na(x))
  n.resid <- sum(!is.na(residuals))
  if(n.resid <=2) {return(NA)} # makes sure we have enogh data points
  return(sqrt(var(residuals,na.rm=TRUE)*(n.resid-1)/(n.resid-2))/sqrt(var(x,na.rm=TRUE)*(nx -1)))
}