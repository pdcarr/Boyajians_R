# this is designed to plot AAVSO magnitudes after some data cleaning, with decreasing magnitude, and to fit the data to a straight line.

# housekeeping
#compile the function libraries
source("data_funcs.R") # a file of supporting functions
source("astro_funcs.R")
source("plot_funcs.R")
# remove old garbage if it's there
rm(allFits)
rm(cleanBand)
rm(binCurve)
# options
options(digits=12) # hard to read JDs without this setting
# load the required packages
library("MASS") # for rlm() and lqs()
library("smooth") # for smoothing
library("earth") # for MARS
library("crayon") # to add color to text

allFits <- list()


##################

# all the inputs are in the sourced file
source("input_files/aavso_bin_input_parameters.R")
source("input_files/observer_edits.R")
source("input_files/missing_airmass.R")
source("input_files/VlineParams.R")

if (includeExclude) {
	inclWord <- "used"
} else {
	inclWord <- "not used"
}

# load light curve
#col.classes = c("numeric","numeric","NULL","character","character","NULL","character","character","character","character","numeric","numeric","character","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL")

lightcurve <- read.csv(file=llightcurve_name,header=TRUE,check.names=TRUE,na.strings="NA")

# lightcurve, which is the AAVSO data read in via read.csv with header=TRUE
totalRec = length(lightcurve$JD)
cat("\n\nTotal records read from file: ",totalRec,"\n\n")

# make a list of every record with a Magnitude reported
goodMags <- !is.na(lightcurve$Magnitude)

# replace missing airmass for JM
cat("replacing missing airmass...\n")
for (index in 1:length(missingAirmass)){
        jmtest = lightcurve$Observer_Code == missingAirmass[index] & is.na(lightcurve$Airmass)
        lightcurve$Airmass[jmtest] <- AirMass(lightcurve$JD[jmtest],missingAMLocs[index,], tabbysLoc)
}

# loop over the desired bands
numBands = length(allBands$bandinQ)
cleanBand <- matrix(nrow=totalRec,ncol=numBands)
allFits <- list()
allQFits <- list()
allResistFit <- list()
tmin <- head(lightcurve$JD,n=1)
cat("earliest time in file: ",tmin,"\n")
index = 1

# edit specific observers over time range(s) and specific passband
cat("editing specific observers\n")
editCurve <- ObserverJDEdit(editUser,lightcurve)

#clean and separate and the bands in question

for (thisBand in allBands$bandinQ) {
	cat("cleaning ",thisBand,"\n")
	# clean the data for this passband
	cleanBand[,index] <- cleanAAVSO3(lightcurve,thisBand,ExclCodes,includeExclude,maxairmass,maxuncertainty,wildsd,earliestJD,okComparison) & editCurve
	index = index + 1
}
	
	#bin the data for each Observer
cat("binning the data\n")
binCurve <- binAAVSO(lightcurve,cleanBand,allBands,deltaJD)

# test for less than max uncertainty
uncertaintyTest <- binCurve$Uncertainty <= maxBinUncertainty

index = 1

# loop over the passbands and do regression for each one
cat("fitting the data \n")
for (thisBand in allBands$bandinQ) {
	# fold in test for passband
	btest <- (binCurve$Band == thisBand) & uncertaintyTest
	desmat <- binCurve$JD[btest] - tmin # subtract off the earliest time to avoid numerical problems
	if(length(desmat) < 2) {
			cat("\n\nWarning: fewer than 2 points in the band:",thisBand,"\n")
			next
		}
	
	# determine the bin weights, if any
	if (userlm & !plotMARS) {
		# use robust algorithm
		thisFit <- rlm(binCurve$Magnitude[btest] ~ desmat,na.action="na.omit",psi=psi.bisquare)
	} else if(!plotMARS){
		if (weightedBins) {
			binWeights <- 1/binCurve[btest,"Uncertainty"]
		} else {
			binWeights <- NULL
		}
		# do the linear regression in the binned data for this band, starting at the earliest time in the file
		thisFit <- lm(binCurve[btest,"Magnitude"] ~ desmat, weights= binWeights)
	}	
	
	# if MARS is selected (plotMARS == TRUE), use it to do the regression
	if(plotMARS) {
		# use MARS algorithm in earth()
		cat("\n Using MARS algorithm\n")
		if(splineRaw) {
			# do the spline fit on the unbinned data
			desmat <- lightcurve$JD[cleanBand[,index]] - tmin
			thisFit <- earth(x=desmat,y=lightcurve$Magnitude[cleanBand[,index]],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)
		} else {
			thisFit <- earth(x=desmat,y=binCurve$Magnitude[btest],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)
	#		print(class(thisFit))
		}
	} 
	
	# output a summary to the console 	
	cat("\n\n Band",thisBand,"summary")
	print(summary(thisFit))
	
	# try resistant regression if selected
	if (tryLQS) {
		resistFit <- lqs(formula = binCurve$Magnitude[btest] ~ desmat)
		cat("\n\n Band",thisBand," resistant fit coefficients")
		cat(resistFit$coefficients,"\n")
		#build up the matrix of resistant fits
		allResistFit <- rbind(allResistFit,resistFit)
	}

	if (plotQuadratic) {
		# fit a quadratic just for fun
		desmat <- outer(binCurve$JD[btest]-tmin,1:2,"^")
		qfit <- lm(binCurve$Magnitude[btest] ~ desmat, weights= binWeights)
		allQFits <- rbind(allQFits,qfit)
	}
	#package the fits into one matrix
	allFits <- rbind(allFits,thisFit)

	if(index==1) {
		allClean <- cleanBand[,index]
	} else {
		allClean <- cleanBand[,index] | allClean
	}
	index = index + 1	
}

##################################### Plot This Stuff ##################################
# axis limits

myxlims = c(startPlot,max(binCurve$JD,na.rm=TRUE))

# calculate pretty y limits

myYlims = c(ceiling(max(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE)*10)/10,floor(10*min(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE))/10) # set up Y limits for reversed Y axis

# set up plot title text
howManyObs = length(unique(binCurve$Observer_Code[uncertaintyTest]))
#ocodesInTitle <- paste("Observer",head(ExclCodes,n=3),inclWord,sep=" ")
#if (length(ExclCodes) > 3){ocodesInTitle <- c(ocodesInTitle,paste("and",howManyObs - 3,"more observer code(s)",sep=" "))}

myBands = paste(allBands$bandinQ,collapse=" ")
titleString <- c(paste("AAVSO",myBands,"Data with",deltaJD,"Day Bins",sep=" "), paste(as.character(howManyObs),"Observers",sep=" "))
myPlotTitle <- paste(titleString,collapse="\n")


# plot the cleaned and binned data, the fit lines and the excluded points
icol=1
if(plotRelTimes) {
	myTimes <- binCurve$JD - tmin
	jdLine <- jdLine - tmin
	myxlims <- myxlims - tmin
	lcTimes = lightcurve$JD - tmin
	myXLabel <- paste("Julian Date -",as.character(tmin),sep=" ")
} else {
	myTimes <- binCurve$JD
	myXLabel <- "Julian Date"
	}

for (thisBand in allBands$bandinQ) {
#	if (icol > 1){par(new=TRUE)}
	ourCleanData <- cleanBand[,icol]
	btest <- (binCurve$Band == thisBand) & uncertaintyTest
		
	if(icol==1) {
		plot(myTimes[btest],binCurve[btest,"Magnitude"],col=allBands$plotColor[icol],xlab= myXLabel,ylab="Magnitude",xlim= myxlims,ylim = myYlims,main=myPlotTitle,pch=3,cex.main=0.7)
	} else {
		points(myTimes[btest],binCurve[btest,"Magnitude"],col=allBands$plotColor[icol],pch=3)
	}
	
	# plot line fit
	if(plot2Lines) {
		cat(red("plt2Lines doesn't really work. Use earth() instead"))
	} else if(!plotMARS) {
		# plot the lm() or rlm() fit
		myslope <- coefficients(allFits[icol,])[2]
		basemag <- coefficients(allFits[icol,])[1]
		curve(basemag + myslope*(x-tmin),from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col="black")
	} 
	# plot the MARS fit if that is selected
	if (plotMARS) {
		mars <- allFits[icol,]
		if (splineRaw) {
			lines(x=lcTimes[cleanBand[,icol]],y=mars$fitted.values,col= "black",lwd=2)
		} else {
			lines(x=myTimes[btest],y=mars$fitted.values,col= "black",lwd=2)
		}
	}
	#optionally plot the LQS fit
	if (tryLQS) {
		myslope <- coefficients(allResistFit[icol,])[2]
		basemag <- coefficients(allResistFit[icol,])[1]
		curve(basemag + myslope*(x-tmin),from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col=lqsColor)
	}
	
	# an option to plot the quadratic fit
	if (plotQuadratic) {
		#quadratic fit
		basemag <- coefficients(allQFits[icol,])[1]
		linterm <- coefficients(allQFits[icol,])[2]
		qterm <- coefficients(allQFits[icol,])[3]
		curve(basemag + linterm*(x-tmin) + qterm*(x -tmin)^2,from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col="red")
	}
	
#plot the excluded points, if option selected
	if(plotExcluded) {
		# plot excluded points in black
		IResid <- !ourCleanData & goodMags & (lightcurve$Band == allBands$bandinQ[icol])		
		points(lcTimes[IResid],lightcurve[IResid,"Magnitude"],col="black",pch=20) # plot removed points in black
	}	
	icol = icol + 1
}

if (!is.na(plotMee)) {
	imSpecial <- grep(plotMee,binCurve$Observer_Code,ignore.case=TRUE)
	points(myTimes[imSpecial],binCurve$Magnitude[imSpecial],col=meeColor,pch=20,cex=1.5)
}

grid(col="black")

##### draw vertical lines at various dates
if(drawDateLine) { verticalDateLines(jdLine, jdLineText, myYlims, jdLineColor)}



# generate a time series and plot if called for

if (generateTS) {
	quartz("Time Series")
	tsMain <- "Time Series"	
	myts <- binAAVSO_ts(lightcurve, cleanBand,allBands, tsBinWidth)
	if (smoothTS) {
		tsMain <- paste("Smoothed",tsMain," - Order = ", tsSmoothOrder)
		tsMain 
		okts <- !sapply(X=myts,FUN=is.na,simplify="logical")
		for(iband in 1:numBands) {
			tsScratch <- myts[,iband]
			tsScratch[okts[,iband]] <- sma(data = myts[okts[,iband],iband],order=tsSmoothOrder)$fitted
			if (iband == 1) {
					bts <- tsScratch
				} else {
					bts  <- cbind(bts,tsScratch,deparse.level=0)
				}
		}
	} else {
		bts <- myts[,1:numBands] # not smoothed
	}
	if(numBands > 1 ) {
		colnames(bts) <- sapply(allBands$bandinQ,FUN=paste,"Mag",sep=" ")
		plot.ts(bts,plot.type="multiple",xlim= myxlims,ylim = myYlims,main=tsMain,cex.axis=0.9,lwd=2)
	} else {
		plot(bts,ylim=myYlims,main=tsMain,lwd=2,col=allBands$plotColor,ylab=paste(allBands$bandinQ,"Mag"))
		points(myts,col="grey",pch=20,cex=0.5)

	}
	grid(col="black")
}

# plot the residuals if desired:
if (plotResiduals & !splineRaw) {
	quartz("Flux Relative to Fit")
	irow <- 1
	for (thisBand in allBands$bandinQ) {
		btest <- (binCurve$Band == thisBand) & uncertaintyTest
		relFlux <- sapply(allFits[irow,]$residual,ReverseMagnitude)
		fluxYLims <-c(min(sapply(allFits[irow,]$residual,ReverseMagnitude))-0.01,
						max(sapply(allFits[irow,]$residual,ReverseMagnitude))+0.01)
		if(irow == 1) {
			plot(myTimes[btest],relFlux,col= allBands$plotColor[irow],xlab=myXLabel,ylab="Relative Flux",xlim= myxlims,
				ylim= fluxYLims,main="Residuals",pch=20,cex.main=1.0)
		} else {
			points(myTimes[btest],relFLux,col=allBands$plotColor[irow],pch=20)
		}
		irow <- irow + 1
	}
	grid(col="black")
}

if(drawDateLine) { verticalDateLines(jdLine, jdLineText, fluxYLims, jdLineColor)}



# who the observers were for cleaned data:
buncha <- AAVSOObsStats(lightcurve,cleanBand,allBands)
cat("\n\nObserver Summary - Raw Observations\n")
print(buncha)

probOK <- logical()

for (thisBand in allBands$bandinQ) {probOK <- cbind(probOK,!is.na(binCurve$JD) & uncertaintyTest)}
binda <-AAVSOObsStats(binCurve,probOK,allBands)
cat("\n\n    Observer Summary - Binned Observations with acceptable scatter\n")
print(binda)


# summaries of the fits 

cat("\n\n    Observer Code: ", ExclCodes," observations",inclWord,"\n")
cat("   ",length(lightcurve$JD),"total observations loaded\n")
#cat("    ",length(lightcurve[cleanI | cleanR | cleanV | cleanB,"JD"]),"raw observations after cleaning\n")
cat("    ",length(binCurve$JD[uncertaintyTest]),"binned observations with",deltaJD,"day bins")

#cat("\n\nObserver Stats for this set of fits")
#AAVSOObsStats(lightcurve,cleanBand,allBands)

for (index in 1:numBands) {
	cat("\n\n",allBands$bandinQ[index],"    Summary\n")
	pctPerYear <- (10^(-1*coefficients(allFits[index,])[2]*365.24/2.5)-1)*100
	cat("    ",coefficients(allFits[index,])[2]*365.24*100,"magnitudes per century","or",pctPerYear,"% per year\n")
}


