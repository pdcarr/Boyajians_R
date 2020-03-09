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
rm(deriv.mat)
rm(resid.mat)
rm(smoove.fit)
#rm(bts)
# options
options(digits=12) # hard to read JDs without this setting
# load the required packages
library("MASS") # for rlm() and lqs()
library("smooth") # for smoothing
library("earth") # for MARS
library("crayon") # to add color to text
library("Hmisc") # for error bars

allFits <- list()


##################

# all the inputs are in the sourced file
source("input_files/aavso_bin_input_parameters.R")
# source("input_files/Betelgeuse_bin_input_parameters.R")
source("input_files/observer_edits.R")
source("input_files/missing_airmass.R")
source("input_files/VlineParams.R")
source("input_files/dip_mask.R")
source("input_files/spline_control.R")

if (includeExclude) {
	inclWord <- "used"
} else {
	inclWord <- "not used"
}


# lightcurve, which is the AAVSO data read in via read.csv with header=TRUE
lightcurve <- read.csv(file=llightcurve_name,header=TRUE,check.names=TRUE,na.strings="NA",stringsAsFactors=FALSE)
totalRec = length(lightcurve$JD)
cat("\n\nTotal records read from file: ",totalRec,"\n\n")

if ((merge.asassn) & (sum(allBands$bandinQ %in% asassn.bands) > 0)){
	cat("\nmerging in ASASSN Data\n")
  cat("\nASAS-SN cameras",asassn.cameras,"\n")
	asassn_data <- read.csv(asassn.csv.file,header=TRUE,stringsAsFactors=FALSE)
	if(!exists("asassn_data")) {
	  cat("\nfailed to open ASASSN_FILE\n")
	  merge.asassn <- FALSE
	}
	if(sum(asassn.cameras %in% asassn_data$Camera) == 0)
	{
	  cat("\nNone of the required ASAS-SN cameras in file\n")
	}
	ok.camera <- !is.na(match(asassn_data$Camera,asassn.cameras))
	lightcurve <- asassn.merge(lightcurve,asassn_data[ok.camera,],asassn.code,g.to.V= convert.asassn,star.BminusV= our.BminusV,V.bias= converted.V.bias) # function merges asassn data into lightcurve
} else {
  cat("\nNot merging ASAS-SN data - no relevant bands\n")
}


latestJD <- max(lightcurve$JD,na.rm=TRUE)
# make a list of every record with a Magnitude reported
goodMags <- !is.na(lightcurve$Magnitude)

# replace missing airmass for observer codes in missing Airmass
cat("replacing missing airmass...\n")
for (whos.missing in missingAirmass){
        in.codes <- whos.missing %in% ExclCodes
        jmtest <- lightcurve$Observer_Code == whos.missing & is.na(lightcurve$Airmass) & 
          ((includeExclude & in.codes) | (!includeExclude & !in.codes))
        if(sum(jmtest) > 0) {
          which.loc <- match(whos.missing,missingAM.data$lax.observers)
          this.loc <- c(missingAM.data$missing.lats[which.loc],missingAM.data$missing.lons[which.loc])
          lightcurve$Airmass[jmtest] <- AirMass(lightcurve$JD[jmtest],this.loc,tabbysLoc)
        }
}

# loop over the desired bands
cat("initializing loop over bands\n")
numBands = length(allBands$bandinQ)
cleanBand <- matrix(nrow=nrow(lightcurve),ncol=numBands)
allFits <- list()
all.smoove <- list()
allQFits <- list()
allResistFit <- list()
tmin <- head(lightcurve$JD,n=1)
cat("earliest time in file: ",tmin,"\n")

#set tmin to some rounder number:
if(exists("pretty.JD.interval")) {
	tmin <- pretty.JD.interval*floor(tmin/pretty.JD.interval)
} else {tmin <- 10*floor(tmin/10)}

# edit specific observers over time range(s) and specific passband
cat("editing specific observers\n")
editCurve <- ObserverJDEdit(editUser,lightcurve)

##################### clean and separate and the bands in question
index = 1
for (thisBand in allBands$bandinQ) {
	cat("cleaning ",thisBand,"\n")
	# clean the data for this passband
	cleanBand[,index] <- cleanAAVSO3(lightcurve,thisBand,ExclCodes,includeExclude,maxairmass,maxuncertainty,wildsd,earliestJD,okComparison,bad.comp.star) &
						 editCurve
	##### NOTE: this is really slow, so only do it if it's important to get the ensemble means for each bin
	if(plot.ensemble & use.static.biases) {
		cat("\nApplying static biases. This could take some time..","\n")
		lightcurve[cleanBand[,index],] <- lightcurve.biases(lightcurve[cleanBand[,index],],biasObserver)
	}
	index = index + 1
}
	
###################### bin the data for each Observer
cat("binning the data\n")
bin.list <- binAAVSO(lightcurve,cleanBand,allBands,deltaJD,weightless,trial.bin,min.population)
binCurve <- as.data.frame(bin.list[1])
ensemble.curve <- as.data.frame(bin.list[2])
#binCurve <- binAAVSO(lightcurve,cleanBand,allBands,deltaJD) # old function call


######### filter bins for excess uncertainty
# test for less than max uncertainty
uncertaintyTest <- binCurve$Uncertainty <= maxBinUncertainty & binCurve$Uncertainty > 0

index = 1
used.in.fit <- matrix(nrow=length(binCurve$JD),ncol=length(allBands$bandinQ))
# loop over the passbands and do regression for each one

cat("fitting the data \n")
bias.vec <- vector(mode="numeric",length=length(binCurve$JD))

for (thisBand in allBands$bandinQ) {
	# fold in test for passband
	btest <- (binCurve$Band == thisBand) & uncertaintyTest
	bias.vec[!btest] <- NaN
	
    # mask out data taken during dips (0 weight)
    dipless <- filter.dips.JD(binCurve$JD[btest],dip.mask)
    dipless <- dipless | !mask.Dips #option to turn off the dip masking
    
	desmat <- binCurve$JD[btest] - tmin # subtract off the earliest time to avoid numerical problems
	if(length(desmat) < 2) {
			cat("\n\nWarning: fewer than 2 points in the band:",thisBand,"\n")
			next
		}
		
	# apply any defined observer biases to the binned magnitude(s) if(use.static.biases)
	if (exists("biasObserver") & exists("use.static.biases")) {
		if(use.static.biases & !plot.ensemble) {
			biasBand <- biasObserver$band == thisBand
			# if there are observer biases for this band, then apply them
			# loop over all the observers in the biases data frame
			for (thisObs in unique(biasObserver$obsCode[biasBand])) {
				cat("\n applying biases for",thisObs,"in",thisBand)
	#			cat("\nApplying biases to",thisBand,"\n")			
				myObs <- binCurve$Observer_Code[btest] == thisObs
				myBias <- as.numeric(biasObserver$bias[biasBand & biasObserver$obsCode==thisObs][1])
				if(sum(myObs) > 0) {
					 mag.x <- binCurve$Magnitude[btest][myObs]
					 mag.x <-  mag.x - myBias
					 binCurve$Magnitude[btest][myObs] <- mag.x
					 bias.vec[btest][myObs] <- myBias
				}
			}
		}		
	}
	
	# determine the bin weights, if any then do the fit(s)
	used.in.fit[,index] <- btest # create the column using btest (test for correct band)
	if (userlm & !plotMARS & !perform.smooth) {
		# use robust algorithm
		thisFit <- rlm(binCurve$Magnitude[btest] ~ desmat,na.action="na.omit",psi=psi.bisquare)
	} else {
		if (weightedBins) {
			de.weight <- dipless
			if(exists("weightless") & !is.na(weightless)) {
			    for (thisObs in weightless) {
				    	# there may be observers whose observations we want to weight zero
			    		de.weight <- de.weight & !(binCurve$Observer_Code[btest] == thisObs)
			    }
		    }
	    	
		    used.in.fit[btest,index] <- de.weight	#  set weights of dips to zero

			binWeights <- (1/binCurve[btest,"Uncertainty"])*as.numeric(de.weight)
		} else {
			binWeights <- NULL
		}
		# do the linear regression in the binned data for this band, starting at the earliest time in the file
#		thisFit <- lm(binCurve[btest,"Magnitude"] ~ desmat, weights= binWeights)
	}	

    if(index==1) {
	    	binCurve <- cbind(binCurve,used.in.fit[,index],deparse.level = 1)
		} else {
		    	binCurve$used.in.fit <- binCurve$used.in.fit | as.logical(used.in.fit[,index])
	    }

	# if MARS is selected (plotMARS == TRUE), use it to do the regression
	if(plotMARS) {
		# use MARS algorithm in earth()
		cat("\n Using MARS algorithm\n")
		if(splineRaw) {
			# do the spline fit on the unbinned data
			desmat <- lightcurve$JD[cleanBand[,index]] - tmin
			# hello earth()
			thisFit <- earth(x=desmat,y=lightcurve$Magnitude[cleanBand[,index]],
							nk= marsOrder,
							pmethod= marsPMethod,
							penalty = marsPenalty,
							thresh = mars.thresh,
							minspan = mars.minspan)
		} else {
			thisFit <- earth(x=desmat,y=binCurve$Magnitude[btest],
								nk= marsOrder,
								pmethod= marsPMethod,
								penalty = marsPenalty,
								weights=binWeights,
								thresh = mars.thresh,
								minspan = mars.minspan)
		}
		allFits <- rbind(allFits,thisFit)
	} 
#	print(length(desmat))
	
	# do a smooth.spline on the bins if called for
	if(perform.smooth) {
		smoove.fit <- smooth.spline(desmat,binCurve$Magnitude[btest],
									w=binWeights,
									all.knots=FALSE,nknots= smooth.n.knots,
									keep.data=TRUE,cv=FALSE,penalty= df.penalty)
		cat("\n\n Smooth Spline Fit: \n")
		print(smoove.fit$call)
		cat("\n")
		print(smoove.fit$fit)
		all.smoove[[index]] <- smoove.fit
	}
	
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
	#package the cleanBands into one matrix
	if(index==1) {
		allClean <- cleanBand[,index]
	} else {
		allClean <- cleanBand[,index] | allClean
	}
	index = index + 1	
}
binCurve <- cbind(binCurve,bias.vec,deparse.level=1)

##################################### Plot This Stuff ##################################
# axis limits

if(exists("stop.plot")) {
	if(!is.na(stop.plot) & stop.plot > startPlot) {
		myxlims <- c(startPlot,stop.plot)
	} else {
		myxlims <- c(startPlot,max(binCurve$JD,na.rm=TRUE))
	}
} else {
	myxlims <- c(startPlot,max(binCurve$JD,na.rm=TRUE))
}

# calculate pretty y limits

if(exists("pretty.interval")) {
	myYlims <- c(ceiling(max(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE)* pretty.interval)/pretty.interval,
				floor(pretty.interval*min(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE))/pretty.interval) # set up Y limits for reversed Y axis
} else {
	myYlims = c(ceiling(max(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE)*10)/10,
				floor(10*min(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE))/10) # set up Y limits for reversed Y axis	
}


# set up plot title text
howManyObs = length(unique(binCurve$Observer_Code[uncertaintyTest]))
#ocodesInTitle <- paste("Observer",head(ExclCodes,n=3),inclWord,sep=" ")
#if (length(ExclCodes) > 3){ocodesInTitle <- c(ocodesInTitle,paste("and",howManyObs - 3,"more observer code(s)",sep=" "))}

myBands = paste(allBands$bandinQ,collapse=" ")
titleString <- c(paste(myBands,"Data with",length(desmat),"bins",sep=" "), paste(as.character(howManyObs),"Observers",sep=" "))
myPlotTitle <- paste(titleString,collapse="\n")

# initialize some data structures
these.resids <- vector(mode="numeric",length=length(binCurve$JD))
my.derivs <- vector(mode="numeric",length=length(binCurve$JD))
resid.mat = matrix(ncol=length(binCurve$JD))
deriv.mat = matrix(ncol=length(binCurve$JD))

# plot the cleaned and binned data, the fit lines and the excluded points
icol=1
if(plotRelTimes) {
	myTimes <- binCurve$JD - tmin
	jdLine <- jdLine - tmin
	if(add.predict > 0) {myxlims[2] <- myxlims[2] + add.predict}
	myxlims <- myxlims - tmin
	lcTimes = lightcurve$JD - tmin
	myXLabel <- paste("Julian Date -",as.character(tmin),sep=" ")
} else {
	myTimes <- binCurve$JD
	myXLabel <- "Julian Date"
}

deriv.bounds <- c(NA,NA)
quartz("AAVSO Magnitude Data")
bin.predict <- vector(mode="numeric",length=length(binCurve$JD))
bin.predict[1:length(binCurve$JD)] <- NaN

for (thisBand in allBands$bandinQ) {
	ourCleanData <- cleanBand[,icol]
	btest <- (binCurve$Band == thisBand) & uncertaintyTest
	my.y.plus <- binCurve$Magnitude[btest] + binCurve$Uncertainty[btest]
	my.y.minus<-  binCurve$Magnitude[btest] - binCurve$Uncertainty[btest]
		
	if(icol==1) {
		errbar(myTimes[btest],binCurve[btest,"Magnitude"],
				yplus=my.y.plus,yminus=my.y.minus,
				col=allBands$plotColor[icol],
				xlab= myXLabel,ylab="Magnitude",xlim= myxlims,ylim = myYlims,
                main=myPlotTitle,pch=3,cex.main=0.7,
				add=FALSE,errbar.col=ebar.color)
		points(myTimes[btest],binCurve[btest,"Magnitude"],col=allBands$plotColor[icol],pch=3)
		title(main=myPlotTitle)
	} else {
		errbar(myTimes[btest],binCurve$Magnitude[btest],yplus=my.y.plus,yminus=my.y.minus,
				col=allBands$plotColor[icol],errbar.col=allBands$plotColor[icol],pch=3,add=TRUE)
	}
	
	# plot line fit
	 if(!plotMARS & userlm) {
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
	

	# plot a smooth.spline on the bins if called for
	if(perform.smooth) {
		smoove.fit <- all.smoove[[icol]]
		these.values <- predict(smoove.fit,myTimes[btest])$y
		# bin.predict[btest] <- these.values

		lines(myTimes[btest],these.values,col=smoove.color,lwd=3) # plot as a line of specified color
		##### add prediction if add.predict > 0
		if(add.predict > 0) {
			predict.days <- seq(from=tail(myTimes[btest],n=1),to=tail(myTimes[btest],n=1)+add.predict,by=1)
			add.prediction <- predict(smoove.fit,predict.days,deriv=0)
			lines(predict.days,add.prediction$y,lty="dashed",lwd=3,col=smoove.color)
		}
		
		these.resids[btest] <- binCurve$Magnitude[btest] - these.values # store residuals
		if(icol==1) {
			resid.mat[icol,] <- these.resids
		} else {
			resid.mat <- rbind(resid.mat,these.resids)
		}
		if(smooth.deriv) { 
				my.derivs[btest]  <- predict(smoove.fit, myTimes[btest],deriv=1)$y
				deriv.bounds[1] <- min(deriv.bounds[1],my.derivs,na.rm=TRUE)
				deriv.bounds[2] <- max(deriv.bounds[2],my.derivs,na.rm=TRUE)
				if (icol==1) {
					deriv.mat[icol,btest] <- my.derivs[btest]
				} else {
					deriv.mat <- rbind(deriv.mat,my.derivs)
				}
		}
		
	}
	
	#optionally plot the LQS fit
	if (tryLQS) {
		myslope <- coefficients(allResistFit[icol,])[2]
		basemag <- coefficients(allResistFit[icol,])[1]
		curve(basemag + myslope*(x-tmin),from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col=lqsColor)
	}
	icol = icol + 1
}	

if (!is.na(plotMee)) {
	imSpecial <- grep(plotMee,binCurve$Observer_Code,ignore.case=TRUE)
	points(myTimes[imSpecial],binCurve$Magnitude[imSpecial],col=meeColor,pch=20,cex=1.5)
}
# tack the predictions onto the binCurve
binCurve <- cbind(binCurve,bin.predict,deparse.level=1)
grid(col="black")

##### draw vertical lines at various dates
if(drawDateLine) { verticalDateLines(jdLine, jdLineText, myYlims, jdLineColor)}


################ plot the residuals if desired:
if (plotMARS & plotResiduals & !splineRaw) {
	quartz("Flux Relative to Fit")
	irow <- 1
	for (thisBand in allBands$bandinQ) {
		btest <- (binCurve$Band == thisBand) & uncertaintyTest
		not.in.dip <- used.in.fit[,irow][btest]
		relFlux <- sapply(allFits[irow,]$residual,ReverseMagnitude)
		fluxYLims <-c(min(sapply(allFits[irow,]$residual,ReverseMagnitude))-0.01,
						max(sapply(allFits[irow,]$residual,ReverseMagnitude))+0.01)
		if(irow == 1) {
			plot(myTimes[btest][not.in.dip],relFlux[not.in.dip],col= allBands$plotColor[irow],xlab=myXLabel,ylab="Relative Flux",xlim= myxlims,
				ylim= fluxYLims,main="Residuals",pch=20,cex.main=1.0,type= res.plot.type)
			points(myTimes[btest][!not.in.dip],relFlux[!not.in.dip],col="grey",pch=20)
		} else {
			points(x=myTimes[btest][not.in.dip],y= relFlux[not.in.dip],col=allBands$plotColor[irow],pch=20)
			points(myTimes[btest][!not.in.dip],relFlux[!not.in.dip],col="grey",pch=20)

		}
		irow <- irow + 1
	}
	##### draw vertical lines at various dates
	if(drawDateLine) { verticalDateLines(jdLine, jdLineText, fluxYLims, jdLineColor)}
	grid(col="black")
} else if(!plotMARS & plotResiduals & perform.smooth) {
	quartz("Flux Relative to Fit")	
	irow <- 1
	# loop over all the bands under consideration
	for (thisBand in allBands$bandinQ) {
		btest <- (binCurve$Band == thisBand) & uncertaintyTest
		not.in.dip <- used.in.fit[,irow][btest]
		relFlux <- sapply(resid.mat[irow,btest],ReverseMagnitude) # convert magnitiude difference to flux ratio
		fluxYLims <-c(min(relFlux)-0.01,max(relFlux)+0.01)	# calcuate y plot limits with a little more room
		if(irow == 1) {
			plot(myTimes[btest][not.in.dip],
                relFlux[not.in.dip],
                col= allBands$plotColor[irow],
                xlab=myXLabel,ylab="Relative Flux",
                xlim=myxlims,
                ylim= fluxYLims,
                main="Residuals",
                pch=20,cex.main=1.0,
                type= res.plot.type)
            
			points(myTimes[btest][!not.in.dip],
                relFlux[!not.in.dip],
                col="grey",pch=20)
		} else {
			points(x=myTimes[btest][not.in.dip],y= relFlux[not.in.dip],col=allBands$plotColor[irow],pch=20)
			points(myTimes[btest][!not.in.dip],relFlux[!not.in.dip],col="grey",pch=20)

		}
		irow <- irow + 1
		grid(col="black")
		##### draw vertical lines at various dates
		if(drawDateLine) { verticalDateLines(jdLine, jdLineText, fluxYLims, jdLineColor)}
	}

}

# smooth fit derivative plot

if(perform.smooth & smooth.deriv) {
			quartz("Derivative of Spline Fit")
			
			irow <- 1
			for (thisBand in allBands$bandinQ) {
				btest <- (binCurve$Band == thisBand) & uncertaintyTest	# the subset of the bins to use
				myXlabel <- "time"
				myYlabel <- "1st Derivative of Magnitude wrt Time"
				if(irow == 1) {
					plot(myTimes[btest],deriv.mat[irow,btest],col= allBands$plotColor[irow],
						xlab=myXLabel,ylab=myYlabel,
						main="Derivative of Smooth Spline",
						cex.main=1.0,type= "l",lwd=2,
						xlim= myxlims,ylim=deriv.margin*deriv.bounds)
				} else {
					lines(myTimes[btest],deriv.mat[irow,btest],col= allBands$plotColor[irow],lwd=2)
				}
				if(add.predict > 0) {
					add.times <- seq(from=tail(myTimes[btest],n=1),to=tail(myTimes[btest],n=1) + add.predict,by=1)
					add.derivs <-predict(all.smoove[[irow]],add.times,deriv=1)
					lines(add.times,add.derivs$y,col=allBands$plotColor[irow],lty="dashed",lwd=2)
				}
				irow <- irow + 1
			}
			grid(col="black")
			##### draw vertical lines at various dates
			if(drawDateLine) { verticalDateLines(jdLine, jdLineText, myYlims, jdLineColor)}

		}

if(drawDateLine) { verticalDateLines(jdLine, jdLineText, deriv.bounds, jdLineColor)}



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

########## plot ensemble averages

if(plot.ensemble) {
	ensemble.title = "Unweighted Ensemble averages"
    titleString <- c(ensemble.title, paste(as.character(howManyObs),"Observers",sep=" "))
    myPlotTitle <- paste(titleString,collapse="\n")
	plot.times <- ensemble.curve$JD - tmin # same time offset as the other plots
	x.string <- paste("Julian Date - ",tmin) # plot x axis label
	quartz("ensemble average bins")
	i.band <- 1
	for (thisBand in allBands$bandinQ) {
		# first band in list
		band.test <- ensemble.curve$Band == thisBand
		if(i.band==1) {	
			plot(x=plot.times[band.test],y=ensemble.curve$Magnitude[band.test],
			xlim= myxlims,ylim = myYlims,
			xlab=x.string,ylab="Magnitude",
			col=allBands$plotColor[i.band],
			pch=20,main=myPlotTitle)
		} else { 		# subsequent bands in the list, if any
			points(x=plot.times[band.test],y=ensemble.curve$Magnitude[band.test],
			col=allBands$plotColor[i.band],
			pch=20)
		}
	
		# plot a smooth.spline on the bins if called for
		if(perform.smooth) {
#			cat("\nstarting plot of smooth fit\n")
			smoove.fit <- all.smoove[[i.band]]
			these.values <- predict(smoove.fit,plot.times)$y
			bin.predict[btest] <- these.values
	
			lines(plot.times,these.values,col=allBands$plotColor[i.band],lwd=3) # plot as a line of specified color
			
		}
		i.band <- i.band + 1	# point to next band
	}
	grid(col="black")
	if(drawDateLine) { verticalDateLines(jdLine, jdLineText, myYlims, jdLineColor)}
}

