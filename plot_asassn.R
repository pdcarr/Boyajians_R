# read in and plot ASASSN data
# a little housekeeping
options(digits=12)
library("earth") # for MARS
library("crayon") # to add color to text

source("plot_funcs.R")
source("data_funcs.R")

rm("superObs")
rm("allSuperObs")


# the name of the file
asassn.csv.file <- "data/ASASSN_latest.csv"
# cull limits to reject really wild points
mag.cull.limit = 0.05
# mag margin for expanding plot limits (positive number)
mag.margin <- 0.025
tmargin <- 1 # positive number, days extra on xscale. 
##################################### binning control
bin.width = 1 # days
bin.it <- TRUE

######### plot controls
plot.col = "darkgreen"
plot.pch <- 20
t.epsilon = 1.0 # days
mag.epsilon <- 0.025 # magnitudes



####### MARS
mars.fit <- FALSE
marsOrder <- 21 # maximum number of knots
marsPenalty <- 5 # set to 0 to avoid penalizing knots in pruning pass
mars.thresh <- 0.000002
#marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "exhaustive" # set to "none" to avoid pruning
min.span <- 2

######## Smooth Spline
asassn.n.knots <- 7

##########
earliestJD <- 2457449
earliestJD <- 2457880
#earliestJD <- 2457990
source("input_files/VlineParams.R")
source("input_files/dip_mask.R")


#####################################################################
# load the files
asassn_data <- read.csv(asassn.csv.file,header=TRUE)
# clean obvious wild points
mag.median <- median(asassn_data$mag)
good_data <- (asassn_data$mag < (mag.median + abs(mag.cull.limit))) & (asassn_data$mag > (mag.median- abs(mag.cull.limit)))

#fit a spline after filtering out the dips
# which points are NOT in the known dips?
not.dips <- filter.dips.JD(asassn_data$HJD[good_data],dip.mask)
# set a nice round tmin
tmin <- max(min(asassn_data$HJD[good_data]),earliestJD)
tmin <- 20*floor(tmin/20)
tmax <- max(asassn_data$HJD[good_data])

x.limits <- c(0,(tmax-tmin) + t.epsilon)
y.limits <- c(max(asassn_data$mag[good_data],na.rm=TRUE) + mag.epsilon,
			 min(asassn_data$mag[good_data],na.rm=TRUE) - mag.epsilon)
x.label <- paste("JD - ",tmin)


# if binning, do this:
if (bin.it) {
	plot.title <- paste("ASAS-SN Data binned with width",bin.width,"Days")
	y.label <- paste("Mean magnitude over bin")
	allSuperObs <- data.frame(HJD=numeric(),V.mag=numeric(),Uncertainty=numeric(),stringsAsFactors=FALSE)
	superObs <- data.frame(HJD=numeric(),V.mag=numeric(),Uncertainty=numeric(),stringsAsFactors=FALSE)
	# bin the data
	for (t in seq(tmin,tmax,by=bin.width)) {
		in.bin <- !is.na(asassn_data$HJD[good_data]) & (asassn_data$HJD[good_data] >= t) & (asassn_data$HJD[good_data] < (t + bin.width))
		n.in.bin <- length(asassn_data$HJD[good_data][in.bin])
#		print(n.in.bin)
		if (n.in.bin > 0) {
			superObs[1,"HJD"] <- mean(asassn_data$HJD[good_data][in.bin],na.rm=TRUE)
			superObs[1,"V.mag"] <- mean(asassn_data$mag[good_data][in.bin],na.rm=TRUE)
			superObs[1,"Uncertainty"] <- mean(asassn_data$mag_err[good_data][in.bin],na.rm=TRUE)/sqrt(n.in.bin)
			allSuperObs <- rbind(allSuperObs,superObs)
		} else {
			next
		}
	}
	# mask dips from the fit
    dipless <- filter.dips.JD(allSuperObs$HJD,dip.mask)
    binWeights <- as.numeric(dipless) # weight of 1 if not in a known dip, 0 otherwise.
	# calculate robust linear fit
	desmat <- allSuperObs$HJD - tmin
	if(mars.fit) {
	    theFit <- earth(x=desmat,
	                    y=allSuperObs$V.mag,
	                    weights=binWeights,
	                    pmethod="exhaustive",
	                    nk= marsOrder,
	                    penalty= marsPenalty,
	                    minspan=min.span)
	
		print(summary(theFit))
		fitpoints <- theFit$fitted.values 
	} else {
		smoove.fit <- smooth.spline(x=desmat,
									y= allSuperObs$V.mag,
									w= binWeights,
									all.knots=FALSE,nknots= asassn.n.knots,
									keep.data=TRUE,
									cv=TRUE,
									penalty= 1)
									
		cat("\n\n Smooth Spline Fit: \n")
		print(smoove.fit$call)
		cat("\n")
		print(smoove.fit$fit)
		fit.points <- predict(smoove.fit,x=desmat)$y
		
	}
	quartz("binned Data Plot")
	plot(x=desmat[dipless],
		y=allSuperObs$V.mag[dipless],
		xlim=x.limits,ylim=y.limits,xlab=x.label,ylab=y.label,
		main=plot.title,pch= plot.pch,col=plot.col)
	points(x=desmat[!dipless],
			y=allSuperObs$V.mag[!dipless],
			col="darkgrey",pch=plot.pch)
	# plot the fit
	lines(x=desmat,y=fit.points,col="darkgrey",lwd=2)
############### unbinned #################
} else {
	desmat <- asassn_data$HJD[good_data] - tmin
	#calculate weights for data outside the dips
	my.weights <- (1/asassn_data$mag_err[good_data])*(as.numeric(not.dips))

	if(mars.fit) {
		#mars <- earth(x=desmat,y= asassn_data$mag[good_data],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty,subset=not.dips)
		mars <- earth(x=desmat,y= asassn_data$mag[good_data],
						nk= marsOrder,
						pmethod= marsPMethod,
						penalty = marsPenalty,
						weights=my.weights,
						thresh= mars.thresh,
						minspan=10)
		
		print(summary(mars))
	} else {
		smoove.fit <- smooth.spline(x=desmat,
									y= asassn_data$mag[good_data],
									w=my.weights,
									all.knots=FALSE,nknots= asassn.n.knots,
									keep.data=TRUE,
									cv=TRUE,
									penalty= 1)
		cat("\n\n Smooth Spline Fit: \n")
		print(smoove.fit$call)
		cat("\n")
		print(smoove.fit$fit)
	}
	# plot limits
	y.lims <- c(max(asassn_data$mag[good_data]) + mag.margin,min(asassn_data$mag[good_data]) - mag.margin)
	my.xlab <- paste("Julian Date - ",tmin)
	x.lims <- c(0,max(desmat) + tmargin)
	# plot
	quartz("asas-sn")
	plot.times <- asassn_data$HJD[good_data] - tmin
	# plot the points used in the fit
	plot(plot.times[not.dips],asassn_data$mag[good_data][not.dips],
		col="darkgreen",pch=20,
		main="ASASSN V Data",
		xlab=my.xlab,ylab="V magnitude",
		ylim=y.limits,xlim= x.limits)
	# add the points not used in the fit in grey
	points(plot.times[!not.dips],asassn_data$mag[good_data][!not.dips],col="grey",pch=20)
	# add grid
	
	if(mars.fit) {
		# plot MARS fit as a line
		#lines(x=desmat[not.dips],y=mars$fitted.values,col= "black",lwd=2)
		lines(x=desmat,y=mars$fitted.values,col= "black",lwd=2)
	} else {
		smoove.values <- predict(smoove.fit,desmat)$y
		lines(x=desmat,y=smoove.values,col="black",lwd=2)
	}

}

grid(col="black")

# plot vertical lines if any
jdLine <- jdLine - tmin
verticalDateLines(jdLine, jdLineText, y.limits, jdLineColor)
