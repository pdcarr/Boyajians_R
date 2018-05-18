#### housekeeping
library("MASS") # for rlm() and lqs()
library("earth")
library("Hmisc")
options(digits=12) # hard to read JDs without this setting
rm("superObs")
rm("allSuperObs")
source("data_funcs.R")
source("input_files/VlineParams.R")
source("plot_funcs.R")
source("special_funcs.R")
source("input_files/dip_mask.R")


######################## control parameters
maxAirmass <- 2.0 # data with airmass higher than this will not be included
bin.width = 60/1440 # days
#bin.width = 1440/1440 # days
#################################
t.epsilon = 0.0 # days
mag.epsilon <- 0.01 # magnitudes
bin.it <- TRUE
plot.col = "darkgreen"
plot.pch <- 20
#dfile.name <- "data/BruceGary.csv"
#dfile.name <- "data/BG_gprime_latest.txt"
##### bands, colors, files, knots
allBands <- data.frame(data.type=c("g'","r'","i'"),plotColor=c("blue","red","darkviolet"), 
				dfile.name = c("data/BG_gprime_latest.txt","data/BG_rprime_latest.csv","data/BG_iprime_latest.csv"),
				bg.n.knots = c(6,4,4),stringsAsFactors=FALSE)
				
n.pts.tail <- 12 # number of points in tail of super observations to print
#data.type <- "G prime"
#data.type <- "V"
source("input_files/dip_mask.R")
#earliest.MJD <- 58197
earliest.MJD <- NA
earliest.plot.MJD <- 58170
#latest.MJD <- 58098
latest.MJD <- NA
pretty.days <- 10
#pretty.days <- 1
########## MARS control
n.knots <- 33
knot.penalty <- 0
min.span <- 2
earth.thresh <- 0.00001
######## fitting parameters
plot.spline <- TRUE
#bg.n.knots <- 6
#######
plot.SG <- FALSE
###### read in the data
band.index <- 1
for (data.type in allBands$data.type) {
	my.cols <- c("numeric","numeric","numeric","numeric")
	bg.data <- read.csv(allBands$dfile.name[band.index],header=TRUE, stringsAsFactors=FALSE,colClasses=my.cols)
	
	# filter out the high airmass data
	ok.airmass <- bg.data$air.mass <= maxAirmass
	if(is.numeric(earliest.MJD)) {
		ok.airmass <- ok.airmass & (bg.data$MJD >= earliest.MJD)
	}
	if(is.numeric(latest.MJD)) {
		ok.airmass <- ok.airmass & (bg.data$MJD <= latest.MJD)
	}
	
	# find earliest and latest times and highest and lowest magnitudes
	
	# get the bin boundaries to be somewhere sensible
	
	tmin <-  min(bg.data$MJD[ok.airmass],na.rm=TRUE)
	tmin <- floor(tmin/pretty.days)*pretty.days
	
	tmax = max(bg.data$MJD[ok.airmass],na.rm=TRUE)
	x.limits <- c(0,(tmax-tmin) + t.epsilon)
	y.limits <- c(max(bg.data$V.mag[ok.airmass],na.rm=TRUE) - mag.epsilon,
				 min(bg.data$V.mag[ok.airmass],na.rm=TRUE) + mag.epsilon)
	x.label <- paste("MJD - ",tmin)
	
	# if binning, do this:
	if (bin.it) {
		plot.title <- paste("Data binned with width",bin.width,"Days")
		y.label <- paste("Mean",data.type,"magnitude over bin")
		allSuperObs <- data.frame(MJD=numeric(),V.mag=numeric(),Uncertainty=numeric(),stringsAsFactors=FALSE)
		superObs <- data.frame(MJD=numeric(),V.mag=numeric(),Uncertainty=numeric(),stringsAsFactors=FALSE)
		
	
	
		for (t in seq(tmin,tmax,by=bin.width)) {
			in.bin <- !is.na(bg.data$MJD) & bg.data$MJD >= t & bg.data$MJD < (t + bin.width) & ok.airmass
			n.in.bin <- length(bg.data$MJD[in.bin])
			if (n.in.bin > 0) {
				superObs[1,"MJD"] <- mean(bg.data$MJD[in.bin],na.rm=TRUE)
				superObs[1,"V.mag"] <- mean(bg.data$V.mag[in.bin],na.rm=TRUE)
				superObs[1,"Uncertainty"] <- mean(bg.data$SE[in.bin],na.rm=TRUE)/sqrt(n.in.bin)
				allSuperObs <- rbind(allSuperObs,superObs)
			} else {
				next
			}
		}
		# mask dips from the fit
		bg.binned.JDs <- allSuperObs$MJD  + 2400000.5
	    dipless <- filter.dips.JD(bg.binned.JDs,dip.mask)
	    binWeights <- as.numeric(dipless) # weight of 1 if not in a known dip, 0 otherwise.
		# calculate robust linear fit
		desmat <- allSuperObs$MJD - tmin
	    #theFit <- rlm(allSuperObs$V.mag ~ desmat,na.action="na.omit",psi=psi.bisquare,subset=dipless)
	#    theFit <- earth(x=desmat,
	#                    y=allSuperObs$V.mag,
	#                    weights=binWeights,
	#                    pmethod="exhaustive",
	 #                   nk=n.knots,
	#                    penalty=knot.penalty,
	#                    minspan=min.span,
	#                    thresh=earth.thresh)
	
		# open a window to plot in
		q.string <- paste("binned data plot in ", data.type)
		quartz(q.string)
		# error bars
		my.y.plus <- allSuperObs$V.mag[dipless]  + allSuperObs$Uncertainty[dipless]
		my.y.minus <- allSuperObs$V.mag[dipless] -  allSuperObs$Uncertainty[dipless]
		# plot bins used in the fit
		errbar(desmat[dipless],allSuperObs$V.mag[dipless],
			yplus=my.y.plus,yminus=my.y.minus,
			xlim=x.limits,ylim=y.limits,
			xlab=x.label,ylab=y.label,
			main=plot.title,
			pch= plot.pch,col=allBands$plotColor[band.index])
			
			title(plot.title)
		
		# plot omitted bins
		points(x=desmat[!dipless],
			y=allSuperObs$V.mag[!dipless],
			pch= plot.pch,col="grey")
	# not binned
	} else {
		quartz("Data Plot")
		y.label = "V Magnitude"
		plot.title <- paste("Bruce Gary Data with Airmass <=",maxAirmass)
		desmat <- bg.data$MJD[ok.airmass] - tmin
		theFit <- rlm(bg.data$V.mag[ok.airmass] ~ desmat,na.action="na.omit",psi=psi.bisquare)
		plot(desmat,bg.data$V.mag[ok.airmass],
			xlim=x.limits,
			ylim=y.limits,
			xlab=x.label,
			ylab=y.label,
			main=plot.title,
			pch= plot.pch,
			col=allBands$plotColor[band.index])
		grid(col="black")
	}
	
	#### plot the fit if selected
	if(plot.spline) {
		theFit <- smooth.spline(x=desmat,
			y= allSuperObs$V.mag,
			w= binWeights,
			all.knots=FALSE,nknots= allBands$bg.n.knots[band.index],
			keep.data=TRUE,
			cv=TRUE,
			penalty= 1)
		# let's see what the fit is
		cat("\n\n Smooth Spline Fit: \n")
		print(theFit$call)
		cat("\n")
		print(theFit$fit)
		fit.points <- predict(theFit,x=desmat)$y
	
		#myslope <- coefficients(theFit)[2]
		#basemag <- coefficients(theFit)[1]
		#curve(basemag + myslope*x,from=0,to=tmax-tmin, add=TRUE,col="black")
		lines(x=desmat,y=fit.points,col="black")
	
	}
	##### draw vertical lines at various dates
	jdLine <- jdLine - 2400000.5 - tmin
	if(drawDateLine) { verticalDateLines(jdLine, jdLineText, y.limits, jdLineColor)}
	
	# add a grid
	grid(col="black")
	allSuperObs <- cbind(allSuperObs,"spline fit" = fit.points)
	print(tail(allSuperObs,n=12))
	summary(theFit)

	band.index <- band.index + 1
}

if(plot.SG) {
	source("SG_plot_AAVSO.R")
}
